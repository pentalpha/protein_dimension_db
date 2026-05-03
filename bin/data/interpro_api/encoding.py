import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import os

from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os
import pickle
import json
from collections import Counter
import multiprocessing

import torch.multiprocessing as mp

# Force PyTorch to use spawn instead of fork to prevent UCX segfaults
try:
    mp.set_start_method("spawn", force=True)
except RuntimeError:
    pass

# ==========================================
# 1. THE MODEL ARCHITECTURE
# ==========================================

home_dir = os.path.expanduser("~")
cache_dir = os.path.join(home_dir, ".cache")
os.environ["TRITON_CACHE_DIR"] = cache_dir + "/triton"
os.environ["TORCH_HOME"] = cache_dir + "/torch"

intermediary_len = 4096


class InterProAutoencoder(nn.Module):
    def __init__(self, input_dim, embedding_dim, dropout_rate=0.2):
        super(InterProAutoencoder, self).__init__()
        # 16k -> intermediary_len -> 1200-256
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, intermediary_len),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(intermediary_len, embedding_dim),
        )
        # 16k -> intermediary_len -> 1200-256
        self.decoder = nn.Sequential(
            nn.Linear(embedding_dim, intermediary_len),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(intermediary_len, input_dim),
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded, encoded


# ==========================================
# 2. DATASET & WRAPPER CLASS
# ==========================================


class InterProDataset(Dataset):
    def __init__(self, family_lists, vocab_map, input_dim):
        self.input_dim = input_dim
        # Pre-compute integer indices once to avoid dictionary lookups in the training loop
        print("Pre-mapping dataset to integer indices...")
        self.data_indices = [
            [vocab_map[f] for f in fam_list if f in vocab_map]
            for fam_list in family_lists
        ]

    def __len__(self):
        return len(self.data_indices)

    def __getitem__(self, idx):
        # Use fast tensor indexing instead of a Python for-loop
        vector = torch.zeros(self.input_dim, dtype=torch.float32)
        indices = self.data_indices[idx]
        if indices:
            vector[indices] = 1.0
        return vector


class AutoEncoderWrapper:
    def __init__(self, input_dim=20000, embedding_dim=2000, predefined_vocab=None):
        self.input_dim = input_dim
        self.embedding_dim = embedding_dim
        self.vocab_map = None
        self.model = None
        if predefined_vocab is not None:
            if type(predefined_vocab) == list:
                self.vocab_map = {token: i for i, token in enumerate(predefined_vocab)}
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.using_cpu = not torch.cuda.is_available()
        print(f"Using device: {self.device}")
        self.history = []
        self.actual_input_dim = None

        self.loader_procs = multiprocessing.cpu_count() // 3
        self.processing_procs = round(multiprocessing.cpu_count() * 0.6)
        torch.set_num_threads(self.processing_procs)
        print(
            f"Using {self.loader_procs} loader cores and {self.processing_procs} processing cores"
        )

    def _build_vocab(self, family_lists):
        """Creates a mapping for the top N most frequent InterPro families."""
        if self.vocab_map is None:
            print(f"Building vocabulary for top {self.input_dim} families...")
            all_fams = [f for sublist in family_lists for f in sublist]
            counter = Counter(all_fams)
            print(f"Found {len(counter)} unique families.")
            most_common = counter.most_common(self.input_dim)
            print(f"Using top {self.input_dim} families.")
            print(f"Top 5 families: {most_common[:5]}")
            print(f"Bottom 5 families: {most_common[-5:]}")
            self.vocab_map = {token: i for i, (token, count) in enumerate(most_common)}
        fams_set = set(self.vocab_map.keys())
        # calculate percentage of proteins with at least one family in the vocabulary
        fams_count = 0
        for fam_list in family_lists:
            if len(set(fam_list).intersection(fams_set)) > 0:
                fams_count += 1
        perc = fams_count / len(family_lists) * 100
        print(
            f"Percentage of proteins with at least one family in the vocabulary: {perc:.2f}%"
        )
        # Update input_dim if data is smaller than 20k
        self.actual_input_dim = len(self.vocab_map)

    def fit(
        self,
        family_lists,
        epochs=10,
        batch_size=64,
        lr=8e-4,
        max_epochs_no_improve=5,
        directory=None,
        optimize_cpu=False,
        eval_perc=0.05,
    ):
        pos_weight_value = 250.0
        # torch.set_num_threads(self.processing_procs)
        if optimize_cpu:
            import os

            # Force PyTorch to use the system C++ compiler instead of the broken Conda one
            os.environ["CXX"] = "/usr/bin/g++"
            os.environ["CC"] = "/usr/bin/gcc"
        """Trains the model on a list of lists of InterPro families."""
        print(
            f"Training params: epochs={epochs}, batch_size={batch_size}, "
            f"lr={lr}, max_epochs_no_improve={max_epochs_no_improve}, "
            f"optimize_cpu={optimize_cpu}, pos_weight={pos_weight_value}"
        )
        bar = tqdm(total=epochs + 5)

        print("Building Vocab")
        self._build_vocab(family_lists)
        bar.update(1)

        print("Initializing Model", file=sys.stderr)
        self.model = InterProAutoencoder(
            self.actual_input_dim, self.embedding_dim, dropout_rate=0.2
        ).to(self.device)
        # Tell PyTorch to optimize the execution graph for the CPU
        """if optimize_cpu:
            try:
                compiled = torch.compile(self.model)
                self.model = compiled
            except Exception as e:
                print(f"Could not compile model: {e}", file=sys.stderr)"""
        bar.update(1)

        if eval_perc > 0:
            n_eval = int(len(family_lists) * eval_perc)
            eval_indexes = np.random.choice(len(family_lists), n_eval, replace=False)
            eval_data = [family_lists[i] for i in eval_indexes]

        else:

            eval_data = None
        train_data = family_lists

        print("Creating Dataset", file=sys.stderr)
        dataset = InterProDataset(train_data, self.vocab_map, self.actual_input_dim)

        if eval_data is not None:
            eval_dataset = InterProDataset(
                eval_data, self.vocab_map, self.actual_input_dim
            )
        else:
            eval_dataset = None
        bar.update(1)

        print("Creating DataLoader", file=sys.stderr)
        if optimize_cpu:
            dataloader = DataLoader(
                dataset,
                batch_size=batch_size,
                shuffle=True,
                num_workers=0,  # Avoid seg fault
                pin_memory=False,  # Disable pinning for pure CPU execution
                # prefetch_factor=2,
            )
        else:
            dataloader = DataLoader(
                dataset, batch_size=batch_size, shuffle=True, pin_memory=True
            )

        if eval_dataset is not None:
            eval_dataloader = DataLoader(
                eval_dataset,
                batch_size=1600,
                shuffle=False,
                num_workers=0,
                pin_memory=False,
            )
        else:
            eval_dataloader = None

        bar.update(1)

        print("Initializing training", file=sys.stderr)
        pos_weight = torch.tensor([pos_weight_value]).to(self.device)

        criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        optimizer = optim.Adam(self.model.parameters(), lr=lr)

        self.model.train()
        bar.update(1)

        best_loss = float("inf")
        best_loss_epoch = 0

        n_batchs = len(dataset) // batch_size

        print("Training Autoencoder", file=sys.stderr)
        for epoch in range(epochs):
            total_loss = 0
            batch_n = 0
            sub_bar = tqdm(total=n_batchs)
            for batch in dataloader:
                batch = batch.to(self.device)
                # print(f"Epoch {epoch} - Batch {batch_n} - Allocated", file=sys.stderr)
                reconstruction, _ = self.model(batch)
                # print(
                #    f"Epoch {epoch} - Batch {batch_n} - Runned model", file=sys.stderr
                # )
                loss = criterion(reconstruction, batch)
                # print(f"Epoch {epoch} - Batch {batch_n} - Criterion", file=sys.stderr)

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                # print(f"Epoch {epoch} - Batch {batch_n} - Step", file=sys.stderr)
                total_loss += loss.item()
                batch_n += 1
                sub_bar.update(1)
            sub_bar.close()
            # Flush the lingering training tensors from VRAM
            if self.device.type == "cuda":
                torch.cuda.empty_cache()

            if eval_dataloader is not None:
                """The dataset lines are binary vectors, so we can calculate
                the precision and recall of the reconstruction. We dont need eval_loss
                """

                self.model.eval()
                with torch.no_grad():
                    num_thresholds = 30
                    thresholds = torch.linspace(
                        0.05, 0.95, steps=num_thresholds, device=self.device
                    )

                    # Accumulators for each threshold
                    sum_precision = torch.zeros(num_thresholds, device=self.device)
                    sum_recall = torch.zeros(num_thresholds, device=self.device)

                    epsilon = 1e-8

                    for batch in eval_dataloader:
                        batch = batch.to(self.device)

                        # 1. Get raw logits from the model
                        raw_reconstruction, _ = self.model(batch)

                        probs = torch.sigmoid(raw_reconstruction)

                        actual_pos = batch.sum(dim=1)

                        # Calculate metrics for all 30 thresholds
                        for i, thresh in enumerate(thresholds):
                            preds = (probs > thresh).float()

                            tp = (preds * batch).sum(dim=1)
                            pred_pos = preds.sum(dim=1)

                            # Add the batch mean to our global accumulators
                            sum_precision[i] += (tp / (pred_pos + epsilon)).mean()
                            sum_recall[i] += (tp / (actual_pos + epsilon)).mean()

                        # Aggressively clear VRAM for the next batch
                        del (
                            batch,
                            raw_reconstruction,
                            probs,
                            preds,
                            tp,
                            pred_pos,
                            actual_pos,
                        )

                    # 1. Average the accumulated metrics across all batches
                    avg_precision = sum_precision / len(eval_dataloader)
                    avg_recall = sum_recall / len(eval_dataloader)

                    # 2. Calculate F1 Score for all thresholds simultaneously
                    # F1 = 2 * (Precision * Recall) / (Precision + Recall)
                    f1_scores = (2 * avg_precision * avg_recall) / (
                        avg_precision + avg_recall + epsilon
                    )

                    # 3. Find the threshold that produced the absolute best F1 Score
                    best_idx = torch.argmax(f1_scores)

                    best_f1 = float(f1_scores[best_idx].item())
                    best_thresh = thresholds[best_idx].item()
                    best_prec = float(avg_precision[best_idx].item())
                    best_rec = float(avg_recall[best_idx].item())

                    print(
                        f"\nEval Fmax: {best_f1:.4f} (at threshold {best_thresh:.2f}) | "
                        f"Precision: {best_prec:.4f} | Recall: {best_rec:.4f}"
                    )
            else:
                print(
                    f"Epoch {epoch+1}/{epochs} - Loss: {total_loss/len(dataloader):.6f}"
                )
                best_prec = None
                best_rec = None
                best_f1 = None
            loss_rounded = round(total_loss, 6)
            self.history.append(
                {
                    "epoch": epoch,
                    "loss": loss_rounded,
                    "best_loss": best_loss,
                    "precision": best_prec,
                    "recall": best_rec,
                }
            )
            if loss_rounded < best_loss:
                best_loss = loss_rounded
                best_loss_epoch = epoch

                # save model
                if directory:
                    import os

                    if os.path.exists(directory):
                        self.save(directory)

            elif epoch - best_loss_epoch > max_epochs_no_improve:
                print(f"Early stopping at epoch {epoch+1}")
                print(
                    f"Last improvement: epoch {best_loss_epoch} with loss {best_loss}"
                )
                break
            bar.update(1)

        bar.close()

    def predict(self, family_lists):
        """Encodes a list of lists into embedding vectors (Latent Space)."""
        if self.model is None or self.vocab_map is None:
            raise ValueError("Model not trained or loaded.")

        self.model.eval()
        dataset = InterProDataset(family_lists, self.vocab_map, self.actual_input_dim)
        dataloader = DataLoader(dataset, batch_size=128, shuffle=False)

        embeddings = []
        with torch.no_grad():
            for batch in dataloader:
                batch = batch.to(self.device)
                _, encoded = self.model(batch)
                embeddings.append(encoded.cpu().numpy())

        return np.vstack(embeddings)

    def save(self, directory):
        """Saves both the model state and the vocabulary mapping."""
        os.makedirs(directory, exist_ok=True)
        # Save PyTorch Model
        torch.save(self.model.state_dict(), os.path.join(directory, "model.pth"))
        # Save Metadata (Vocab and Config)
        meta = {
            "vocab_map": self.vocab_map,
            "input_dim": self.input_dim,
            "actual_input_dim": self.actual_input_dim,
            "embedding_dim": self.embedding_dim,
        }
        with open(os.path.join(directory, "metadata.pkl"), "wb") as f:
            pickle.dump(meta, f)
        print(f"Model and metadata saved to {directory}")

        # Save history as json
        with open(os.path.join(directory, "history.json"), "w") as f:
            json.dump(self.history, f, indent=4)
        print(f"History saved to {directory}")

    @classmethod
    def load(cls, directory):
        """Loads a saved instance."""
        with open(os.path.join(directory, "metadata.pkl"), "rb") as f:
            meta = pickle.load(f)

        instance = cls(input_dim=meta["input_dim"], embedding_dim=meta["embedding_dim"])
        instance.vocab_map = meta["vocab_map"]
        instance.actual_input_dim = meta["actual_input_dim"]

        instance.model = InterProAutoencoder(
            instance.actual_input_dim, instance.embedding_dim
        )
        # Load the raw state dictionary from disk
        raw_state_dict = torch.load(
            os.path.join(directory, "model.pth"), map_location=instance.device
        )

        # Strip out the '_orig_mod.' prefix added by torch.compile()
        clean_state_dict = {
            key.replace("_orig_mod.", ""): value
            for key, value in raw_state_dict.items()
        }

        # Load the cleaned dictionary into the standard model
        instance.model.load_state_dict(clean_state_dict)

        instance.model.to(instance.device)
        instance.model.eval()

        # Load history
        with open(os.path.join(directory, "history.json"), "r") as f:
            instance.history = json.load(f)
        return instance
