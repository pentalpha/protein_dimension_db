import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
import sys
import os
import pickle
import json
from collections import Counter

# ==========================================
# 1. THE MODEL ARCHITECTURE
# ==========================================


class InterProAutoencoder(nn.Module):
    def __init__(self, input_dim, embedding_dim):
        super(InterProAutoencoder, self).__init__()
        # 20k -> 4096 -> 2000
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 4096),
            nn.ReLU(),
            nn.Linear(4096, embedding_dim),
        )
        # 2000 -> 4096 -> 20k
        self.decoder = nn.Sequential(
            nn.Linear(embedding_dim, 4096), nn.ReLU(), nn.Linear(4096, input_dim)
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
        self.data = family_lists
        self.vocab_map = vocab_map
        self.input_dim = input_dim

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        vector = torch.zeros(self.input_dim, dtype=torch.float32)
        for f in self.data[idx]:
            if f in self.vocab_map:
                vector[self.vocab_map[f]] = 1.0
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
        self.history = []

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
        max_epochs_no_improve=4,
        directory=None,
    ):
        """Trains the model on a list of lists of InterPro families."""
        self._build_vocab(family_lists)

        self.model = InterProAutoencoder(self.actual_input_dim, self.embedding_dim).to(
            self.device
        )
        dataset = InterProDataset(family_lists, self.vocab_map, self.actual_input_dim)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(self.model.parameters(), lr=lr)

        self.model.train()

        best_loss = float("inf")
        best_loss_epoch = 0

        for epoch in range(epochs):
            total_loss = 0
            for batch in dataloader:
                batch = batch.to(self.device)
                reconstruction, _ = self.model(batch)
                loss = criterion(reconstruction, batch)

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                total_loss += loss.item()

            print(f"Epoch {epoch+1}/{epochs} - Loss: {total_loss/len(dataloader):.6f}")
            loss_rounded = round(total_loss, 6)
            self.history.append(
                {"epoch": epoch, "loss": loss_rounded, "best_loss": best_loss}
            )
            if loss_rounded < best_loss:
                best_loss = loss_rounded
                best_loss_epoch = epoch

                # save model
                if directory:
                    if os.path.exists(directory):
                        self.save(directory)

            elif epoch - best_loss_epoch > max_epochs_no_improve:
                print(f"Early stopping at epoch {epoch+1}")
                print(
                    f"Last improvement: epoch {best_loss_epoch} with loss {best_loss}"
                )
                break

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
        instance.model.load_state_dict(
            torch.load(
                os.path.join(directory, "model.pth"), map_location=instance.device
            )
        )
        instance.model.to(instance.device)
        instance.model.eval()

        # Load history
        with open(os.path.join(directory, "history.json"), "r") as f:
            instance.history = json.load(f)
        return instance
