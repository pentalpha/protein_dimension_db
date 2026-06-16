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
import scipy

import torch.multiprocessing as mp

import matplotlib.pyplot as plt
from collections import Counter
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import glob
from PIL import Image


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

intermediary_len = 3600


def generate_pca_gif(directory, output_filename="latent_evolution.gif", duration=400):
    """
    Compiles all saved PCA scatter plots into an animated GIF.
    duration: milliseconds per frame (400ms = 2.5 frames per second)
    """
    if directory is None:
        return

    print(f"\nStitching PCA frames into GIF...")

    # Find all pca images and sort them alphabetically (which matches chronologically due to the 000 padding)
    search_pattern = os.path.join(directory, "pca_epoch_*.png")
    image_paths = sorted(glob.glob(search_pattern))

    if not image_paths:
        print("No PCA images found to make a GIF.")
        return

    # Open all images
    frames = [Image.open(img_path) for img_path in image_paths]

    # Save as GIF
    output_path = os.path.join(directory, output_filename)

    # Save the first frame, and append the rest
    frames[0].save(
        output_path,
        format="GIF",
        append_images=frames[1:],
        save_all=True,
        duration=duration,
        loop=0,  # 0 means infinite loop
    )
    print(f"GIF successfully updated at: {output_path}")


def prepare_pca_clusters(eval_data, top_n=14):
    """
    Finds the top N terms in the evaluation set and assigns proteins
    to strictly mutually exclusive clusters for clean visualization.
    """
    print(f"Preparing PCA clusters for top {top_n} terms...")

    # 1. Find the top N most common terms in the eval set
    all_eval_terms = [term for protein in eval_data for term in protein]
    term_counts = Counter(all_eval_terms)
    top_terms = [term for term, _ in term_counts.most_common(top_n)]
    top_terms_set = set(top_terms)

    target_indices = []
    cluster_labels = []

    # 2. Filter proteins: must contain EXACTLY ONE of the top N terms
    for i, protein in enumerate(eval_data):
        if len(protein) >= 3:
            overlap = top_terms_set.intersection(protein)
            if len(overlap) == 1:
                label = list(overlap)[0]
                target_indices.append(i)
                cluster_labels.append(label)

    print(f"Selected {len(target_indices)} purely clustered proteins for PCA.")
    different_labels = len(set(cluster_labels))
    print(f"Different labels: {different_labels}")
    return target_indices, cluster_labels, top_terms


def plot_epoch_pca(
    eval_embeddings, target_indices, cluster_labels, top_terms, epoch, directory
):
    """
    Calculates 2D PCA for the selected embeddings and saves plot.
    """
    # 1. Subset the embeddings to only our mutually exclusive proteins
    X = eval_embeddings[target_indices]

    # 2. Compute 2D PCA
    # (PCA on a few thousand rows is virtually instant on CPU)
    # pca = PCA(n_components=2)
    # X_pca = pca.fit_transform(X)

    # Use t-SNE for better visualization of local structure
    tsne = TSNE(
        n_components=2, perplexity=30, learning_rate="auto", init="pca", random_state=42
    )
    X_pca = tsne.fit_transform(X)

    fig, ax = plt.subplots(figsize=(4.2, 4.2), dpi=180)

    # Assign distinct colors using the tab20 colormap
    cmap = plt.get_cmap("tab20")
    colors = {term: cmap(i) for i, term in enumerate(top_terms)}

    # Scatter plot by cluster
    for term in top_terms:
        # Create a boolean mask for the current term
        term_mask = [label == term for label in cluster_labels]

        if any(term_mask):
            ax.scatter(
                X_pca[term_mask, 0],
                X_pca[term_mask, 1],
                label=term,
                color=colors[term],
                alpha=0.7,
                s=8,  # Small point size to prevent overlapping blobs
                edgecolors="none",
            )

    # 4. Formatting
    ax.set_title(
        f"Compact Representation of Protein Families - {epoch}",
        fontsize=10,
    )
    # Remove axis ticks for a cleaner look
    ax.set_xticks([])
    ax.set_yticks([])

    # Place a tiny legend outside the main plot area so it doesn't cover data
    ax.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)

    # 5. Save the plot
    if directory is None:
        directory = "."
    os.makedirs(directory, exist_ok=True)

    # bbox_inches='tight' ensures the legend isn't cut off
    filepath = os.path.join(directory, f"pca_epoch_{epoch:03d}.png")
    plt.savefig(filepath, format="png", bbox_inches="tight")
    plt.close(fig)  # Free memory


def centroid_correlation_metric(eval_targets, eval_embeddings, term_sim_matrix):
    # 1. Calculate the sum of embeddings for each term
    # Matrix Mult: [num_terms, N_eval] @ [N_eval, embedding_dim] -> [num_terms, embedding_dim]
    term_sums = torch.mm(eval_targets.t(), eval_embeddings)

    # 2. Count how many proteins contain each term
    term_counts = eval_targets.sum(dim=0).unsqueeze(1)  # Shape: [num_terms, 1]

    # 3. Filter out terms that didn't appear in the evaluation set
    valid_terms_mask = term_counts.squeeze() > 0

    valid_term_sums = term_sums[valid_terms_mask]
    valid_term_counts = term_counts[valid_terms_mask]

    # 4. Divide by counts to get the true mathematical Centroid (mean)
    term_centroids = (
        valid_term_sums / valid_term_counts
    )  # Shape: [Valid_Terms, embedding_dim]

    # 5. Calculate pairwise Euclidean distances between the centroids
    centroid_distances = torch.cdist(term_centroids, term_centroids, p=2).numpy()

    # Subset the ground truth matrix to only the terms present in this eval batch
    valid_indices = torch.where(valid_terms_mask)[0].cpu().numpy()
    valid_term_sims = term_sim_matrix[np.ix_(valid_indices, valid_indices)]

    # --- Spearman Correlation ---
    n_valid = len(valid_indices)
    triu_indices = np.triu_indices(n_valid, k=1)

    flat_centroid_dist = centroid_distances[triu_indices]
    flat_term_sim = valid_term_sims[triu_indices]

    spearman_corr, p_value = scipy.stats.spearmanr(flat_centroid_dist, flat_term_sim)

    # Since Distance and Similarity are inverse, a perfect score is -1.0
    spearman_corrected = abs(spearman_corr)
    return float(spearman_corrected), float(p_value)


def embedding_neighborhood_score(eval_targets, eval_embeddings, k=5):
    """
    Computes the average fraction of k-nearest neighbours that share
    at least one term with the query protein. Higher is better.

    Args:
        eval_targets: torch.Tensor of shape (N, num_terms), binary.
        eval_embeddings: torch.Tensor of shape (N, dim).
        k: number of neighbours to consider.

    Returns:
        float: average neighbour label agreement.
    """
    # Normalise embeddings for cosine similarity (optional but often helps)
    normed_emb = eval_embeddings / (eval_embeddings.norm(dim=1, keepdim=True) + 1e-8)

    # Compute pairwise cosine similarity matrix
    sim = torch.mm(normed_emb, normed_emb.t())  # (N, N)

    # For each protein, we want its k nearest neighbours excluding itself.
    # Set diagonal to -inf so the query itself isn't chosen
    sim.fill_diagonal_(-float("inf"))

    # Get top‑k indices for each row
    _, knn_indices = torch.topk(sim, k, dim=1)  # (N, k)

    # For each query, check which neighbours share at least one term.
    # targets[i] has shape (num_terms,). We need to check if
    # (targets[i] * targets[neighbour]) sum > 0 for each neighbour.
    N = eval_targets.size(0)
    agreement_fracs = torch.zeros(N, device=eval_targets.device)

    for i in range(N):
        query = eval_targets[i].unsqueeze(0)  # (1, num_terms)
        neighbours = eval_targets[knn_indices[i]]  # (k, num_terms)
        # overlap[i] = 1 if sum(query * neighbour) > 0 else 0
        overlaps = (torch.sum(query * neighbours, dim=1) > 0).float()
        agreement_fracs[i] = overlaps.mean()

    return agreement_fracs.mean().item()


class InterProAutoencoder(nn.Module):
    def __init__(self, input_dim, embedding_dim, dropout_rate=0.1):
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
        self.predefined_vocab = predefined_vocab
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

    def _calc_ont_sim_matrix(self, ann_sim_matrix):
        import numpy as np

        zero_sim_proteins_perc = 0.25
        n = len(ann_sim_matrix)

        # 1. Get indices for the strictly lower triangle (i > j)
        row_idx, col_idx = np.tril_indices(n, k=-1)

        # 2. Extract values for these indices
        lower_tri_vals = ann_sim_matrix[row_idx, col_idx]

        # 3. Create boolean masks for > 0.01 and <= 0.01
        min_val = 0.15
        non_zero_mask = lower_tri_vals > min_val
        zero_mask = lower_tri_vals <= min_val

        # 4. Filter for non-zero pairs and their values
        non_zero_rows = row_idx[non_zero_mask]
        non_zero_cols = col_idx[non_zero_mask]
        non_zero_vals = lower_tri_vals[non_zero_mask]

        n_non_zero = len(non_zero_vals)
        print(f"Number of non-zero pairs: {n_non_zero}")

        # 5. Determine how many zeros to sample
        n_zero_pairs_to_sample = int(n_non_zero * zero_sim_proteins_perc)
        print(f"Number of zero pairs to sample: {n_zero_pairs_to_sample}")

        # Filter out the zero pairs
        zero_rows = row_idx[zero_mask]
        zero_cols = col_idx[zero_mask]
        zero_vals = lower_tri_vals[zero_mask]

        # Safety check: ensure we don't try to sample more zeros than exist
        n_zero_pairs_to_sample = min(n_zero_pairs_to_sample, len(zero_vals))

        # 6. Randomly sample the zero pairs efficiently
        sample_indices = np.random.choice(
            len(zero_vals), size=n_zero_pairs_to_sample, replace=False
        )

        sampled_zero_rows = zero_rows[sample_indices]
        sampled_zero_cols = zero_cols[sample_indices]
        sampled_zero_vals = zero_vals[sample_indices]

        print(f"Number of zero pairs sampled: {len(sample_indices)}")

        # 7. Combine the arrays
        combined_rows = np.concatenate([non_zero_rows, sampled_zero_rows])
        combined_cols = np.concatenate([non_zero_cols, sampled_zero_cols])
        pairs_for_sim_calculation = list(
            zip(combined_rows.tolist(), combined_cols.tolist())
        )

        # Final values array
        pair_sims = np.concatenate([non_zero_vals, sampled_zero_vals])

        return pairs_for_sim_calculation, pair_sims

    def fit(
        self,
        family_lists,
        ia_weights: dict,
        obo_path,
        epochs=10,
        batch_size=64,
        lr=8e-4,
        max_epochs_no_improve=5,
        directory=None,
        optimize_cpu=False,
        eval_perc=0.3333,
    ):
        assert eval_perc > 0
        pos_weight_value = 5.0
        dropout_rate = 0.2
        # torch.set_num_threads(self.processing_procs)
        if optimize_cpu:
            import os

            # Force PyTorch to use the system C++ compiler instead of the broken Conda one
            os.environ["CXX"] = "/usr/bin/g++"
            os.environ["CC"] = "/usr/bin/gcc"
        """Trains the model on a list of lists of InterPro families."""
        self.metaparams_archive = {
            "epochs": epochs,
            "batch_size": batch_size,
            "lr": lr,
            "max_epochs_no_improve": max_epochs_no_improve,
            "optimize_cpu": optimize_cpu,
            "pos_weight": pos_weight_value,
            "eval_perc": eval_perc,
            "intermediary_len": intermediary_len,
            "dropout_rate": dropout_rate,
        }
        print(f"Training params: {self.metaparams_archive}")

        bar = tqdm(total=epochs + 5)

        print("Building Vocab")
        self._build_vocab(family_lists)
        bar.update(1)

        print("Initializing Model", file=sys.stderr)
        self.model = InterProAutoencoder(
            self.actual_input_dim, self.embedding_dim, dropout_rate=dropout_rate
        ).to(self.device)
        bar.update(1)

        n_eval = int(len(family_lists) * eval_perc)
        eval_indexes_tmp_file = (
            f"/tmp/eval_interpro_autoencoder_{n_eval}_of_{len(family_lists)}.json"
        )
        ann_sim_matrix_file = f"/tmp/ann_sim_matrix_interpro_autoencoder_{n_eval}_of_{len(family_lists)}.npy"
        from bioinfo_utils.sim import OntologyTermSimilarity

        """sim_calculator = OntologyTermSimilarity(
            self.predefined_vocab, ia_weights, obo_path
        )"""
        import os

        if os.path.exists(eval_indexes_tmp_file):
            print(f"Using cached evaluation indexes from {eval_indexes_tmp_file}")
            eval_indexes = json.load(open(eval_indexes_tmp_file))["indexes"]
            eval_data = [family_lists[i] for i in eval_indexes]
            print(
                f"Loading cached ontology similarity matrix from {ann_sim_matrix_file}"
            )
            # ann_sim_matrix = np.load(ann_sim_matrix_file)
        else:
            print("Creating new evaluation indexes")
            eval_indexes = np.random.choice(len(family_lists), n_eval, replace=False)
            eval_data = [family_lists[i] for i in eval_indexes]
            json.dump(
                {"indexes": eval_indexes.tolist()},
                open(eval_indexes_tmp_file, "w"),
            )

            # ann_sim_matrix = sim_calculator.calculate_bma_matrix(eval_data)
            # save matrix
            # np.save(ann_sim_matrix_file, ann_sim_matrix)
        print(f"Loaded eval data: {len(eval_data)}")

        pca_indices, pca_labels, pca_terms = prepare_pca_clusters(eval_data, top_n=14)

        # pairs_for_sim_calculation, pair_sims = self._calc_ont_sim_matrix(ann_sim_matrix)

        train_data = family_lists

        print("Creating Dataset", file=sys.stderr)
        dataset = InterProDataset(train_data, self.vocab_map, self.actual_input_dim)

        eval_dataset = InterProDataset(eval_data, self.vocab_map, self.actual_input_dim)
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

        eval_dataloader = DataLoader(
            eval_dataset,
            batch_size=1600,
            shuffle=False,
            num_workers=0,
            pin_memory=False,
        )

        bar.update(1)

        print("Initializing training", file=sys.stderr)
        if pos_weight_value is not None:
            pos_weight = torch.tensor([pos_weight_value]).to(self.device)
            criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        else:
            criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(self.model.parameters(), lr=lr)

        self.model.train()
        bar.update(1)

        # best_loss = float("inf")
        best_score = float("-inf")
        best_loss_epoch = 0

        n_batchs = len(dataset) // batch_size

        print("Training Autoencoder", file=sys.stderr)
        for epoch in range(epochs):
            total_loss = 0
            batch_n = 0
            # sub_bar = tqdm(total=n_batchs)
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
                # sub_bar.update(1)
            # sub_bar.close()
            # Flush the lingering training tensors from VRAM
            if self.device.type == "cuda":
                torch.cuda.empty_cache()

            """The dataset lines are binary vectors, so we can calculate
            the precision and recall of the reconstruction. We dont need eval_loss
            """

            self.model.eval()
            with torch.no_grad():
                num_thresholds = 30
                thresholds = torch.linspace(
                    0.05, 0.95, steps=num_thresholds, device=self.device
                )

                # --- 1. Accumulators for Original Instance-centric metrics ---
                sum_precision = torch.zeros(num_thresholds, device=self.device)
                sum_recall = torch.zeros(num_thresholds, device=self.device)

                # --- 2. NEW: Accumulators for Term-centric metrics ---
                # We need TP and Pred_Pos for each threshold, for each term.
                num_terms = self.actual_input_dim
                term_tp = torch.zeros((num_thresholds, num_terms), device=self.device)
                term_pred_pos = torch.zeros(
                    (num_thresholds, num_terms), device=self.device
                )
                term_actual_pos = torch.zeros(num_terms, device=self.device)

                epsilon = 1e-8
                all_embeddings = []
                all_targets = []

                for batch in eval_dataloader:
                    batch = batch.to(self.device)

                    # Get raw logits and latent space
                    raw_reconstruction, encoded = self.model(batch)
                    all_embeddings.append(encoded.cpu())
                    all_targets.append(batch.cpu())

                    probs = torch.sigmoid(raw_reconstruction)

                    # Actual positives per protein (Instance)
                    actual_pos_inst = batch.sum(dim=1)
                    # NEW: Actual positives per term (Term)
                    term_actual_pos += batch.sum(dim=0)

                    # Calculate metrics for all 30 thresholds
                    for i, thresh in enumerate(thresholds):
                        preds = (probs > thresh).float()

                        # --- Instance Metrics (sum over terms: dim=1) ---
                        tp_inst = (preds * batch).sum(dim=1)
                        pred_pos_inst = preds.sum(dim=1)

                        sum_precision[i] += (tp_inst / (pred_pos_inst + epsilon)).mean()
                        sum_recall[i] += (tp_inst / (actual_pos_inst + epsilon)).mean()

                        # --- NEW: Term Metrics (sum over batch: dim=0) ---
                        term_tp[i] += (preds * batch).sum(dim=0)
                        term_pred_pos[i] += preds.sum(dim=0)

                    # Aggressively clear VRAM for the next batch
                    del (
                        batch,
                        raw_reconstruction,
                        probs,
                        preds,
                        tp_inst,
                        pred_pos_inst,
                        actual_pos_inst,
                    )

                # --- Original Instance Fmax & AUPRC ---
                avg_precision = sum_precision / len(eval_dataloader)
                avg_recall = sum_recall / len(eval_dataloader)

                f1_scores = (2 * avg_precision * avg_recall) / (
                    avg_precision + avg_recall + epsilon
                )
                best_idx = torch.argmax(f1_scores)
                best_f1 = float(f1_scores[best_idx].item())
                best_thresh = thresholds[best_idx].item()

                sorted_indices = torch.argsort(avg_recall)
                rec_sorted = avg_recall[sorted_indices]
                prec_sorted = avg_precision[sorted_indices]
                auprc = torch.trapz(prec_sorted, rec_sorted).item()

                # --- NEW: Term-Centric Fmax Calculation ---
                # Precision/Recall shape: (num_thresholds, num_terms)
                # We unsqueeze actual_pos to broadcast it to all thresholds
                term_precision = term_tp / (term_pred_pos + epsilon)
                term_recall = term_tp / (term_actual_pos.unsqueeze(0) + epsilon)

                term_f1s = (2 * term_precision * term_recall) / (
                    term_precision + term_recall + epsilon
                )

                # Find the best F1 score for EACH term across the 30 thresholds
                # best_term_f1s shape: (num_terms,)
                best_term_f1s, best_term_thresh_idx = torch.max(term_f1s, dim=0)

                # Optional: Retrieve the exact thresholds the model chose for each term
                # optimal_thresholds_per_term = thresholds[best_term_thresh_idx]

                # Average the best F1s to get a single global score.
                # CRITICAL: We only average over terms that actually appeared in the eval set,
                # otherwise terms with 0 true positives drag the score down unfairly.
                valid_terms_mask = term_actual_pos > 0
                mean_term_fmax = best_term_f1s[valid_terms_mask].mean().item()

                # --- Latent Distance & Spearman (Unchanged) ---
                eval_embeddings = torch.cat(all_embeddings, dim=0)
                eval_targets = torch.cat(all_targets, dim=0)

                if directory:  # Only plot if we have a directory to save to
                    plot_epoch_pca(
                        torch.sigmoid(eval_embeddings).numpy(),
                        pca_indices,
                        pca_labels,
                        pca_terms,
                        epoch,
                        directory,
                    )

                    if (epoch + 1) % 5 == 0 or epoch == (epochs - 1):
                        generate_pca_gif(directory)

                print(f"Calculating Term-Term Spearman Correlation...")

                """spearman_corr, p_value = centroid_correlation_metric(
                    eval_targets, eval_embeddings, sim_calculator.term_sim_matrix
                )"""
                knn_value_5 = embedding_neighborhood_score(
                    eval_targets, eval_embeddings, k=5
                )

                knn_value_3 = embedding_neighborhood_score(
                    eval_targets, eval_embeddings, k=3
                )

                knn_value_15 = embedding_neighborhood_score(
                    eval_targets, eval_embeddings, k=15
                )

                loss_rounded = round(total_loss, 6)

                # Update Epoch score to use your new Term-centric metric instead of Instance F1
                metric_ws = [2, 2, 1, 0.5, 0.5]
                metric_vals = [
                    (1.0 - loss_rounded),
                    mean_term_fmax,
                    knn_value_5,
                    knn_value_3,
                    knn_value_15,
                ]
                epoch_score = float(np.average(metric_vals, weights=metric_ws))
                score_rounded = round(epoch_score, 4)

                print(
                    f"\nInst Fmax: {best_f1:.4f} (th={best_thresh:.2f}) | "
                    f"Term Fmax: {mean_term_fmax:.4f} | "
                    f"AUPRC: {auprc:.4f} | "
                    f"Loss: {total_loss/len(dataloader):.6f} | "
                    f"kNN-sim5: {knn_value_5:.4f} | "
                    f"kNN-sim3: {knn_value_3:.4f} | "
                    f"kNN-sim15: {knn_value_15:.4f} | "
                    f"Epoch Score: {score_rounded:.4f}"
                )

            self.history.append(
                {
                    "epoch": epoch,
                    "loss": loss_rounded,
                    "best_score": best_score,
                    "auprc": auprc,
                    "f1_max": best_f1,
                    "term_fmax": mean_term_fmax,
                    "knn_value_5": knn_value_5,
                    "knn_value_3": knn_value_3,
                    "knn_value_15": knn_value_15,
                    "score_rounded": score_rounded,
                }
            )
            if score_rounded > best_score:
                best_score = score_rounded
                best_loss_epoch = epoch

                # save model
                if directory:
                    import os

                    if os.path.exists(directory):
                        self.save(directory)

            elif epoch - best_loss_epoch > max_epochs_no_improve:
                print(f"Early stopping at epoch {epoch+1}")
                print(
                    f"Last improvement: epoch {best_loss_epoch} with loss {best_score}"
                )
                break
            bar.update(1)

        bar.close()

        if directory:
            generate_pca_gif(directory)

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
            "input_dim": self.input_dim,
            "actual_input_dim": self.actual_input_dim,
            "embedding_dim": self.embedding_dim,
            "metaparams_archive": self.metaparams_archive,
            "vocab_map": self.vocab_map,
        }
        with open(os.path.join(directory, "metadata.pkl"), "wb") as f:
            pickle.dump(meta, f)
        json.dump(
            meta,
            open(os.path.join(directory, "metadata.json"), "w"),
            indent=4,
        )
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
        instance.metaparams_archive = meta["metaparams_archive"]

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
