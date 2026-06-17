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
from sklearn.metrics import (
    silhouette_score,
    davies_bouldin_score,
    calinski_harabasz_score,
)
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

default_intermediary_len = 3600


def plot_history(directory, history_list, output_filename, metric_keys=["fmax"]):
    x_value = "epoch"
    epochs = [hist[x_value] for hist in history_list]
    data_for_plot = {key: [hist[key] for hist in history_list] for key in metric_keys}

    # Create the plot
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot all metrics, which can be many, but tend to be between 0 and 1.
    for key in metric_keys:
        ax1.plot(epochs, data_for_plot[key], "-o", label=key)
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Metric")
    ax1.legend()
    plt.title("Training History")
    output_path = os.path.join(directory, output_filename)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"History plot saved to {output_path}")


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
            if len(overlap) <= 2 and len(overlap) > 0:
                usable_labels = list(overlap)
                usable_labels.sort(key=lambda l: term_counts[l])
                label = usable_labels[0]
                target_indices.append(i)
                cluster_labels.append(label)

    print(f"Selected {len(target_indices)} purely clustered proteins for PCA.")
    different_labels = len(set(cluster_labels))
    print(f"Different labels: {different_labels}")
    return target_indices, cluster_labels, top_terms


"""taxid_pca_labels = {
    "3398": {"name": "Flowering plants", "color": "darkgreen"},
    "3193": {"name": "Land Plants", "color": "forestgreen"},
    "3041": {"name": "Green Algae", "color": "limegreen"},
    "9263": {"name": "Mammals - Marsupials", "color": "red"},
    "9347": {"name": "Mammals - Placentals", "color": "darkred"},
    "8504": {"name": "Lizards and Snakes", "color": "hotpink"},
    "7898": {"name": "Fishes - Actinopterygii", "color": "darkblue"},
    "1783272": {"name": "Bacteria - Bacillati", "color": "yellow"},
    "3379134": {"name": "Bacteria - Pseudomonadati", "color": "gold"},
    "2": {"name": "Bacteria - Others", "color": "khaki"},
    "2732396": {"name": "Viruses - RNA", "color": "darkcyan"},
    "10239": {"name": "Viruses - Others", "color": "black"},
    "other": {"name": "Other Organisms", "color": "magenta"},
}"""
"""taxid_pca_labels = {
    # --- REINO ANIMAL & EUKARYOTA DA ALTA CIÊNCIA ---
    "40674": {
        "name": "Mammals",
        "color": "#B83B5E",
    },  # Tom Terracota/Framboesa profundo
    "7742": {
        "name": "Other Vertebrates",
        "color": "#E23E57",
    },  # Vermelho vivo (Peixes, Aves, Anfíbios, Répteis)
    "6656": {
        "name": "Arthropods & Insects",
        "color": "#FF9A00",
    },  # Laranja/Ouro (Drosophila, etc.)
    "4751": {
        "name": "Fungi",
        "color": "#6A2C70",
    },  # Roxo Escuro/Plum (Leveduras e fungos)
    "33090": {
        "name": "Plants & Algae",
        "color": "#227043",
    },  # Verde Floresta Fechado (Arabidopsis e afins)
    # --- MUNDO PROCARIOTO (BACTÉRIAS) ---
    # Focando nas duas super-potências bacterianas do Swiss-Prot
    "3379134": {
        "name": "Bacteria - Pseudomonadati",
        "color": "#005F73",
    },  # Azul Petróleo (Antiga Proteobacteria / E. coli)
    "1783272": {
        "name": "Bacteria - Bacillati",
        "color": "#0A9396",
    },  # Verde Água Escuro (Antiga Firmicutes / Bacillus)
    "2": {
        "name": "Bacteria - Others",
        "color": "#94D2BD",
    },  # Menta Suave (Demais bactérias)
    # --- VÍRUS ---
    "10239": {
        "name": "Viruses",
        "color": "#2F3E46",
    },  # Cinza Grafite Escuro (Todos os vírus unificados)
    # --- RUÍDO / OUTROS ---
    "other": {
        "name": "Other Organisms",
        "color": "#D3D3D3",
    },  # Cinza Claro Neutro (Evita a poluição visual)
}"""
taxid_pca_labels = {
    # --- EUKARYOTA ---
    "40674": {"name": "Mammals", "color": "darkred"},
    "7742": {"name": "Other Vertebrates", "color": "hotpink"},
    "6656": {"name": "Arthropods & Insects", "color": "#FF9A00"},
    "6231": {"name": "Nematodes", "color": "#E9C46A"},  # C. elegans & worms
    "4751": {"name": "Fungi", "color": "magenta"},
    "33090": {"name": "Plants & Algae", "color": "forestgreen"},
    "5794": {"name": "Apicomplexans", "color": "#A8DADC"},  # Plasmodium / Parasites
    # --- BACTERIA ---
    "3379134": {"name": "Bacteria - Pseudomonadati", "color": "cornflowerblue"},
    "1783272": {"name": "Bacteria - Bacillati", "color": "dodgerblue"},
    "201174": {
        "name": "Bacteria - Actinomycetota",
        "color": "midnightblue",
    },  # Mycobacterium, etc.
    "2": {"name": "Bacteria - Others", "color": "#E0E1DD"},
    # --- VIRUSES ---
    "2732396": {"name": "Viruses - RNA", "color": "#457B9D"},  # Distinct cool blue
    "10239": {
        "name": "Viruses - DNA & Others",
        "color": "#1D3557",
    },  # Deep navy/graphite
    # --- NOISE ---
    "other": {"name": "Other Organisms", "color": "silver"},  # Very light gray
}


def prepare_pca_clusters_taxo(eval_data, top_n=16):
    """
    Finds the top N terms in the evaluation set and assigns proteins
    to strictly mutually exclusive clusters for clean visualization.
    """
    print(f"Preparing PCA clusters for top {top_n} terms...")

    # 1. Find the top N most common terms in the eval set
    all_eval_terms = [term for protein in eval_data for term in protein]

    target_indices = []
    cluster_labels = []

    other_label = "other"

    for i, protein in enumerate(eval_data):
        label = "other"
        for taxid, info in taxid_pca_labels.items():
            if taxid in protein:
                label = taxid
                break
        # pretty_label = taxid_pca_labels[label]["name"]
        cluster_labels.append(label)
        target_indices.append(i)

    undefined_indexes = [i for i in target_indices if cluster_labels[i] == other_label]

    print("Undefined organisms:", len(undefined_indexes) / len(eval_data) * 100, "%")

    print(eval_data[0])
    print(eval_data[-1])
    undefined_eval_data = [
        [x for x in eval_data[i] if x not in taxid_pca_labels]
        for i in undefined_indexes
    ]

    term_list_for_freq = []
    print(undefined_eval_data[0])
    print(undefined_eval_data[-1])
    for lineage in undefined_eval_data:
        for taxid in lineage:
            term_list_for_freq.append(taxid)

    term_counts = Counter(term_list_for_freq)
    top_terms = [term for term, _ in term_counts.most_common(32)]
    print("Top most common terms in undefined organisms:")
    for x in top_terms:
        print(x, term_counts[x])

    print(f"Selected {len(target_indices)} different lineages for PCA.")
    different_labels = len(set(cluster_labels))
    print(f"Different labels: {different_labels}")

    label_counts = Counter(cluster_labels)

    print("Frequencies of used labels:")
    for label, count in label_counts.items():
        print(taxid_pca_labels[label]["name"], count)

    labels_used = list(set(cluster_labels))

    return target_indices, cluster_labels, labels_used


def plot_epoch_pca(
    eval_embeddings, target_indices, cluster_labels, top_terms, epoch, directory
):
    """
    Calculates 2D t-SNE for the selected embeddings and saves a clean, ordered plot.
    """
    import matplotlib.lines as mlines

    X = eval_embeddings[target_indices]

    tsne = TSNE(
        n_components=2, perplexity=30, learning_rate="auto", init="pca", random_state=42
    )
    X_pca = tsne.fit_transform(X)

    fig, ax = plt.subplots(figsize=(4.8, 4.5), dpi=180)  # Slightly wider for the legend

    is_taxid = all([t in taxid_pca_labels for t in cluster_labels])
    if is_taxid:
        colors = {term: taxid_pca_labels[term]["color"] for term in top_terms}
    else:
        cmap = plt.get_cmap("tab20")
        colors = {term: cmap(i) for i, term in enumerate(top_terms)}

    # 1. Scatter plot by cluster
    for term in top_terms:
        term_mask = [label == term for label in cluster_labels]

        if any(term_mask):
            # We skip adding the 'label' argument here to prevent default legend creation
            ax.scatter(
                X_pca[term_mask, 0],
                X_pca[term_mask, 1],
                color=colors[term],
                alpha=0.95,
                s=12,  # Keep the actual plot points small and distinct
                edgecolors="none",
            )

    # 2. Build Custom Legend
    legend_handles = []

    # Iterate through the dictionary to maintain exact taxonomic order
    for taxid, info in taxid_pca_labels.items():
        # Only add to the legend if the organism actually exists in this evaluation batch
        if taxid in top_terms or (taxid == "other" and "other" in cluster_labels):
            handle = mlines.Line2D(
                [],
                [],
                color=info["color"],
                marker="o",
                linestyle="None",
                markersize=9,  # LARGER CIRCLES in the legend
                label=info["name"],
            )
            legend_handles.append(handle)

    # 3. Formatting
    title = (
        f"Compact Representation of Taxonomy - {epoch}"
        if is_taxid
        else f"Compact Representation of Protein Families - {epoch}"
    )
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xticks([])
    ax.set_yticks([])

    # Apply the custom legend
    ax.legend(
        handles=legend_handles,
        fontsize=7,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),  # Push legend entirely outside the plot box
        frameon=False,
        title="Taxonomic Clades",
        title_fontsize=8,
        labelspacing=0.8,  # Add a little breathing room between items
    )

    # 4. Save the plot
    if directory is None:
        directory = "."
    os.makedirs(directory, exist_ok=True)

    filepath = os.path.join(directory, f"pca_epoch_{epoch:03d}.png")
    plt.savefig(filepath, format="png", bbox_inches="tight")
    plt.close(fig)


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


def embedding_neighborhood_score(
    eval_targets, eval_embeddings, k=5, min_overlap_ratio=0.6
):
    """
    Computes the average fraction of k-nearest neighbours that share
    at least `min_overlap_ratio` of terms with the query protein.

    Args:
        eval_targets: torch.Tensor of shape (N, num_terms), binary.
        eval_embeddings: torch.Tensor of shape (N, dim).
        k: number of neighbours to consider.
        min_overlap_ratio: float (0.0 to 1.0). The minimum proportion of the query's
                           terms that the neighbour must also possess to be a "match".

    Returns:
        float: average neighbour label agreement.
    """
    # 1. Normalise embeddings for cosine similarity
    normed_emb = eval_embeddings / (eval_embeddings.norm(dim=1, keepdim=True) + 1e-8)

    # 2. Compute pairwise cosine similarity matrix
    sim = torch.mm(normed_emb, normed_emb.t())  # (N, N)

    # Set diagonal to -inf so the query itself isn't chosen
    sim.fill_diagonal_(-float("inf"))

    # 3. Get top‑k indices for each row
    _, knn_indices = torch.topk(sim, k, dim=1)  # (N, k)

    # 4. Fetch the target vectors for all k neighbours of all N queries
    # Shape: (N, k, num_terms)
    knn_targets = eval_targets[knn_indices]

    # Prepare query targets for broadcasting
    # Shape: (N, 1, num_terms)
    query_targets = eval_targets.unsqueeze(1)

    # 5. Calculate intersection size (number of shared terms)
    # Shape: (N, k)
    intersection_sizes = (query_targets * knn_targets).sum(dim=2)

    # 6. Calculate the number of terms the query actually has
    # Shape: (N, 1)
    query_sizes = query_targets.sum(dim=2)

    # Prevent division by zero if a query has absolutely no labels
    # (Though in taxonomic lineages, this shouldn't happen)
    query_sizes = torch.clamp(query_sizes, min=1.0)

    # 7. Calculate the overlap ratio (Intersection / Query Size)
    # Shape: (N, k)
    overlap_ratios = intersection_sizes / query_sizes

    # 8. Check which neighbours meet the minimum threshold
    # Shape: (N, k) -> boolean mask converted to float (1.0 or 0.0)
    matches = (overlap_ratios >= min_overlap_ratio).float()

    # 9. Calculate the fraction of the k neighbours that are matches for each query
    # Shape: (N,)
    agreement_fracs = matches.mean(dim=1)

    # Return the global mean score
    return agreement_fracs.mean().item()


class InterProAutoencoder(nn.Module):
    def __init__(
        self,
        input_dim,
        embedding_dim,
        dropout_rate=0.1,
        intermediary_len=default_intermediary_len,
    ):
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
    def __init__(
        self,
        input_dim=20000,
        embedding_dim=2000,
        predefined_vocab=None,
        intermediary_len=default_intermediary_len,
        data_family="interpro",
    ):
        self.data_family = data_family
        self.intermediary_len = intermediary_len
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

    def make_sillouette_labels(self, taxa_clusters_for_sillouette, eval_targets):
        """
        Receives taxids of different taxonomic levels (e.g. Kingdom, Phylum, etc.)
        and the one-hot encoded samples. It uses the taxids as clusters and produces
        several labelings efficiently using PyTorch tensor operations.
        Assumes all provided taxids are pre-validated against the vocab_map.
        """
        label_lists = []

        for level_name, cluster_taxids in taxa_clusters_for_sillouette.items():
            # 1. Convert set to list to guarantee stable order and allow indexing
            taxid_list = list(cluster_taxids)

            # 2. Directly map all taxids to their indices using the stable list
            clusters_indexes = [self.vocab_map[taxid] for taxid in taxid_list]

            # 3. Slice the targets: Shape (N, num_clusters_in_level)
            level_targets = eval_targets[:, clusters_indexes]

            # 4. Find the max value (1 if present, 0 if not) and its column index for each row
            max_vals, argmaxes = torch.max(level_targets, dim=1)

            current_labels = [None] * len(eval_targets)

            # 5. Extract rows that actually belong to a cluster at this taxonomic level
            valid_mask = max_vals > 0
            valid_row_indices = torch.where(valid_mask)[0].tolist()
            valid_argmaxes = argmaxes[valid_mask].tolist()

            # 6. Map the valid rows to their actual string taxid labels using the indexable list
            for row_idx, col_idx in zip(valid_row_indices, valid_argmaxes):
                current_labels[row_idx] = taxid_list[col_idx]

            label_lists.append((level_name, current_labels))

        return label_lists

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
        taxa_clusters_for_sillouette: dict = None,
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
            "intermediary_len": self.intermediary_len,
            "dropout_rate": dropout_rate,
        }
        print(f"Training params: {self.metaparams_archive}")

        bar = tqdm(total=epochs)

        print("Building Vocab")
        self._build_vocab(family_lists)
        # bar.update(1)

        print("Initializing Model", file=sys.stderr)
        self.model = InterProAutoencoder(
            self.actual_input_dim, self.embedding_dim, dropout_rate=dropout_rate
        ).to(self.device)
        # bar.update(1)

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

        if self.data_family == "interpro":
            pca_indices, pca_labels, pca_terms = prepare_pca_clusters(
                eval_data, top_n=14
            )
            min_proportion1, min_proportion2 = (0.05, 0.01)
        elif self.data_family == "taxid":
            min_proportion1, min_proportion2 = (0.65, 0.55)
            pca_indices, pca_labels, pca_terms = prepare_pca_clusters_taxo(eval_data)
        else:
            raise Exception("Invalid dta family")

        # pairs_for_sim_calculation, pair_sims = self._calc_ont_sim_matrix(ann_sim_matrix)

        train_data = family_lists

        print("Creating Dataset", file=sys.stderr)
        dataset = InterProDataset(train_data, self.vocab_map, self.actual_input_dim)

        eval_dataset = InterProDataset(eval_data, self.vocab_map, self.actual_input_dim)
        # bar.update(1)

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

        # bar.update(1)

        print("Initializing training", file=sys.stderr)
        if pos_weight_value is not None:
            pos_weight = torch.tensor([pos_weight_value]).to(self.device)
            criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        else:
            criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(self.model.parameters(), lr=lr)

        self.model.train()
        # bar.update(1)

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

                # Make sillouette labels
                eval_labels_list = self.make_sillouette_labels(
                    taxa_clusters_for_sillouette, eval_targets
                )

                silhouette_scores = {}
                calinski_scores = {}
                davies_scores = {}

                eval_embeddings_np = eval_embeddings.numpy()

                for tax_level, labeling in eval_labels_list:
                    # print(f"Calculating silhouette score at level {tax_level}...")

                    # Filter out samples that lack a label at this taxonomic level
                    non_null_indexes = [
                        i for i, l in enumerate(labeling) if l is not None
                    ]

                    if len(non_null_indexes) < 2:
                        silhouette_scores[tax_level] = 0.0
                        continue

                    eval_pos = eval_embeddings_np[non_null_indexes]
                    eval_labels = [labeling[i] for i in non_null_indexes]

                    # Silhouette score strictly requires at least 2 distinct clusters
                    if len(set(eval_labels)) > 1:
                        sil_score = (
                            1.0
                            + silhouette_score(
                                eval_pos, eval_labels, metric="euclidean"
                            )
                        ) / 2.0
                        cal_score = (
                            calinski_harabasz_score(eval_pos, eval_labels) / 700.0
                        )
                        dav_score = 1.0 - np.log10(
                            davies_bouldin_score(eval_pos, eval_labels) + 1
                        )
                    else:
                        sil_score = 0.0
                        cal_score = 0.0
                        dav_score = 0.0

                    silhouette_scores[tax_level] = float(sil_score)
                    calinski_scores[tax_level] = float(cal_score)
                    davies_scores[tax_level] = float(dav_score)

                # Calculate the mean silhouette score across all valid taxonomic levels
                if len(silhouette_scores) > 0:
                    mean_silhouette = sum(silhouette_scores.values()) / len(
                        silhouette_scores
                    )
                    mean_calinski = sum(calinski_scores.values()) / len(calinski_scores)
                    mean_davies = sum(davies_scores.values()) / len(davies_scores)
                else:
                    mean_silhouette = 0.0
                    mean_calinski = 0.0
                    mean_davies = 0.0

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
                    eval_targets,
                    eval_embeddings,
                    k=5,
                    min_overlap_ratio=min_proportion2,
                )

                knn_value_3 = embedding_neighborhood_score(
                    eval_targets,
                    eval_embeddings,
                    k=3,
                    min_overlap_ratio=min_proportion1,
                )

                loss_rounded = round(total_loss, 6)
                loss_norm = 1.0 - loss_rounded
                if loss_norm < 0:
                    loss_norm = 0.0

                # Update Epoch score to use your new Term-centric metric instead of Instance F1
                metric_ws = [
                    1,
                    2,
                    1,
                    1,
                ]
                metric_vals = [
                    loss_norm,
                    # mean_term_fmax,
                    best_f1,
                    mean_silhouette,
                    mean_davies,
                ]
                epoch_score = float(np.average(metric_vals, weights=metric_ws))
                score_rounded = round(epoch_score, 4)
                s_scores = "; ".join(
                    f"Level {l}: {round(s, 3)}" for l, s in silhouette_scores.items()
                )
                c_scores = "; ".join(
                    f"Level {l}: {round(s, 3)}" for l, s in calinski_scores.items()
                )
                d_scores = "; ".join(
                    f"Level {l}: {round(s, 3)}" for l, s in davies_scores.items()
                )

                print(
                    f"kNN-3: {knn_value_3:.4f} | "
                    f"kNN-5: {knn_value_5:.4f} | "
                    f"\nSilhouette Scores: {s_scores} | "
                    f"\nCalinski Scores: {c_scores} | "
                    f"\nMean Silh: {mean_silhouette:.4f} | "
                    f"Mean Calinski: {mean_calinski:.4f} | "
                    f"Mean Davies: {mean_davies:.4f} | "
                    f"Term Fmax: {best_f1:.4f} | "
                    f"Loss: {total_loss/len(dataloader):.6f} | "
                    f"\nEpoch Score: {score_rounded:.4f}"
                )

            # Update the history tracker
            self.history.append(
                {
                    "epoch": epoch,
                    "loss": loss_rounded,
                    "loss_norm": loss_norm,
                    "best_score": best_score,
                    "auprc": auprc,
                    "f1_max": best_f1,
                    "term_fmax": best_f1,
                    "knn_value_5": knn_value_5,
                    "knn_value_3": knn_value_3,
                    "mean_silhouette": mean_silhouette,
                    "mean_davies": mean_davies,
                    "mean_calinski": mean_calinski,
                    "score_rounded": score_rounded,
                }
            )

            if directory:
                plot_history(
                    directory,
                    self.history,
                    "history.png",
                    metric_keys=[
                        "loss_norm",
                        "term_fmax",
                        "score_rounded",
                        "mean_silhouette",
                        "mean_davies",
                        "mean_calinski",
                        "auprc",
                    ],
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
