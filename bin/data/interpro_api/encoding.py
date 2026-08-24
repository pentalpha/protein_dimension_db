import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import os

from tqdm import tqdm
import numpy as np
import sys
import os
import pickle
import json
from collections import Counter
import multiprocessing

import torch.multiprocessing as mp

import matplotlib.pyplot as plt
from collections import Counter
from sklearn.metrics import (
    silhouette_score,
    davies_bouldin_score,
    calinski_harabasz_score,
)

from bioinfo_utils.using_torch import embedding_neighborhood_score
from bioinfo_utils.clustering import (
    prepare_pca_clusters_predef,
)
from bioinfo_utils.sim import (
    DictBasedTree,
    find_intersecting_terms,
)
from plotting_lib.scatter import (
    prepare_pca_clusters_taxo,
    plot_epoch_pca,
    generate_pca_gif,
    plot_tsne_scatter_multilabel_pies,
)

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


def validate_clusters_in_vocab(clusterings, vocab_map):
    """
    Clusterings are given as a dictionary of lists of terms.
    Check if the terms exist in the vocab_map.
    """
    filtered_clusterings = {}
    for name, terms in clusterings.items():
        filtered_terms = [term for term in terms if term in vocab_map]
        print(f"Proposed list of labels had {len(terms)} terms")
        print(f"Filtered list of labels has {len(filtered_terms)} terms")

        if len(filtered_terms) > 0:
            filtered_clusterings[name] = filtered_terms
        else:
            print(f"No terms in {name} found in vocab_map")
    return filtered_clusterings


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
        max_epochs_no_improve=4,
        directory=None,
        optimize_cpu=False,
        eval_perc=0.3333,
        taxa_clusters_for_sillouette=None,
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

        obo_tree = DictBasedTree(obo_path)

        import os

        if taxa_clusters_for_sillouette is not None:
            taxa_clusters_for_sillouette = validate_clusters_in_vocab(
                taxa_clusters_for_sillouette, self.vocab_map
            )

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

        vocab_list = sorted(
            self.vocab_map.keys(), key=lambda x: ia_weights[x], reverse=True
        )
        list_for_counting = []
        for l in eval_data:
            list_for_counting += l

        vocab_freqs = Counter(list_for_counting)

        if self.data_family == "interpro":
            # intersections_ref = find_intersecting_terms(eval_data, vocab_list)
            cluster_ids, cluster_names = obo_tree.get_clusters_abrangentes(
                vocab_freqs, vocab_list, n=40
            )
            # Preference to more specific terms
            cluster_ids = sorted(
                cluster_ids, key=lambda x: vocab_freqs.get(x, 0), reverse=False
            )
            pca_indices, pca_labels, pca_mlabels, pca_terms = (
                prepare_pca_clusters_predef(
                    eval_data, cluster_ids, include_unclustered=True
                )
            )
            cluster_names = {x: y for x, y in zip(cluster_ids, cluster_names)}
            cluster_names["Other"] = "Other"
            pca_terms.append("Other")

            label_counts = Counter(pca_labels)

            print("Clusters abrangentes de Interpro:")
            for k, v in cluster_names.items():
                print(f"{k}: {v} ({label_counts[k]} / {vocab_freqs.get(k, 0)} items)")

            print(sorted(pca_terms))
            print(sorted(list(label_counts.keys())))

            pca_terms = sorted(pca_terms, key=lambda x: label_counts[x], reverse=True)

            """if taxa_clusters_for_sillouette is not None:
                pca_indices, pca_labels, pca_terms = prepare_pca_clusters_predef(
                    eval_data,
                    taxa_clusters_for_sillouette["freq_91_to_400"],
                )
            else:
                pca_indices, pca_labels, pca_terms = prepare_pca_clusters(
                    eval_data, top_n=14
                )"""
            min_proportion1, min_proportion2 = (0.05, 0.01)
        elif self.data_family == "taxid":
            intersections_ref = {}
            min_proportion1, min_proportion2 = (0.65, 0.55)
            pca_indices, pca_labels, pca_terms = prepare_pca_clusters_taxo(eval_data)
            pca_mlabels = pca_labels
            cluster_names = {str(x): x for x in pca_terms}
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

        plotting_metainfo = {
            "pca_indices": pca_indices,
            "pca_labels": pca_labels,
            "pca_mlabels": pca_mlabels,
            "pca_terms": pca_terms,
            "pretty_names": cluster_names,
            "eval_data": [";".join(x) for x in eval_data],
        }

        if directory is not None:
            os.makedirs(directory, exist_ok=True)
            json.dump(
                plotting_metainfo,
                open(os.path.join(directory, "plotting_metainfo.json"), "w"),
                indent=2,
                ensure_ascii=False,
            )

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
                    eval_embs_data = torch.sigmoid(eval_embeddings).numpy()
                    epoch_eval_state_path = f"{directory}/eval_embs_epoch_{epoch+1}.npy"
                    np.save(epoch_eval_state_path, eval_embs_data)
                    if self.data_family == "interpro":
                        if (epoch + 1) % 2 == 0 or epoch == (epochs - 1):
                            plot_tsne_scatter_multilabel_pies(
                                eval_embs_data,
                                pca_indices,
                                pca_labels,
                                pca_terms,
                                epoch,
                                directory,
                                pretty_names=cluster_names,
                            )
                    else:
                        plot_epoch_pca(
                            eval_embs_data,
                            pca_indices,
                            pca_labels,
                            pca_terms,
                            epoch,
                            directory,
                        )

                    if (epoch + 1) % 6 == 0 or epoch == (epochs - 1):
                        generate_pca_gif(directory)

                """knn_value_5 = embedding_neighborhood_score(
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
                )"""

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
                    # f"kNN-3: {knn_value_3:.4f} | "
                    # f"kNN-5: {knn_value_5:.4f} | "
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
                    "global_fmax": best_f1,
                    "term_fmax_mean": mean_term_fmax,
                    "mean_silhouette": mean_silhouette,
                    "mean_davies": mean_davies,
                    "mean_calinski": mean_calinski,
                    "score_rounded": score_rounded,
                    # "auprc": auprc,
                    # "knn_value_5": knn_value_5,
                    # "knn_value_3": knn_value_3,
                }
            )

            if directory:
                plot_history(
                    directory,
                    self.history,
                    "history.png",
                    metric_keys=[
                        "loss_norm",
                        "global_fmax",
                        "score_rounded",
                        "mean_silhouette",
                        "mean_davies",
                        "mean_calinski",
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
