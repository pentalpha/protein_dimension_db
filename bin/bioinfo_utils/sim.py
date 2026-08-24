from collections import defaultdict
import os

from tqdm import tqdm
import numpy as np


class DictBasedTree:
    def __init__(self, obo_path):
        self.term_ancestors = defaultdict(set)
        self.term_descendants = defaultdict(set)
        self.parent_to_child = defaultdict(set)
        self.child_to_parents = defaultdict(set)
        self.id_to_names = defaultdict(list)
        self.all_ids = set()
        self._load_obo(obo_path)

    def _load_obo(self, path):
        """Builds a map of all ancestors for every term."""
        child_to_parents = defaultdict(set)
        parent_to_child = defaultdict(set)
        all_ids = set()
        with open(path, "r") as f:
            curr_id = None
            for line in f:
                if line.startswith("id: "):
                    curr_id = line.strip().split("id: ")[1]
                    all_ids.add(curr_id)
                elif line.startswith("name: "):
                    names = line.strip().split("name: ")[1].split(",")
                    names = [name.strip() for name in names]
                    self.id_to_names[curr_id] = names
                elif line.startswith("is_a: ") or line.startswith("parent: "):
                    parent = line.strip().split(": ")[1].split(" ! ")[0]
                    child_to_parents[curr_id].add(parent)
                    parent_to_child[parent].add(curr_id)

        def get_ancestors(tid, visited):
            if tid in self.term_ancestors:
                return self.term_ancestors[tid]
            anc = {tid}
            for p in child_to_parents.get(tid, []):
                if p not in visited:
                    anc.update(get_ancestors(p, visited | {p}))
            return anc

        def get_descendants(tid, visited):
            if tid in self.term_descendants:
                return self.term_descendants[tid]
            desc = {tid}
            for c in parent_to_child.get(tid, []):
                if c not in visited:
                    desc.update(get_descendants(c, visited | {c}))
            return desc

        print("Indexing ontology ancestors...")
        for tid in all_ids:
            self.term_ancestors[tid] = get_ancestors(tid, set())
            self.term_descendants[tid] = get_descendants(tid, set())
        self.child_to_parents = child_to_parents
        self.parent_to_child = parent_to_child
        self.all_ids = all_ids

    def get_id_by_name(self, name):
        for tid, names in self.id_to_names.items():
            if name in names:
                return tid
        return None

    def get_clusters_abrangentes(self, freq_dict: dict, vocab_list: list, n: int):
        vocab = sorted(vocab_list, key=lambda x: freq_dict.get(x, 0), reverse=True)
        vocab_set = set(vocab)
        to_ignore = set()
        clusters = []

        for label in vocab:
            if label in to_ignore:
                continue
            clusters.append(label)
            descendants = self.term_descendants.get(label, set())
            for d in descendants:
                if d in vocab_set and d not in to_ignore:
                    to_ignore.add(d)

            if len(clusters) == n:
                break

        names = [self.id_to_names.get(t, [t])[0] for t in clusters]

        return clusters, names


def find_intersecting_terms(eval_data, vocab):
    print("Finding intersections between terms to draw networks")
    intersections = {}
    for terms in tqdm(eval_data):
        in_vocab = [t for t in terms if t in vocab]
        if len(in_vocab) > 1:
            for t in in_vocab:
                if t not in intersections:
                    intersections[t] = set()
                others = [x for x in in_vocab if x != t]
                intersections[t].update(others)
    return intersections


class OntologyTermSimilarity:
    def __init__(self, vocab_list, ia_data, obo_path, save=True):
        self.ia_data = ia_data  # Dictionary {term: IC}
        self.term_ancestors = defaultdict(set)
        self._load_obo(obo_path)
        self.vocab_list = vocab_list
        self.vocab_map = {v: i for i, v in enumerate(self.vocab_list)}
        self.term_sim_matrix_path = obo_path.replace(".obo", ".npy")
        if os.path.exists(self.term_sim_matrix_path):
            self.term_sim_matrix = np.load(self.term_sim_matrix_path)
        else:
            self.term_sim_matrix = self.get_term_term_matrix(self.vocab_list)
            np.save(self.term_sim_matrix_path, self.term_sim_matrix)

    def _load_obo(self, path):
        """Builds a map of all ancestors for every term."""
        child_to_parents = defaultdict(list)
        all_ids = set()
        with open(path, "r") as f:
            curr_id = None
            for line in f:
                if line.startswith("id: "):
                    curr_id = line.strip().split("id: ")[1]
                    all_ids.add(curr_id)
                elif line.startswith("is_a: ") or line.startswith("parent: "):
                    parent = line.strip().split(": ")[1].split(" ! ")[0]
                    child_to_parents[curr_id].append(parent)

        def get_ancestors(tid, visited):
            if tid in self.term_ancestors:
                return self.term_ancestors[tid]
            anc = {tid}
            for p in child_to_parents.get(tid, []):
                if p not in visited:
                    anc.update(get_ancestors(p, visited | {p}))
            return anc

        print("Indexing ontology ancestors...")
        for tid in all_ids:
            self.term_ancestors[tid] = get_ancestors(tid, set())

    def get_term_term_matrix(self, vocab_list):
        """Calculates Lin Similarity between all terms in the vocabulary."""
        T = len(vocab_list)
        matrix = np.zeros((T, T), dtype=np.float32)

        # Pre-cache IC values and ancestor sets for the specific vocab
        ics = [self.ia_data.get(t, 0.0) for t in vocab_list]
        ancestors = [self.term_ancestors.get(t, {t}) for t in vocab_list]

        print(f"Calculating {T}x{T} Ontology Similarity Matrix...")
        bar = tqdm(total=T)
        for i in range(T):
            for j in range(i, T):
                # Find MICA (Most Informative Common Ancestor)
                common = ancestors[i].intersection(ancestors[j])
                if not common:
                    continue

                mica_ic = max([self.ia_data.get(anc, 0.0) for anc in common])

                # Lin's Formula
                denom = ics[i] + ics[j]
                sim = (2.0 * mica_ic) / denom if denom > 0 else 0.0

                matrix[i, j] = matrix[j, i] = sim
            bar.update(1)
        bar.close()

        return matrix

    def calculate_bma_matrix(self, protein_term_lists):
        """
        Computes Protein-Protein similarity using Best-Match Average.
        protein_term_lists: List of lists of InterPro terms
        """
        N = len(protein_term_lists)
        sim_matrix = np.zeros((N, N), dtype=np.float32)

        # Convert terms to indices for fast lookup
        protein_indices = [
            [self.vocab_map[t] for t in p if t in self.vocab_map]
            for p in protein_term_lists
        ]

        print("Calculating Protein-Protein BMA matrix...")
        bar = tqdm(total=N)
        for i in range(N):
            idxs_a = protein_indices[i]
            if not idxs_a:
                continue
            for j in range(i, N):
                idxs_b = protein_indices[j]
                if not idxs_b:
                    continue

                # Sub-matrix of similarities between terms of Prot A and Prot B
                sub = self.term_sim_matrix[np.ix_(idxs_a, idxs_b)]

                # BMA Formula
                score_ab = np.mean(np.max(sub, axis=1))
                score_ba = np.mean(np.max(sub, axis=0))

                sim_matrix[i, j] = sim_matrix[j, i] = (score_ab + score_ba) / 2.0
            bar.update(1)
        bar.close()

        return sim_matrix


class SemanticSimilarity:
    def __init__(self, ia_data):
        """
        vocab_map: your existing {term: index} dictionary from the Autoencoder.
        """

        # Create a weight vector aligned with your model's input indices
        self.ia_data = ia_data

        print(f"Loaded semantic weights for {len(ia_data)} terms.")

    def calculate_simgic_matrix(self, binary_vectors, weights):
        """
        Calculates the SimGIC similarity matrix for a batch of proteins.
        binary_vectors: Tensor of shape (N_proteins, N_terms)
        """
        # Ensure weights are on the same device as data
        w = weights

        # 1. Weighted Intersection: sum(ia * (A & B))
        # We can do this via matrix multiplication: (N, T) @ (T, N)
        # But we must multiply the rows by weights first.
        weighted_vectors = binary_vectors * w
        intersection = np.dot(weighted_vectors, binary_vectors.T)

        # 2. Weighted Individual Sums: sum(ia * A)
        individual_sums = (binary_vectors * w).sum(axis=1)

        # 3. Weighted Union: sum(ia * A) + sum(ia * B) - intersection
        # Use broadcasting to get all pairs: (N, 1) + (1, N)
        union = (
            individual_sums[:, np.newaxis]
            + individual_sums[np.newaxis, :]
            - intersection
        )

        # 4. SimGIC
        similarity_matrix = intersection / (union + 1e-8)

        return similarity_matrix

    def calculate_simgic_matrix_with_onehot(self, term_lists):

        all_terms = set()
        for ann in term_lists:
            all_terms.update(ann)

        terms_sequence = list(all_terms)
        terms_indexes = {term: i for i, term in enumerate(terms_sequence)}
        weights = np.array([self.ia_data[term] for term in terms_sequence])

        binary_vectors = []
        for ann in term_lists:
            vector = np.zeros(len(terms_sequence))
            for term in ann:
                vector[terms_indexes[term]] = 1
            binary_vectors.append(vector)
        binary_vectors = np.asarray(binary_vectors)

        similarity_matrix = self.calculate_simgic_matrix(binary_vectors, weights)

        return similarity_matrix


def calc_ont_sim_matrix(ann_sim_matrix):

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
