import os
import pickle
import numpy as np
from collections import Counter


class OneHotEncoder:
    # Similar structure to the autoencoder, but simply builds a vocab and returns one-hot vectors of a certain max-length
    # Does not use torch, any machine learning or other heavy dependencies, just simple data structures
    def __init__(self, max_families=10000):
        self.max_families = max_families
        self.vocab_map = None

    def _build_vocab(self, family_lists):
        """Creates a mapping for the top N most frequent InterPro families."""
        print(f"Building vocabulary for top {self.max_families} families...")
        all_fams = [f for sublist in family_lists for f in sublist]
        counter = Counter(all_fams)
        print(f"Found {len(counter)} unique families.")

        most_common = counter.most_common(self.max_families)
        print(f"Using top {len(most_common)} families.")

        # Sort families alphabetically for deterministic saving/loading
        top_families = sorted([token for token, count in most_common])

        if len(top_families) > 0:
            print(
                f"First 5 families (alphabetical): {top_families[:min(5, len(top_families))]}"
            )
            print(
                f"Last 5 families (alphabetical): {top_families[-min(5, len(top_families)):]}"
            )

        self.vocab_map = {token: i for i, token in enumerate(top_families)}

        fams_set = set(self.vocab_map.keys())
        fams_count = sum(
            1
            for fam_list in family_lists
            if len(set(fam_list).intersection(fams_set)) > 0
        )
        perc = fams_count / len(family_lists) * 100 if len(family_lists) > 0 else 0
        print(
            f"Percentage of proteins with at least one family in the vocabulary: {perc:.2f}%"
        )

    def fit(self, family_lists):
        """Builds vocabulary based on family lists, mirroring the fit interface."""
        if self.vocab_map is None:
            self._build_vocab(family_lists)
        return self

    def predict(self, family_lists):
        """Encodes a list of lists into one-hot (multi-hot) vectors."""
        if self.vocab_map is None:
            raise ValueError("Model not trained or loaded. Call fit() first.")

        vocab_size = len(self.vocab_map)
        result = np.zeros((len(family_lists), vocab_size), dtype=np.float32)
        for i, families in enumerate(family_lists):
            for f in families:
                if f in self.vocab_map:
                    result[i, self.vocab_map[f]] = 1.0
        return result

    def save(self, directory):
        """Saves the vocabulary as a strictly sorted text list."""
        if self.vocab_map is None:
            raise ValueError("Model not trained or loaded. Call fit() first.")
        os.makedirs(directory, exist_ok=True)

        # Ordered appropriately based on mapping index
        families = [
            token for token, idx in sorted(self.vocab_map.items(), key=lambda x: x[1])
        ]

        with open(os.path.join(directory, "vocab.txt"), "w") as f:
            for fam in families:
                f.write(f"{fam}\n")
        print(f"Vocabulary saved to {os.path.join(directory, 'vocab.txt')}")

    @classmethod
    def load(cls, directory):
        """Loads a saved instance from a text list."""
        with open(os.path.join(directory, "vocab.txt"), "r") as f:
            families = [line.strip() for line in f if line.strip()]

        instance = cls(max_families=len(families))
        instance.vocab_map = {token: i for i, token in enumerate(families)}
        return instance
