import os
import uuid
from glob import glob

import polars as pl
from torch.utils.data import Dataset

from bioinfo_utils.fasta import read_fasta

def batch_by_tokens(seqs, max_tokens_per_batch=1800, use_padding = True):
    batchs = []
    new_batch = []
    for seq in seqs:
        if len(new_batch) > 0:
            if use_padding:
                current_batch_len = max([len(s) for s in new_batch]) * len(new_batch)
            else:
                current_batch_len = sum([len(s) for s in new_batch])
        else:
            current_batch_len = 0
        if current_batch_len + len(seq) > max_tokens_per_batch:
            if current_batch_len == 0:
                new_batch = [seq]
            else:
                batchs.append(new_batch)
                new_batch = [seq]
        else:
            new_batch.append(seq)
    if len(new_batch) > 0:
        batchs.append(new_batch)
    return batchs

class FastaDataset(Dataset):
    def __init__(self, fasta_path: str):
        """
        A class for loading a FASTA file.

        Args:
            fasta_path: String specifying the path of the FASTA file
        """
        non_embedded_seqs = read_fasta(fasta_path)
        non_embedded_seqs = [seq for _, seq in non_embedded_seqs]

        self.sequences = [
            list(line) for line in non_embedded_seqs
        ]

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, index):
        return self.sequences[index]

class EmbCache:
    def __init__(self, cache_path, model_name, full=False):
        self.cache_path = cache_path
        self.model_name = model_name
        self.model_safe_name = model_name.replace("/", "__").lower()
        if full:
            self.model_safe_name += ".full"
        self.save_paths_expr_pq = cache_path + "/" + self.model_safe_name + ".*.parquet"

    def list_caches(self):
        pqts = glob(self.save_paths_expr_pq)
        valid_caches = []
        for p in pqts:
            txt_path = p.replace(".parquet", ".txt")
            if os.path.exists(txt_path):
                try:
                    # Extremely fast integrity check (reads footer only)
                    pl.scan_parquet(p)
                    valid_caches.append((p, txt_path))
                except Exception as e:
                    print(f"⚠️ Corrupted cache detected: {p} - Error: {e}")
                    print("Deleting broken files to trigger re-computation...")
                    os.remove(p)
                    if os.path.exists(txt_path):
                        os.remove(txt_path)
        return valid_caches

    def next_cache_name(self):
        """index_name = 0
        path_full = (
            self.cache_path + "/" + self.model_safe_name + "." + str(index_name) + "."
        )
        while os.path.exists(path_full + ".txt") and os.path.exists(
            path_full + ".parquet"
        ):
            index_name += 1
            path_full = (
                self.cache_path
                + "/"
                + self.model_safe_name
                + "."
                + str(index_name)
                + "."
            )"""
        uid = uuid.uuid4().hex[:8]
        path_base = f"{self.cache_path}/{self.model_safe_name}.{uid}"
        return path_base + ".parquet", path_base + ".txt"

    def list_embedded(self):
        caches = self.list_caches()
        embedded_seqs = set()
        for pq_path, txt_path in caches:
            with open(txt_path, "r") as f:
                for line in f:
                    protein_seq = line.strip()
                    embedded_seqs.add(protein_seq)
        return embedded_seqs

    def list_non_embedded(self, fasta_path):
        fasta_content = read_fasta(fasta_path)
        embedded_seqs = self.list_embedded()
        non_embedded_seqs = set()
        embedded_count = 0
        non_embedded_count = 0
        for seq_id, seq in fasta_content:
            if seq in embedded_seqs:
                embedded_count += 1
            else:
                non_embedded_seqs.add(seq)
                non_embedded_count += 1
        non_embedded_seqs = list(non_embedded_seqs)
        non_embedded_seqs.sort(key=len, reverse=False)
        print(f"Embedded: {embedded_count}, Non-embedded: {non_embedded_count}")
        perc_done = embedded_count / (embedded_count + non_embedded_count) * 100
        print(f"Percentage done: {perc_done:.2f}%")
        return non_embedded_seqs