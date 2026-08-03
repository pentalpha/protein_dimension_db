#!/usr/bin/env python3
import sys
import os
import time
from glob import glob
import uuid

from transformers import AutoModel
import polars as pl
import torch
from tqdm import tqdm

original_setattr = type(torch._dynamo.config).__setattr__


def custom_setattr(self, name, value):
    if name == "recompile_limit":
        name = "cache_size_limit"
    return original_setattr(self, name, value)


type(torch._dynamo.config).__setattr__ = custom_setattr
# --- END PATCH ---

# List cuda available devices
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        print(f"CUDA Device {i}: {torch.cuda.get_device_name(i)}")
    device = "cuda"
else:
    print("CUDA is not available")
    device = "cpu"
    quit(1)

from bioinfo_utils.fasta import fasta_equal_split_by_len, read_fasta

"""
#TODO
Split fasta into N chunks and run embed_dataset on each chunk separately.
Create final merge with same order of sequences as fasta file.
Each chunk should be run with it's own max_len (actual max len in that chunk).
"""

nvidia3050_max_tokens = {
    "Synthyra/ANKH_base": 9000,
    "Synthyra/ANKH_large": 680,
    "Synthyra/ANKH2_large": 680,
    "Synthyra/ANKH3_large": 680,
    "Synthyra/ANKH3_xl": 320,
    "default": 12000,
}

nvidiav100_max_tokens = {
    "Synthyra/ANKH_base": 18000,
    "Synthyra/ANKH_large": 14000,
    "Synthyra/ANKH2_large": 14000,
    "Synthyra/ANKH3_large": 14000,
    "Synthyra/ANKH3_xl": 6000,
    "default": 12000,
}

poolings = ["mean", "std", "max", "parti", "cls", "norm"]

default_max_tokens = nvidiav100_max_tokens
N_PARTS = 1337
CACHE_DESIRED_N_TOKENS = 1800 * 300
max_processing_time_secs = 60 * 25
min_processing_time_secs = 60 * 10


class EmbCache:
    def __init__(self, cache_path, model_name):
        self.cache_path = cache_path
        self.model_name = model_name
        self.model_safe_name = model_name.replace("/", "__").lower()
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
        non_embedded_seqs.sort(key=len, reverse=True)
        print(f"Embedded: {embedded_count}, Non-embedded: {non_embedded_count}")
        perc_done = embedded_count / (embedded_count + non_embedded_count) * 100
        print(f"Percentage done: {perc_done:.2f}%")
        return non_embedded_seqs


def embed_to_cache(
    seqs: list,
    model,
    max_tokens_per_batch: int,
    cache_pqt: str,
    cache_txt: str,
):
    n_poolings = len(poolings)
    time_start = time.time()

    max_len = max(len(s) for s in seqs)

    optimal_batch_len_float = max_tokens_per_batch / max_len
    processing_batch_size = max(int(optimal_batch_len_float), 1)
    print(f"Setting max len to actual max len in the dataset: {max_len}")
    print(f"Processing batch size: {processing_batch_size} ({optimal_batch_len_float})")

    embeddings = model.embed_dataset(
        sequences=seqs,
        batch_size=processing_batch_size,
        pooling_types=poolings,
        max_len=max_len,
        save=False,
        # save_path=save_path,
    )

    emb_example = next(iter(embeddings.values()))
    fully_pooled_len = emb_example.shape[0]
    assert fully_pooled_len % n_poolings == 0
    actual_emb_width = fully_pooled_len // n_poolings
    pooling_starts = {}
    print(f"Model output shape: {emb_example.shape}")
    for pooling_i in range(n_poolings):
        pooling_name = poolings[pooling_i]
        start = pooling_i * actual_emb_width
        end = pooling_i * actual_emb_width + actual_emb_width - 1
        pooling_starts[pooling_name] = (start, end)
        print(f"Pooling {poolings[pooling_i]} starts at {start} and ends at {end}")

    seqs_embedded = []
    embs = {p: [] for p in poolings}
    for seq, emb in embeddings.items():
        # print(seq)
        # print(emb.shape)
        # print(emb[0].shape)
        # print(emb[-1].shape)
        # torch.Tensor
        # print(type(emb))
        concat_emb_np = emb.cpu().numpy()
        for pooling_name in poolings:
            start, end = pooling_starts[pooling_name]
            emb_np = concat_emb_np[start : end + 1]
            # print(pooling_name, emb_np.shape)
            embs[pooling_name].append(emb_np)
        seqs_embedded.append(seq)

    df = pl.DataFrame({"seq": seqs_embedded, **embs})
    df.write_parquet(cache_pqt)
    with open(cache_txt, "w") as f:
        for seq in seqs_embedded:
            f.write(f"{seq}\n")
    print(f"Saved cache to {cache_pqt} and {cache_txt}")

    duration = time.time() - time_start
    print(f"Duration: {duration}")

    return duration


if __name__ == "__main__":
    fasta_path = sys.argv[1]
    cache_path = sys.argv[2]
    model_name = sys.argv[3]  # Synthyra/ANKH_base
    parquet_name = sys.argv[4]
    MAX_TOKENS_PER_BATCH = default_max_tokens.get(
        model_name, default_max_tokens["default"]
    )
    print(f"MAX_TOKENS_PER_BATCH: {MAX_TOKENS_PER_BATCH}, N_PARTS: {N_PARTS}")

    print(f"Using model: {model_name}")

    model = AutoModel.from_pretrained(model_name, trust_remote_code=True).to("cuda")
    cache = EmbCache(cache_path, model_name)
    durations = []
    non_embedded_seqs = cache.list_non_embedded(fasta_path)
    while len(non_embedded_seqs) > 0:
        next_to_embed = []
        while (
            sum([len(s) for s in next_to_embed]) < CACHE_DESIRED_N_TOKENS
            and len(non_embedded_seqs) > 0
        ):
            next_to_embed.append(non_embedded_seqs.pop(0))
        next_parquet, next_txt = cache.next_cache_name()
        try:
            duration = embed_to_cache(
                next_to_embed, model, MAX_TOKENS_PER_BATCH, next_parquet, next_txt
            )
        except Exception as e:
            print(f"Error embedding sequences: {e}")
            if os.path.exists(next_parquet):
                os.remove(next_parquet)
            if os.path.exists(next_txt):
                os.remove(next_txt)
            raise (e)

        success = os.path.exists(next_parquet) and os.path.exists(next_txt)
        if success:
            if duration > max_processing_time_secs:
                print(f"Slow batch! {duration}s > {max_processing_time_secs}s")
                CACHE_DESIRED_N_TOKENS = int(CACHE_DESIRED_N_TOKENS * 0.6)
                print(f"New CACHE_DESIRED_N_TOKENS: {CACHE_DESIRED_N_TOKENS}")
            elif duration < min_processing_time_secs:
                print(f"Fast batch! {duration}s < {min_processing_time_secs}s")
                CACHE_DESIRED_N_TOKENS = int(CACHE_DESIRED_N_TOKENS * 1.5)
                print(f"New CACHE_DESIRED_N_TOKENS: {CACHE_DESIRED_N_TOKENS}")
            n_tokens_processed = sum([len(s) for s in next_to_embed])
            print(f"Processed {n_tokens_processed} tokens in {duration} seconds")
            tokens_per_sec = n_tokens_processed / duration
            durations.append(tokens_per_sec)
            last_durations = durations[-10:]
            mean_tokens_per_sec = sum(last_durations) / len(last_durations)
            non_embedded_seqs = cache.list_non_embedded(fasta_path)

            remaining_tokens = sum([len(s) for s in non_embedded_seqs])
            pred_total_time = remaining_tokens / mean_tokens_per_sec
            pred_total_time_min = pred_total_time / 60

            print(
                f"Predicted remaining time: {pred_total_time_min:.2f}min ({remaining_tokens} tokens at {mean_tokens_per_sec:.2f} tokens/sec)"
            )
        else:
            print("Error embedding sequences: cache files not found")

    print("Reading original fasta file")
    fasta_content = read_fasta(fasta_path)
    ids = [seq_id for seq_id, _ in fasta_content]
    seqs = [seq for _, seq in fasta_content]

    # 1. Create a Polars DataFrame of the exact original FASTA order
    # .with_row_index gives us an integer column we can use to restore the order later
    df_fasta = pl.DataFrame({"id": ids, "seq": seqs}).with_row_index("original_order")

    all_cache_files = cache.list_caches()
    valid_pqs = [pq_path for pq_path, txt_path in all_cache_files]

    print("Lazy loading and merging embeddings via Polars...")

    # 2. Lazily scan all parquet chunks simultaneously (virtually zero memory)
    df_embs = pl.scan_parquet(valid_pqs)

    # Deduplicate embeddings just in case a failed job left duplicate sequences
    df_embs = df_embs.unique(subset=["seq"], keep="first")

    # 3. Build the computation graph: Join the FASTA frame with the embeddings, then sort
    df_final_lazy = (
        df_fasta.lazy()
        .join(df_embs, on="seq", how="left")
        .sort("original_order")
        .drop("original_order")
    )

    print("Executing join, sort, and collecting into memory...")
    # 4. .collect() executes the highly optimized Rust code.
    # This will use a fraction of the memory compared to native Python lists.
    df_final = df_final_lazy.collect()

    print("Checking for missing embeddings...")
    # If a sequence from the FASTA wasn't found in the caches, the left-join will leave nulls
    n_missing = df_final.filter(pl.col(poolings[0]).is_null()).height
    assert n_missing == 0, f"Missing embeddings for {n_missing} sequences!"

    print(f"Writing final parquet file: {parquet_name}")
    df_final.write_parquet(parquet_name)

    print("Testing reading the final parquet:")
    df_test = pl.read_parquet(parquet_name)
    print(df_test)
