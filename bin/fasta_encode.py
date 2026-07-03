#!/usr/bin/env python3
import sys
import os
import time

from transformers import AutoModel
import polars as pl
import torch

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
    "default": 12000,
}

nvidiah100_max_tokens = {
    "Synthyra/ANKH_base": 14000,
    "Synthyra/ANKH_large": 14000,
    "default": 12000,
}

default_max_tokens = nvidiah100_max_tokens

fasta_path = sys.argv[1]
cache_path = sys.argv[2]
model_name = sys.argv[3]  # Synthyra/ANKH_base
parquet_name = sys.argv[4]

MAX_TOKENS_PER_BATCH = default_max_tokens.get(model_name, default_max_tokens["default"])
N_PARTS = 1337
print(f"Using model: {model_name}")
print(f"MAX_TOKENS_PER_BATCH: {MAX_TOKENS_PER_BATCH}, N_PARTS: {N_PARTS}")

fasta_part_paths = fasta_equal_split_by_len(fasta_path, N_PARTS)

ids = []
seqs = []

model_safe_name = model_name.replace("/", "__").lower()
save_path = cache_path + "/" + model_safe_name + ".pt"
model = AutoModel.from_pretrained(model_name, trust_remote_code=True).to("cuda")
poolings = ["mean", "std", "max"]
embs = {p: [] for p in poolings}
n_poolings = len(poolings)

fasta_i = 0
time_start = time.time()
durations = []

for fasta_part_path in reversed(fasta_part_paths):
    print(f"Processing fasta part {fasta_i} out of {len(fasta_part_paths)}")
    mean_duration = sum(durations) / len(durations) if fasta_i > 0 else 0
    print(f"Mean duration: {mean_duration:.2f}s")
    remaining = (len(fasta_part_paths) - fasta_i) * mean_duration
    print(f"Remaining: {remaining / 60:.2f}min")
    fasta_i += 1

    part_time_start = time.time()
    fasta_content = read_fasta(fasta_part_path)
    max_len = max(len(seq) for _, seq in fasta_content)
    if max_len > MAX_TOKENS_PER_BATCH:
        print(
            f"Max len in chunk ({max_len}) is greater than "
            f"MAX_TOKENS_PER_BATCH ({MAX_TOKENS_PER_BATCH})"
        )
        # truncate to max_len for this chunk if it exceeds the limit, however this will truncate the protein
        max_len = MAX_TOKENS_PER_BATCH
        # quit(1)
    fasta_content_shortened = []
    for seq_id, seq in fasta_content:
        if len(seq) > max_len:
            # first max_len tokens
            seq = seq[:max_len]
        fasta_content_shortened.append((seq_id, seq))
    fasta_content = fasta_content_shortened
    optimal_batch_len_float = MAX_TOKENS_PER_BATCH / max_len
    processing_batch_size = max(int(optimal_batch_len_float), 1)
    print(f"Setting max len to actual max len in the dataset: {max_len}")
    print(f"Processing batch size: {processing_batch_size} ({optimal_batch_len_float})")

    embeddings = model.embed_dataset(
        fasta_path=fasta_part_path,
        batch_size=processing_batch_size,
        pooling_types=poolings,
        max_len=max_len,
        save=True,
        save_path=save_path,
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

    for seq_id, seq in fasta_content:
        emb = embeddings[seq]
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
        seqs.append(seq)
        ids.append(seq_id)

    part_duration = time.time() - part_time_start
    durations.append(part_duration)
    print(f"Part {fasta_i} took {part_duration} seconds")

df = pl.DataFrame({"id": ids, "seq": seqs, **embs})
df.write_parquet(f"{parquet_name}")

for fasta_part_path in fasta_part_paths:
    os.remove(fasta_part_path)

df = pl.read_parquet(f"{parquet_name}")
print(df)
