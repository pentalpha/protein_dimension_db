#!/usr/bin/env -S python3 -u
import sys
import time

import pandas as pd
import numpy as np
import polars as pl

from data.interpro_api.encoding import AutoEncoderWrapper

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(
            "Usage: python interpro_encode.py <model_dir> <id_col> <interpros_col> <input_df>"
        )
        sys.exit(1)

    model_dir = sys.argv[1]
    id_col = sys.argv[2]
    interpros_col = sys.argv[3]
    annots_df_path = sys.argv[4]

    fmts = ["tsv", "csv", "parquet"]

    end_fmt = annots_df_path.split(".")[-1]
    assert end_fmt in fmts, f"Unsupported format {end_fmt}"

    if end_fmt == "tsv":
        read_df = pd.read_csv(annots_df_path, sep="\t")
    elif end_fmt == "csv":
        read_df = pd.read_csv(annots_df_path)
    elif end_fmt == "parquet":
        read_df = pl.read_parquet(annots_df_path)
    else:
        raise ValueError(f"Unsupported format {end_fmt}")

    # to_list from pandas and tolist from polars
    if end_fmt in ["tsv", "csv"]:
        ids = read_df[id_col].to_list()
        terms = read_df[interpros_col].to_list()
    else:
        ids = read_df[id_col].tolist()
        terms = read_df[interpros_col].tolist()

    print(f"Number of IDs: {len(ids)}")
    print(f"Number of taxid lists: {len(terms)}")
    print(f"IDs: {ids[:10]}")
    print(f"Taxids: {terms[:10]}")

    terms = [line if type(line) == str else "" for line in terms]

    raw_data = [";".join([tx.strip() for tx in line.split(";")]) for line in terms]

    original_index_to_nr_index = {}
    non_redundant = []
    non_redundant_set = set()
    for idx, r in enumerate(raw_data):
        if r not in non_redundant_set:
            non_redundant.append(r)
            non_redundant_set.add(r)
        original_index_to_nr_index[idx] = len(non_redundant) - 1

    print(f"Original length: {len(raw_data)}")
    print(f"Non-redundant length: {len(non_redundant)}")
    print(f"Mapping length: {len(original_index_to_nr_index)}")

    clean_data = [line.split(";") for line in non_redundant]
    print(f"Unique lineages: {len(clean_data)}")

    print("Loading model")
    wrapper = AutoEncoderWrapper.load(model_dir)
    encoded = []

    batch_len = 50
    last_prot_per_sec = 0.0
    last_batch_len = 50
    keep_increasing = True
    max_batch_len = 10000
    batch_increase_factor = 1.33333
    next_batch_start_idx = 0
    round = 0

    while next_batch_start_idx < len(clean_data):
        batch = clean_data[next_batch_start_idx : next_batch_start_idx + batch_len]
        processing_start = time.time()
        encoded += wrapper.predict(batch).tolist()
        proc_secs = time.time() - processing_start
        proteins_per_sec = len(batch) / proc_secs
        next_batch_start_idx += batch_len

        if keep_increasing:
            if proteins_per_sec > last_prot_per_sec:
                new_batch_len = int(batch_len * batch_increase_factor)
                if new_batch_len > max_batch_len:
                    new_batch_len = max_batch_len
                print(
                    f"Speed increase from {last_prot_per_sec} to {proteins_per_sec}.\n"
                    f"New batch length: {new_batch_len}"
                )
                batch_len = new_batch_len
            else:
                keep_increasing = False
                print(
                    f"Speed decrease from {last_prot_per_sec} to {proteins_per_sec}.\n"
                    f"Back to batch length: {last_batch_len}"
                )
                batch_len = last_batch_len

        last_prot_per_sec = proteins_per_sec
        last_batch_len = batch_len
        perc = next_batch_start_idx / len(clean_data) * 100
        round += 1
        print(f"Round {round}: {perc:.2f}%")

    embs_list = [np.nan] * len(ids)
    for original_idx, nr_idx in original_index_to_nr_index.items():
        embs_list[original_idx] = np.array(encoded[nr_idx])

    print(embs_list[0])
    print(type(embs_list[0]))
    print(embs_list[0].shape)
    print(type(embs_list))

    print(ids[0])
    print(terms[0])

    print(f"Lens: {len(ids), len(embs_list), len(terms)}")

    new_pqt_df = pl.DataFrame(
        {
            "id": ids,
            "emb": embs_list,
            "interpros": terms,
        }
    )
    new_pqt_df.write_parquet("emb.interpro_autoencoded.parquet")
