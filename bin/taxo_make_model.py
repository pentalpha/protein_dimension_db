#!/usr/bin/env -S python3 -u
import sys
import os
import pandas as pd
import json
import numpy as np
import multiprocessing
import polars as pl
from collections import Counter
from shutil import copytree, move


def make_model_with_size(embedding_size, model_dir, taxid_tsv, vocab_json, epochs, lr):
    obo_path = "taxid.obo"
    taxallnomy_path = "taxallnomy.parquet"

    input_dim_max = 22000
    max_samples = 600000
    batch_size = 12000

    low_cpu_mode = multiprocessing.cpu_count() <= 60
    if low_cpu_mode:
        max_samples = 600000
        input_dim_max = 17000
        batch_size = 50000
    vocab_json = json.load(open(vocab_json, "r"))
    vocab = vocab_json["vocab"]
    ia_map = vocab_json["ic_map"]
    print(f"Vocab size: {len(vocab)}")
    print(f"IA map size: {len(ia_map)}")
    if len(vocab) > input_dim_max:
        print(f"Vocab size is greater than {input_dim_max}, truncating...")
        sorted_terms = vocab[:input_dim_max]
    else:
        print(f"Vocab size is less than {input_dim_max}, using all terms...")
        sorted_terms = vocab
        input_dim_max = len(vocab)

    df = pd.read_csv(taxid_tsv, sep="\t")
    # Convert semicolon string to lists

    print("Splitting raw data")
    uniprot_ids = [x for x in df["uniprot_id"].to_list()]
    raw_data = [line.split(";") for line in df["lineage"].to_list()]
    print(f"Cleaning {len(raw_data)} lines raw data")
    clean_data = [
        ";".join(sorted([f.strip() for f in sublist if f.strip()]))
        for sublist in raw_data
    ]
    print(f"Removing duplicates from {len(clean_data)} lines")
    clean_data = set(clean_data)
    print(f"Converting to list of {len(clean_data)} lines")
    clean_data = [s.split(";") for s in clean_data]
    if max_samples is not None:
        if len(clean_data) > max_samples:
            print(f"Taking {max_samples} random samples from {len(clean_data)} lines")
            random_indexes = np.random.choice(
                len(clean_data), max_samples, replace=False
            )
            clean_data = [clean_data[i] for i in random_indexes]

    print(f"Unique lineages: {len(clean_data)}")

    # Removing taxids not in vocab
    clean_data2 = set()
    vocab_in_swissprot = set()
    for taxids in clean_data:
        in_vocab = [x for x in taxids if x in vocab]
        vocab_in_swissprot.update(in_vocab)
        clean_data2.add(";".join(sorted(in_vocab)))
    # clean_data2 = [x.split(";") for x in clean_data2]
    print(f"Unique lineages (reduced vocabulary): {len(clean_data2)}")
    print(f"Unique taxids found (reduced vocabulary): {len(vocab_in_swissprot)}")
    # clean_data2 = set(clean_data2)

    # add more generic lineages
    for taxids_str in list(clean_data2):
        taxids = taxids_str.split(";")
        # iterate over successive sublists of taxids:
        for last_i in range(len(taxids) - 1):
            sublist = taxids[: last_i + 1]
            clean_data2.add(";".join(sorted(sublist)))
    clean_data2 = [s.split(";") for s in clean_data2]
    print(f"After adding more generic lineages: {len(clean_data2)}")

    values_for_freq = []
    for taxids in clean_data2:
        values_for_freq += taxids

    values_freq = Counter(values_for_freq)

    taxallnomy_df = pl.read_parquet(taxallnomy_path)
    # taxid (int), level_1(float), level_2(float), level_3(float), ...
    vocab_int = [int(x) for x in vocab]
    sub_df = taxallnomy_df.filter(pl.col("taxid").is_in(vocab_int))
    print(f"Unique taxids found (reduced vocabulary): {len(sub_df)}")

    taxid_levels = list(range(1, 42))
    taxids_by_level = {x: set() for x in taxid_levels}
    for row in sub_df.iter_rows(named=True):
        node_id = int(row["taxid"])
        for level_v in taxid_levels:
            level_name = f"level_{level_v}"
            lineage_item = int(float(row[level_name]))
            if lineage_item in vocab_int:
                taxids_by_level[level_v].add(str(lineage_item))
            # if lineage_item == node_id:
            #    break

    for level_v in taxid_levels:
        level_taxids = taxids_by_level[level_v]

        next_levels = [v for v in taxid_levels if v > level_v]
        for level_v2 in next_levels:
            taxids_by_level[level_v2] = taxids_by_level[level_v2] - level_taxids

    print("Grouped taxids by level")

    min_cluster_freq = 4
    taxa_clusters_for_sillouette = {
        "2": set(),
        "10": set(),
        "16": set(),
        "21": set(),
        "25": set(),
    }
    for level_v, items in taxids_by_level.items():
        possible_clusters = [x for x in items if values_freq[x] >= min_cluster_freq]
        if str(level_v) in taxa_clusters_for_sillouette.keys():
            taxa_clusters_for_sillouette[str(level_v)] = set(possible_clusters)
        print(
            f"Level {level_v}: {len(items)} taxids, {len(possible_clusters)} clusters"
        )

    # Initialize and Train
    from data.interpro_api.encoding import AutoEncoderWrapper

    wrapper = AutoEncoderWrapper(
        input_dim=input_dim_max,
        embedding_dim=embedding_size,
        predefined_vocab=sorted_terms,
        intermediary_len=1048,
        data_family="taxid",
    )
    os.makedirs(model_dir, exist_ok=True)
    if not wrapper.using_cpu:
        batch_size = 1200
        low_cpu_mode = False
    wrapper.fit(
        clean_data,
        ia_map,
        obo_path,
        epochs=epochs,
        lr=lr,
        directory=model_dir,
        batch_size=batch_size,
        optimize_cpu=low_cpu_mode,
        taxa_clusters_for_sillouette=taxa_clusters_for_sillouette,
    )

    # Quick test of the predict method
    wrapper = AutoEncoderWrapper.load(model_dir)

    test_sample = [[sorted_terms[0], sorted_terms[1]], [sorted_terms[3]]]
    emb = wrapper.predict(test_sample)
    print(f"Test Prediction Shape: {emb.shape}")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python taxo_make_model.py <taxid_tsv> <vocab_json> <epochs> <lr>")
        sys.exit(1)

    embedding_sizes = [8, 16, 32, 64, 128, 256, 512]
    # embedding_sizes = [8, 16]
    base_model_dir = "model_"
    taxid_tsv = sys.argv[1]
    vocab_json = sys.argv[2]
    obo_path = "taxid.obo"
    taxallnomy_path = "taxallnomy.parquet"
    if len(sys.argv) > 3:
        epochs = int(sys.argv[3])
    else:
        epochs = 12

    if len(sys.argv) > 4:
        lr = float(sys.argv[4])
    else:
        lr = 8e-4

    history_jsons = []

    for embedding_size in embedding_sizes:
        model_dir = base_model_dir + str(embedding_size)
        make_model_with_size(
            embedding_size, model_dir, taxid_tsv, vocab_json, epochs, lr
        )
        history_jsons.append(os.path.join(model_dir, "history.json"))

    lines = []

    for p in history_jsons:
        dict = json.load(open(p, "r"))

        best_score = -1
        best_line = {}

        for line in dict:
            new_score = line["score_rounded"]
            if new_score > best_score:
                best_score = new_score
                best_line = line
        best_line["model"] = os.path.dirname(p)
        lines.append(best_line)

    lines.sort(key=lambda x: x["score_rounded"])
    print(f"Results for models with embedding sizes: {embedding_sizes}")
    for line in lines:
        print(f"Model: {line['model']}, Score: {line['score_rounded']}")

    best_model_path = lines[-1]["model"]
    print(f"Moving best model to model_final")
    print(f"{best_model_path} -> model_final")
    move(best_model_path, "model_final")

    assert os.path.exists("model_final")

    print("Done")
