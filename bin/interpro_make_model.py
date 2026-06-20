#!/usr/bin/env -S python3 -u
import sys
import os
import pandas as pd
import json
import numpy as np
import multiprocessing
from shutil import copytree, move

import matplotlib.pyplot as plt
from collections import Counter

from plotting_lib.histograms import plot_term_freq_histogram
from plotting_lib.video_maker import generate_pca_video
from bioinfo_utils.clustering import (
    find_intepro_clusterings,
)


def make_model_with_size(
    model_type, embedding_size, model_dir, interproscan_tsv, vocab_json, epochs, lr
):
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

    df = pd.read_csv(interproscan_tsv, sep="\t", header=None)
    # Convert semicolon string to lists

    print("Splitting raw data")
    raw_data = [str(row).split(";") for row in df[1].values]
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
    print(f"Unique protein annotations: {len(clean_data)}")

    values_for_freq = []
    for taxids in clean_data:
        values_for_freq += taxids

    print("Counting term frequencies")
    values_freq = Counter(values_for_freq)

    plot_term_freq_histogram(values_for_freq, model_dir)

    print("Creating term clusters")
    clusterings = find_intepro_clusterings(values_freq, clean_data)

    print("Training model")

    # Initialize and Train
    if model_type == "autoencoder":
        from data.interpro_api.encoding import AutoEncoderWrapper

        wrapper = AutoEncoderWrapper(
            input_dim=input_dim_max,
            embedding_dim=embedding_size,
            predefined_vocab=sorted_terms,
            data_family="interpro",
        )
        os.makedirs(model_dir, exist_ok=True)
        if not wrapper.using_cpu:
            batch_size = 5000
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
            taxa_clusters_for_sillouette=clusterings,
            eval_perc=0.5,
        )

        generate_pca_video(model_dir, fps=2.5)
    elif model_type == "onehot":
        from data.interpro_api.onehot_encoder import OneHotEncoder

        wrapper = OneHotEncoder(
            max_families=embedding_size, predefined_vocab=sorted_terms
        )
        os.makedirs(model_dir, exist_ok=True)
        wrapper.fit(clean_data)
        # Save the result
        wrapper.save(model_dir)
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    # Quick test of the predict method
    if model_type == "autoencoder":
        from data.interpro_api.encoding import AutoEncoderWrapper

        wrapper = AutoEncoderWrapper.load(model_dir)
    elif model_type == "onehot":
        from data.interpro_api.onehot_encoder import OneHotEncoder

        wrapper = OneHotEncoder.load(model_dir)
    else:
        raise ValueError(f"Unknown model type: {model_type}")
    test_sample = [[sorted_terms[0], sorted_terms[1]], [sorted_terms[3]]]
    emb = wrapper.predict(test_sample)
    print(f"Test Prediction Shape: {emb.shape}")


if __name__ == "__main__":
    # singularity run --nv ~/repos/protein_dimension_db/singularity/sif/torch_frieren.sif python src/interpro_make_model.py autoencoder 32 model_v5_32 concatenated.tsv interpro_vocab_ia.json 90 8e-3
    # concatenated.tsv, interpro_vocab_ia.json, interpro.obo
    if len(sys.argv) < 6:
        print(
            "Usage: python interpro_make_model.py [autoencoder|onehot] <interproscan_tsv> <vocab_json> <epochs> <lr>"
        )
        sys.exit(1)

    embedding_sizes = [32, 8, 16, 64, 128]

    model_type = sys.argv[1]
    assert model_type in ["autoencoder", "onehot"]
    interproscan_tsv = sys.argv[2]
    vocab_json = sys.argv[3]
    obo_path = "interpro.obo"
    if len(sys.argv) > 4:
        epochs = int(sys.argv[4])
    else:
        epochs = 12

    if len(sys.argv) > 5:
        lr = float(sys.argv[5])
    else:
        lr = 8e-4

    base_model_dir = "model_"

    # make_model_with_size(model_type, embedding_size, model_dir, interproscan_tsv, vocab_json, epochs, lr)

    history_jsons = []

    for embedding_size in embedding_sizes:
        model_dir = base_model_dir + str(embedding_size)
        make_model_with_size(
            model_type,
            embedding_size,
            model_dir,
            interproscan_tsv,
            vocab_json,
            epochs,
            lr,
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
