#!/usr/bin/env -S python3 -u
import sys
import os
import pandas as pd
import json
import numpy as np
import multiprocessing

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(
            "Usage: python interpro_autoencoder.py <model_type> <embedding_size> <model_dir> <raw_interpro_output> <vocab_json>"
        )
        sys.exit(1)

    model_type = sys.argv[1]
    assert model_type in ["autoencoder", "onehot"]
    embedding_size = int(sys.argv[2])  # current best: 800
    model_dir = sys.argv[3]
    interproscan_tsv = sys.argv[4]
    vocab_json = sys.argv[5]
    if len(sys.argv) > 6:
        epochs = int(sys.argv[6])
    else:
        epochs = 12

    if len(sys.argv) > 7:
        lr = float(sys.argv[7])
    else:
        lr = 8e-4

    input_dim_max = 22000
    max_samples = 600000
    batch_size = 12000

    low_cpu_mode = multiprocessing.cpu_count() <= 60
    if low_cpu_mode:
        max_samples = 600000
        input_dim_max = 17000
        batch_size = 50000

    vocab = json.load(open(vocab_json, "r"))["vocab"]
    print(f"Vocab size: {len(vocab)}")
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
    print(f"Proteins: {len(clean_data)}")

    # Initialize and Train
    if model_type == "autoencoder":
        from data.interpro_api.encoding import AutoEncoderWrapper

        wrapper = AutoEncoderWrapper(
            input_dim=input_dim_max,
            embedding_dim=embedding_size,
            predefined_vocab=sorted_terms,
        )
        os.makedirs(model_dir, exist_ok=True)
        if not wrapper.using_cpu:
            batch_size = 5000
            low_cpu_mode = False
        wrapper.fit(
            clean_data,
            epochs=epochs,
            lr=lr,
            directory=model_dir,
            batch_size=batch_size,
            optimize_cpu=low_cpu_mode,
        )
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
