#!/usr/bin/env python
import sys
import os
import pandas as pd

if __name__ == "__main__":
    model_type = sys.argv[1]
    assert model_type in ["autoencoder", "onehot"]
    embedding_size = int(sys.argv[2])  # current best: 800
    model_dir = sys.argv[3]
    interproscan_tsvs = sys.argv[4:]

    if len(sys.argv) < 5:
        print(
            "Usage: python interpro_autoencoder.py <model_type> <embedding_size> <model_dir> <raw_interpro_output>"
        )
        sys.exit(1)

    # input_file = raw_interproscan_files[0]

    # interpro_parsed_tsv = "interpro_clfs_parsed.tsv"
    # parse_interpro_raw(input_file, interpro_parsed_tsv)

    # Load data from TSV
    dfs = [
        pd.read_csv(interpro_parsed_tsv, sep="\t", header=None)
        for interpro_parsed_tsv in interproscan_tsvs
    ]
    df = pd.concat(dfs, ignore_index=True)
    # Convert semicolon string to lists
    raw_data = [str(row).split(";") for row in df[1].values]
    clean_data = [[f.strip() for f in sublist if f.strip()] for sublist in raw_data]

    # Initialize and Train
    if model_type == "autoencoder":
        from data.interpro_api.encoding import AutoEncoderWrapper

        wrapper = AutoEncoderWrapper(input_dim=22000, embedding_dim=embedding_size)
        wrapper.fit(clean_data, epochs=12)
    elif model_type == "onehot":
        from data.interpro_api.onehot_encoder import OneHotEncoder

        wrapper = OneHotEncoder(max_families=embedding_size)
        wrapper.fit(clean_data)
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    # Save the result
    os.makedirs(model_dir, exist_ok=True)
    wrapper.save(model_dir)

    # Quick test of the predict method
    test_sample = [["IPR034100", "IPR001163"], ["IPR013965"]]
    emb = wrapper.predict(test_sample)
    print(f"Test Prediction Shape: {emb.shape}")
