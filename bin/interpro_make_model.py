#!/usr/bin/env python
import sys
import os
import pandas as pd
import json

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
    vocab = json.load(open(vocab_json, 'r'))
    sorted_terms = vocab[:embedding_size]
    
    df = pd.read_csv(interproscan_tsv, sep="\t", header=None)
    # Convert semicolon string to lists
    raw_data = [str(row).split(";") for row in df[1].values]
    clean_data = [[f.strip() for f in sublist if f.strip()] for sublist in raw_data]

    # Initialize and Train
    if model_type == "autoencoder":
        from data.interpro_api.encoding import AutoEncoderWrapper

        wrapper = AutoEncoderWrapper(input_dim=22000, embedding_dim=embedding_size, predefined_vocab=sorted_terms)
        wrapper.fit(clean_data, epochs=12)
    elif model_type == "onehot":
        from data.interpro_api.onehot_encoder import OneHotEncoder

        wrapper = OneHotEncoder(max_families=embedding_size, predefined_vocab=sorted_terms)
        wrapper.fit(clean_data)
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    # Save the result
    os.makedirs(model_dir, exist_ok=True)
    wrapper.save(model_dir)

    # Quick test of the predict method
    test_sample = [[sorted_terms[0], sorted_terms[1]], [sorted_terms[3]]]
    emb = wrapper.predict(test_sample)
    print(f"Test Prediction Shape: {emb.shape}")
