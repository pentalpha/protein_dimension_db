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
            "Usage: python interpro_autoencoder.py <embedding_size> <model_dir> <raw_interpro_output> <vocab_json>"
        )
        sys.exit(1)

    model_type = "autoencoder"
    embedding_size = int(sys.argv[1])  # current best: 64
    model_dir = sys.argv[2]
    taxid_tsv = sys.argv[3]
    vocab_json = sys.argv[4]
    obo_path = "taxid.obo"
    if len(sys.argv) > 5:
        epochs = int(sys.argv[5])
    else:
        epochs = 12

    if len(sys.argv) > 6:
        lr = float(sys.argv[6])
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
    clean_data2 = [x.split(";") for x in clean_data2]
    print(f"Unique lineages (reduced vocabulary): {len(clean_data2)}")
    print(f"Unique taxids found (reduced vocabulary): {len(vocab_in_swissprot)}")

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
    )

    # Quick test of the predict method
    wrapper = AutoEncoderWrapper.load(model_dir)

    test_sample = [[sorted_terms[0], sorted_terms[1]], [sorted_terms[3]]]
    emb = wrapper.predict(test_sample)
    print(f"Test Prediction Shape: {emb.shape}")
