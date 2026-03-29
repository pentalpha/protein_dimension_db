import sys
import os
import pandas as pd
from data.interpro_api.onehot_encoder import OneHotEncoder
from data.interpro_api.parsing import parse_interpro_raw

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: python interpro_onehot.py <embedding_size> <model_dir> <raw_interpro_output>"
        )
        sys.exit(1)

    embedding_size = int(sys.argv[1])  # current best: 800
    model_dir = sys.argv[2]
    raw_interproscan_files = sys.argv[3:]

    input_file = raw_interproscan_files[0]

    interpro_parsed_tsv = "interpro_clfs_parsed.tsv"
    parse_interpro_raw(input_file, interpro_parsed_tsv)

    # Load data from TSV
    print(f"Loading data from {input_file} -> {interpro_parsed_tsv}...")
    df = pd.read_csv(interpro_parsed_tsv, sep="\t", header=None)
    # Convert semicolon string to lists
    raw_data = [str(row).split(";") for row in df[1].values]
    clean_data = [[f.strip() for f in sublist if f.strip()] for sublist in raw_data]

    # Initialize and Train
    wrapper = OneHotEncoder(max_families=embedding_size)
    wrapper.fit(clean_data)

    # Save the result
    os.makedirs(model_dir, exist_ok=True)
    wrapper.save(model_dir)

    # Quick test of the predict method
    test_sample = [["IPR034100", "IPR001163"], ["IPR013965"]]
    emb = wrapper.predict(test_sample)
    print(f"Test Prediction Shape: {emb.shape}")
