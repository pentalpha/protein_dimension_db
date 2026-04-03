#!/usr/bin/env python
# Usage: python encode_proteins.py <encoder_type:interpro_autoencoder|interpro_onehot>
#   <model_dir> <output_parquet> <protein_ids_sorted> <input_files...>

import sys
from tqdm import tqdm
import polars as pl

if __name__ == "__main__":
    encoder_type = sys.argv[1]
    model_dir = sys.argv[2]
    output_parquet = sys.argv[3]
    protein_ids_sorted = sys.argv[4]
    input_files = sys.argv[5:]

    if encoder_type == "interpro_autoencoder":
        from src.data.interpro_api.encoding import AutoEncoderWrapper
        from data.interpro_api.parsing import parse_interpro_raw as prot_class_parser
        from data.interpro_api.parsing import (
            parsed_interpro_to_data as prot_class_to_data,
        )

        encoder = AutoEncoderWrapper.load(model_dir)

    elif encoder_type == "interpro_onehot":
        from src.data.interpro_api.onehot_encoder import OneHotEncoder
        from data.interpro_api.parsing import parse_interpro_raw as prot_class_parser
        from data.interpro_api.parsing import (
            parsed_interpro_to_data as prot_class_to_data,
        )

        encoder = OneHotEncoder.load(model_dir)
    else:
        raise ValueError(f"Unknown encoder type: {encoder_type}")

    protein_ids_sorted = open(protein_ids_sorted, "r").read().strip().split("\n")
    parsed_input_path = input_files[0] + ".tsv"
    prot_class_parser(",".join(input_files), parsed_input_path)
    ids_list, clean_data = prot_class_to_data(parsed_input_path)

    original_order_map = {pid: i for i, pid in enumerate(protein_ids_sorted)}
    sorted_data = [encoder.predict([]) for _ in range(len(protein_ids_sorted))]
    for i, pid in tqdm(enumerate(ids_list)):
        classes = clean_data[i]
        original_pos = original_order_map[pid]
        sorted_data[original_pos] = encoder.predict([classes])[0]

    df = pl.DataFrame({"id": protein_ids_sorted, "emb": sorted_data})
    df.write_parquet(output_parquet)
