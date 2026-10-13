import os
import sys
import glob
import math
import json

import polars as pl

# Configuration
CHUNK_SIZE = 80_000
TARGET_POOLINGS = ["max", "mean", "std", "parti"]

if len(sys.argv) >= 2:
    DATA_DIR = sys.argv[1]
else:
    DATA_DIR = "."

if __name__ == "__main__":
    # 1. Find all embedding parquets in the directory
    search_pattern = os.path.join(DATA_DIR, "emb.*.parquet")
    files = glob.glob(search_pattern)
    
    # 2. Group files by their model name
    # Example: emb.ankh3_large_s2s_std.parquet -> model: "ankh3_large_s2s", pooling: "std"
    models = {}
    for filepath in files:
        basename = os.path.basename(filepath)
        name = basename.replace("emb.", "").replace(".parquet", "")
        
        for pooling in TARGET_POOLINGS:
            if name.endswith(f"_{pooling}"):
                model_name = name[:-len(f"_{pooling}")]
                if model_name not in models:
                    models[model_name] = {}
                models[model_name][pooling] = filepath
                break
    
    if "ankh3_large_nlu" in models:
        models["ankh3_large"] = models["ankh3_large_nlu"]
        del models["ankh3_large_nlu"]
    
    if "ankh3_large_s2s" in models:
        del models["ankh3_large_s2s"]
                
    # 3. Filter for models that have all 4 target poolings
    complete_models = {m: files for m, files in models.items() if len(files) == len(TARGET_POOLINGS)}
    print(f"Found {len(complete_models)} models with all 4 poolings: {list(complete_models.keys())}\n")

    print(json.dumps(complete_models, indent=4))

    # 4. Process each model
    for model_name, pooling_files in complete_models.items():
        print(f"--- Processing model: {model_name} ---")
        
        # Create output directory for this specific model
        out_dir = os.path.join(DATA_DIR, model_name)
        os.makedirs(out_dir, exist_ok=True)
        
        # Use 'max' as the base to determine the total number of rows
        base_file = pooling_files["max"]
        
        # scan_parquet reads metadata without loading data into memory
        print("columns:", pl.scan_parquet(base_file).collect_schema().names)
        
        total_rows = pl.scan_parquet(base_file).select(pl.len()).collect().item()
        total_chunks = math.ceil(total_rows / CHUNK_SIZE)
        
        print(f"Total rows: {total_rows:,} | Chunks to create: {total_chunks}")
        
        for chunk_idx in range(total_chunks):
            offset = chunk_idx * CHUNK_SIZE
            chunk_filename = f"train-{chunk_idx:05d}-of-{total_chunks:05d}.parquet"
            out_path = os.path.join(out_dir, chunk_filename)

            if os.path.exists(out_path):
                print(f"  [Skipping] {chunk_filename} already exists.")
                continue

            print(f"  [Building] {chunk_filename} (rows {offset} to {offset + CHUNK_SIZE})...")

            chunk_dfs = {}
            for pooling, fpath in pooling_files.items():
                lf = pl.scan_parquet(fpath).slice(offset, CHUNK_SIZE)
                if pooling != "max":
                    lf = lf.drop("seq")   # keep only id + the embedding column
                chunk_dfs[pooling] = lf

            merged_lf = chunk_dfs["max"]
            for pooling in ["mean", "std", "parti"]:
                merged_lf = merged_lf.join(chunk_dfs[pooling], on="id", how="left")

            merged_lf.sink_parquet(out_path)
            
        print(f"Finished {model_name}.\n")