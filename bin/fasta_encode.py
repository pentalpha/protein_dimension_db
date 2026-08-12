#!/usr/bin/env python3
import sys
import os
import time

import polars as pl
from tqdm import tqdm
import numpy as np

from plm_runner.plm_model_lib import plm_master_loader

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    cache_path = sys.argv[2]
    model_name = sys.argv[3]  # Synthyra/ANKH_base
    parquet_name = sys.argv[4]
    poolings_list = ','.split(sys.argv[5]) if len(sys.argv) > 5 else None
    if poolings_list is None:
        poolings_list = ["mean", "max", "softmax", "std", "norm",
                "k4p_max", "k8p_max", "k16p_max", "k32p_max", "parti","full"]
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)
    
    print(f"Using model: {model_name}")

    model = plm_master_loader(model_name, cache_path)
    model.embed_saving_progress(fasta_path, parquet_name.replace('.parquet', ''), 
        poolings_list)