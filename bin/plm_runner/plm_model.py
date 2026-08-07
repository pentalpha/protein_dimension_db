from typing import List, Tuple
import sys
import os
import time

import numpy as np
import torch
from transformers import (
    T5EncoderModel,
    #T5ForConditionalGeneration,
    AutoTokenizer,
    #TFT5EncoderModel,
    #TFT5ForConditionalGeneration,
    T5Tokenizer,
)
from tqdm import tqdm
import polars as pl

from plm_runner.poolings import POOLERS
from plm_runner.caching import batch_by_tokens, EmbCache, read_fasta


"""
Class to load any of the following PLMs from huggingface: ESM, ANKH, DPLM, Profluent E1
Loading, embedding, batched embedding by token count, on-the-fly pooling
"""

def get_device(use_gpu: bool = True) -> torch.device:
    if torch.cuda.is_available() and use_gpu:
        device = torch.device("cuda:0")
    else:
        device = torch.device("cpu")
    return device

def get_vram_gb() -> int:
    """
    Checks for a CUDA device and returns its total VRAM in GB as an integer.
    Defaults to 4 if no CUDA device is available or if an error occurs.
    """
    try:
        if torch.cuda.is_available():
            # Get total memory of the current default device (GPU 0) in bytes
            total_bytes = torch.cuda.get_device_properties(0).total_memory
            
            # Convert bytes to Gigabytes (1 GB = 1024^3 bytes)
            total_gb = total_bytes / (1024 ** 3)
            
            # Round to the nearest whole number to get standard sizes (e.g., 7.93 GB -> 8)
            return int(round(total_gb))
        else:
            return 4
            
    except Exception:
        # Fallback to 4 for any unexpected errors
        return 4

max_tokens_by_model = [
    {
        "VRAM": 6,
        "Synthyra/ANKH_base": 400,
        "ElnaggarLab/ankh-base": 9600,
        "Synthyra/ANKH_large": 680,
        "Synthyra/ANKH2_large": 680,
        "Synthyra/ANKH3_large": 680,
        "Synthyra/ANKH3_xl": 320,
        "default": 600,
    },
    {
        "VRAM": 16,
        "Synthyra/ANKH_base": 1800,
        "Synthyra/ANKH_large": 14000,
        "Synthyra/ANKH2_large": 14000,
        "Synthyra/ANKH3_large": 14000,
        "Synthyra/ANKH3_xl": 6000,
        "default": 1800,
    }
]

def get_max_tokens_dict():
    total_vram = get_vram_gb()
    print(f"Total VRAM: {total_vram} GB")
    best_match = max_tokens_by_model[0]
    for d in max_tokens_by_model:
        if d["VRAM"] <= total_vram and d["VRAM"] > best_match["VRAM"]:
            best_match = d
    return best_match

def define_plm_class(model_name: str):
    n = model_name.lower()
    if "esm" in n:
        return "ESM"
    elif "ankh" in n:
        return "ANKH"
    elif "dplm" in n:
        return "DPLM"
    elif "profluent" in n:
        return "PROFLUENT"
    else:
        raise ValueError(f"Unknown model name: {model_name}")

AVAILABLE_MODELS = {
    "ElnaggarLab/ankh-base": {"type": "ANKH"},
    "ElnaggarLab/ankh-large": {"type": "ANKH"},
    "ElnaggarLab/ankh3-large": {"type": "ANKH"},
    "ElnaggarLab/ankh3-xl": {"type": "ANKH"}
}

class PLMModel():
    def __init__(self, model_name: str, cache_path: str):
        self.model_name = model_name
        self.token = None
        if 'HFTOKEN' in os.environ:
            print(f"Loading key", os.environ['HFTOKEN'][12:])
            self.token = os.environ['HFTOKEN']
        if model_name not in AVAILABLE_MODELS:
            raise ValueError(f"Invalid model name: {model_name}. Valid models are: {list(AVAILABLE_MODELS.keys())}")
        self.model_type = define_plm_class(model_name)
        self.device = get_device()
        self.cache = EmbCache(cache_path, model_name)
        self.max_tokens_dict = get_max_tokens_dict()
        if model_name in self.max_tokens_dict:
            self.max_tokens = self.max_tokens_dict[model_name]
        else:
            self.max_tokens = self.max_tokens_dict["default"]

        print(f"Using max_tokens: {self.max_tokens} for model: {model_name}")

    def extract(self, seqs: List[str]):
        raise NotImplementedError("Please use the specific model class for embedding.")
        
    def embed(self, seqs, 
            max_tokens_per_batch = None, 
            poolings: List[str] = list(POOLERS.keys()),
            tqdm_bar=None):
        if max_tokens_per_batch is None:
            max_tokens_per_batch = self.max_tokens
        batchs = batch_by_tokens(seqs, max_tokens_per_batch, use_padding=True)
        embeddings = {p: [] for p in poolings}
        print(f"Processing {len(seqs)} sequences in {len(batchs)} batches of max {max_tokens_per_batch} tokens each.")
        if tqdm_bar is not None:
            to_iter = batchs
        else:
            to_iter = tqdm(batchs, desc="Embedding sequences")
        for batch in to_iter:
            lens = [len(seq) for seq in batch]
            #print(lens)
            full_batch = self.extract(batch)
            
            for p in poolings:
                new_pooled = POOLERS[p](full_batch)
                for emb_fixed_size in new_pooled:
                    embeddings[p].append(emb_fixed_size)
                last_emb = new_pooled[-1]
                #print(f"Embeddings desc {p}: shape={last_emb.shape}, dtype={last_emb.dtype}")
            if tqdm_bar is not None:
                tqdm_bar.update(1)
        return embeddings

    def write_all_embeddings(self, fasta_path, parquet_path, poolings: List[str]):
        print("Reading original fasta file")
        fasta_content = read_fasta(fasta_path)
        ids = [seq_id for seq_id, _ in fasta_content]
        seqs = [seq for _, seq in fasta_content]

        # 1. Create a Polars DataFrame of the exact original FASTA order
        df_fasta = pl.DataFrame({"id": ids, "seq": seqs}).with_row_index("original_order")

        all_cache_files = self.cache.list_caches()
        valid_pqs = [pq_path for pq_path, txt_path in all_cache_files]

        print("Lazy loading and merging embeddings via Polars...")

        # We must load "seq" to perform the join, plus whatever poolings were requested
        columns_to_load = ["seq"] + poolings

        # 2. Lazily scan all parquet chunks simultaneously.
        # .select() pushes the column filter down to the parquet reader, 
        # meaning unrequested poolings are NEVER loaded into memory.
        df_embs = pl.scan_parquet(valid_pqs).select(columns_to_load)

        # Deduplicate embeddings just in case a failed job left duplicate sequences
        df_embs = df_embs.unique(subset=["seq"], keep="first")

        # 3. Build the computation graph: Join the FASTA frame with the embeddings, then sort
        df_final_lazy = (
            df_fasta.lazy()
            .join(df_embs, on="seq", how="left")
            .sort("original_order")
            .drop("original_order")
        )

        print("Executing join, sort, and collecting into memory...")
        # 4. .collect() executes the highly optimized Rust code.
        df_final = df_final_lazy.collect()

        print("Checking for missing embeddings...")
        # If a sequence from the FASTA wasn't found in the caches, the left-join will leave nulls
        n_missing = df_final.filter(pl.col(poolings[0]).is_null()).height
        assert n_missing == 0, f"Missing embeddings for {n_missing} sequences!"

        print(f"Writing final parquet file: {parquet_path}")
        df_final.write_parquet(parquet_path)

        print("Testing reading the final parquet:")
        df_test = pl.read_parquet(parquet_path)
        print(df_test)
    
    def embed_saving_progress(self, fasta_path: str, parquet_path, poolings: List[str] = list(POOLERS.keys()), max_tokens_per_batch = None):
        if max_tokens_per_batch is None:
            max_tokens_per_batch = self.max_tokens
        caching_batch_len = max_tokens_per_batch * 100
        not_embedded_seqs = self.cache.list_non_embedded(fasta_path)
        macro_batchs = batch_by_tokens(not_embedded_seqs, caching_batch_len, use_padding=False)
        small_batchs = [batch_by_tokens(seqs, max_tokens_per_batch, use_padding=True) 
            for seqs in macro_batchs]
        n_micro_batchs = [len(b2) for b2 in small_batchs]
        total_micro_batchs = sum(n_micro_batchs)
        bar = tqdm(total=total_micro_batchs, desc="Embedding sequences")
        for i, macro_batch in enumerate(macro_batchs):
            time_start = time.time()
            new_embeddings = self.embed(macro_batch, poolings=poolings, tqdm_bar=bar)
            next_parquet, next_txt = self.cache.next_cache_name()
            df_dict = {"seq": macro_batch}
            for p_name in poolings:
                if p_name == "full":
                    full_list = [[ [float(e) for e in residue_emb] for residue_emb in emb ] 
                        for emb in new_embeddings[p_name]]
                    df_dict[p_name] = full_list
                else:
                    df_dict[p_name] = np.asarray(new_embeddings[p_name])
            '''print("Schema:")
            for key, vals_list in df_dict.items():
                print(f"{key}: {len(vals_list)} ({type(vals_list[0])})")'''
            df = pl.DataFrame(df_dict)
            df.write_parquet(next_parquet)
            with open(next_txt, "w") as f:
                for seq in macro_batch:
                    f.write(f"{seq}\n")
            print(f"Saved cache to {next_parquet} and {next_txt}")
            duration = time.time() - time_start
            print(f"Duration: {duration}")

        self.write_all_embeddings(fasta_path, parquet_path, poolings)