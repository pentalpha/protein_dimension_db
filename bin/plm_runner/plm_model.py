from typing import List, Tuple
import sys
import os
import time
import gc

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
        "VRAM": 4,
        "default": 600,
    },
    {
        "VRAM": 6,
        "Profluent-Bio/E1-150m": 10000,
        "Profluent-Bio/E1-300m": 8000,
        "Profluent-Bio/E1-600m": 1800,
        "Synthyra/ANKH_base": 8000,
        "ElnaggarLab/ankh-base": 8000,
        "Synthyra/ANKH_large": 3000,
        "Synthyra/ANKH2_large": 3000,
        "Synthyra/ANKH3_large": 3000,
        "Synthyra/ANKH3_xl": 3000,
        "default": 2500,
    },
    {
        "VRAM": 16,
        "Profluent-Bio/E1-150m": 10000,
        "Profluent-Bio/E1-300m": 8000,
        "Profluent-Bio/E1-600m": 6000,
        "Synthyra/ANKH_base": 8000,
        "ElnaggarLab/ankh-base": 8000,
        "Synthyra/ANKH_large": 6000,
        "Synthyra/ANKH2_large": 6000,
        "Synthyra/ANKH3_large": 6000,
        "Synthyra/ANKH3_xl": 6000,
        "default": 2500,
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
    "Synthyra/ANKH_base": {"type": "ANKH"},
    "Synthyra/ANKH_large": {"type": "ANKH"},
    "Synthyra/ANKH2_large": {"type": "ANKH"},
    "Synthyra/ANKH3_large": {"type": "ANKH"},
    "Synthyra/ANKH3_xl": {"type": "ANKH"},
    "ElnaggarLab/ankh-base": {"type": "ANKH"},
    "ElnaggarLab/ankh-large": {"type": "ANKH"},
    "ElnaggarLab/ankh3-large": {"type": "ANKH"},
    "ElnaggarLab/ankh3-xl": {"type": "ANKH"},
    "Profluent-Bio/E1-150m": {"type": "PROFLUENT"},
    "Profluent-Bio/E1-300m": {"type": "PROFLUENT"},
    "Profluent-Bio/E1-600m": {"type": "PROFLUENT"},
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
            try:
                full_batch, full_batch_attn = self.extract(batch)
            except torch.OutOfMemoryError as err:
                print(err)
                print("Switching to per_sec because of OOM.")
                full_batch = []
                full_batch_attn = []
                for seq in batch:
                    new_batch, new_attn = self.extract([seq])
                    full_batch.append(new_batch)
                    full_batch_attn.append(new_attn)
            
            for p in poolings:
                pool_func = POOLERS[p]
                if p == 'parti':
                    new_pooled = pool_func(full_batch, full_batch_attn)
                else:
                    new_pooled = pool_func(full_batch)
                for emb_fixed_size in new_pooled:
                    embeddings[p].append(emb_fixed_size)
                last_emb = new_pooled[-1]
                #print(f"Embeddings desc {p}: shape={last_emb.shape}, dtype={last_emb.dtype}")
            if tqdm_bar is not None:
                tqdm_bar.update(1)
        return embeddings

    def write_all_embeddings(self, fasta_path: str, output_prefix: str, 
            poolings: List[str], allow_missing: bool = False, subset_seq: set = None):
        print("Reading original fasta file...")
        fasta_content = read_fasta(fasta_path)
        if subset_seq is not None:
            fasta_content = [(seq_id, seq) 
                for seq_id, seq in fasta_content if seq in subset_seq]
        ids = [seq_id for seq_id, _ in fasta_content]
        seqs = [seq for _, seq in fasta_content]

        # 1. Create a Polars DataFrame of the exact original FASTA order
        df_fasta = pl.DataFrame({"id": ids, "seq": seqs}).with_row_index("original_order")

        all_cache_files = self.cache.list_caches()
        valid_pqs = [pq_path for pq_path, txt_path in all_cache_files]

        #print("Lazy scanning base embeddings via Polars...")
        #df_embs_base = pl.scan_parquet(valid_pqs)

        # 2. Iterativamente constrói e salva UM arquivo parquet por pooling
        for pooling in poolings:
            target_parquet = f"{output_prefix}_{pooling}.parquet"
            print(f"\n--- Iniciando o merge do pooling: [{pooling}] ---")

            print("Lazy scanning base embeddings via Polars...")
            #df_embs_base = pl.scan_parquet(valid_pqs)
            paths_with_pooling = []
            for p in valid_pqs:
                col_list = pl.scan_parquet(p).collect_schema().names()
                if pooling in col_list:
                    paths_with_pooling.append(p)
            if len(paths_with_pooling) > 0:
                lfs = [
                    pl.scan_parquet(p).select(["seq", pooling]) 
                    for p in paths_with_pooling
                ]
                #df_embs_base = pl.scan_parquet(paths_with_pooling)
                print(f"Found {len(paths_with_pooling)} parquet files for pooling {pooling}.")
                print("Seleciona APENAS a sequência e este pooling específico")
                df_embs = pl.concat(lfs, how="vertical")
                df_embs = df_embs.unique(subset=["seq"], keep="first")

                print("Build the computation graph")
                df_final_lazy = (
                    df_fasta.lazy()
                    .join(df_embs, on="seq", how="left")
                    .sort("original_order")
                    .drop("original_order")
                )
                #print(df_final_lazy)

                try:
                    print(f"Streaming direto para o disco: {target_parquet}")
                    df_final_lazy.sink_parquet(target_parquet)
                except Exception as e:
                    print(f"Erro ao processar o pooling {pooling}: {e}")
                if not allow_missing:
                    print(f"Checando por sequências perdidas em {pooling}...")
                    missing_count = (
                        pl.scan_parquet(target_parquet)
                        .filter(pl.col(pooling).is_null())
                        .select(pl.len())
                        .collect()
                        .item()
                    )
                    assert missing_count == 0, f"Missing embeddings for {missing_count} sequences in {pooling}!"
    
    def embed_saving_progress(self, fasta_path: str, output_prefix: str, poolings: List[str] = list(POOLERS.keys()), max_tokens_per_batch = None):
        if max_tokens_per_batch is None:
            max_tokens_per_batch = self.max_tokens
        caching_batch_len = max_tokens_per_batch * 200
        not_embedded_seqs, embedded_seqs = self.cache.list_non_embedded(fasta_path, poolings=poolings)
        macro_batchs = batch_by_tokens(not_embedded_seqs, caching_batch_len, use_padding=False)
        small_batchs = [batch_by_tokens(seqs, max_tokens_per_batch, use_padding=True) 
            for seqs in macro_batchs]
        n_micro_batchs = [len(b2) for b2 in small_batchs]
        total_micro_batchs = sum(n_micro_batchs)
        bar = tqdm(total=total_micro_batchs, desc="Embedding sequences")

        no_full = [p for p in poolings if p != "full"]
        
        # Salva o arquivo iterativo inicial
        #try:
        self.write_all_embeddings(fasta_path, f"{output_prefix}_incomplete", 
            no_full, allow_missing=True, subset_seq = embedded_seqs)
        #except Exception as ex:
        #    print(ex)
        #    pass
            
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
            
            print(f"Saving batch to {next_parquet}...")
            df = pl.DataFrame(df_dict)
            df.write_parquet(next_parquet)
            with open(next_txt, "w") as f:
                for seq in macro_batch:
                    f.write(f"{seq}\n")
            print(f"Saved cache to {next_parquet} and {next_txt}")
            duration = time.time() - time_start
            print(f"Duration: {duration}")
            del df
            del df_dict
            del new_embeddings
            gc.collect()
            
            # Atualiza os parquets separados a cada macro-batch
            #self.write_all_embeddings(fasta_path, f"{output_prefix}_incomplete", poolings, allow_missing=True)

        # Escrita final garantindo que nada está faltando
        self.write_all_embeddings(fasta_path, output_prefix, no_full)