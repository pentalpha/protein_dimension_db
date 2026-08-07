import numpy as np
from typing import List, Union

def mean_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    # Iterate over list, mean along axis 0 (sequence length), stack into (Batch, 768)
    return np.vstack([np.mean(emb, axis=0) for emb in emb_list])

def max_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return np.vstack([np.max(emb, axis=0) for emb in emb_list])

def std_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return np.vstack([np.std(emb, axis=0) for emb in emb_list])

def norm_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return np.vstack([np.linalg.norm(emb, ord=2, axis=0) for emb in emb_list])

def softmax_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    pooled_results = []
    for emb in emb_list:
        # Embs are (L, 768). Find max across L for stability
        max_emb = np.max(emb, axis=0, keepdims=True) 
        exp_emb = np.exp(emb - max_emb)
        pooled = np.squeeze(max_emb, axis=0) + np.log(np.mean(exp_emb, axis=0))
        pooled_results.append(pooled)
    return np.vstack(pooled_results)

# --- Updated K-Max Pooling Functions ---

def _k_max_pooling_core(emb_list: List[np.ndarray], percentage: float) -> np.ndarray:
    pooled_results = []
    # DYNAMICALLY get the embedding dimension (axis 1) from the first sequence
    hidden_dim = emb_list[0].shape[1]
    
    # Calculate K based on the actual model's dimension
    k = int(round(hidden_dim * percentage))
    k = max(1, k) # Safety catch: ensure we pool at least 1 token
    
    for emb in emb_list:
        seq_len = emb.shape[0]  # seq_len is now axis 0
        actual_k = min(k, seq_len)
        
        if actual_k < seq_len:
            kth = seq_len - actual_k
            # Partition along axis 0 (sequence length)
            partitioned = np.partition(emb, kth, axis=0)
            # Take the top K elements (from kth index to the end)
            top_k_values = partitioned[kth:, :]
        else:
            top_k_values = emb
            
        pooled_results.append(np.mean(top_k_values, axis=0))
        
    return np.vstack(pooled_results)

def k4p_max_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return _k_max_pooling_core(emb_list, 0.04)

def k8p_max_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return _k_max_pooling_core(emb_list, 0.08)

def k16p_max_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return _k_max_pooling_core(emb_list, 0.16)

def k32p_max_pooling(emb_list: List[np.ndarray]) -> np.ndarray:
    return _k_max_pooling_core(emb_list, 0.32)

# --- Updated Dictionary ---

POOLERS = {
    "full": lambda x: x,
    "mean": mean_pooling,
    "max": max_pooling,
    "softmax": softmax_pooling,
    "std": std_pooling,
    "norm": norm_pooling,
    "k4p_max": k4p_max_pooling,
    "k8p_max": k8p_max_pooling,
    "k16p_max": k16p_max_pooling,
    "k32p_max": k32p_max_pooling,
}
