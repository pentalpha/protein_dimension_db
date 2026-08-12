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

def _pagerank(
    A: np.ndarray,
    damping: float = 0.85,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> np.ndarray:
    """Computes PageRank centrality for a given directed attention graph."""
    # Ensure FP64 for iterative numerical stability
    A = np.asarray(A, dtype=np.float64)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"Attention matrix must be square, got {A.shape}")

    L = A.shape[0]

    if L == 0:
        return np.array([])
    if L == 1:
        return np.ones(1, dtype=np.float64)

    # Defensive catch against negative attention values
    A = np.maximum(A, 0.0)

    # Calculate row sums to normalize into transition probabilities
    row_sums = A.sum(axis=1)
    
    # Identify "dangling nodes" (residues that pay zero attention to anything else)
    dangling = row_sums == 0

    # Row-normalize, using a tiny float to prevent division by zero
    A = A / np.maximum(row_sums[:, None], np.finfo(np.float64).tiny)

    # Fix dangling nodes: they distribute probability uniformly across the sequence
    if np.any(dangling):
        A[dangling, :] = 1.0 / L

    # Initialize uniform probability mass
    p = np.full(L, 1.0 / L)
    teleport = (1.0 - damping) / L

    # Power Iteration
    for _ in range(max_iter):
        # Transpose A to match the stochastic transition logic (columns = incoming edges)
        p_next = teleport + damping * (A.T @ p)

        # Check L1 norm for convergence
        if np.linalg.norm(p_next - p, ord=1) < tol:
            p = p_next
            break
        p = p_next

    # Final exact normalization to ensure probabilities sum perfectly to 1
    p /= p.sum()

    return p

def parti_pooling(
    emb_list: List[np.ndarray],
    attn_list: List[np.ndarray], # Now required by your pipeline
    damping: float = 0.85,
) -> np.ndarray:
    """
    Applies PoolPaRTI over a batch of extracted sequence embeddings.
    """
    if len(emb_list) != len(attn_list):
        raise ValueError("Embeddings and attention lists must have the same length.")

    results = []

    for X, A in zip(emb_list, attn_list):
        if X.shape[0] == 0:
            continue
            
        # Get PageRank importances
        weights = _pagerank(A, damping=damping)

        # Weighted sum: weights @ X creates a (D,) vector by broadcasting
        pooled = weights @ X
        results.append(pooled)

    return np.vstack(results)

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
    "parti": parti_pooling,
}
