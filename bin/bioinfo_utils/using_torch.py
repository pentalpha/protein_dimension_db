import torch


def embedding_neighborhood_score(
    eval_targets, eval_embeddings, k=5, min_overlap_ratio=0.6
):
    """
    Computes the average fraction of k-nearest neighbours that share
    at least `min_overlap_ratio` of terms with the query protein.

    Args:
        eval_targets: torch.Tensor of shape (N, num_terms), binary.
        eval_embeddings: torch.Tensor of shape (N, dim).
        k: number of neighbours to consider.
        min_overlap_ratio: float (0.0 to 1.0). The minimum proportion of the query's
                           terms that the neighbour must also possess to be a "match".

    Returns:
        float: average neighbour label agreement.
    """
    # 1. Normalise embeddings for cosine similarity
    normed_emb = eval_embeddings / (eval_embeddings.norm(dim=1, keepdim=True) + 1e-8)

    # 2. Compute pairwise cosine similarity matrix
    sim = torch.mm(normed_emb, normed_emb.t())  # (N, N)

    # Set diagonal to -inf so the query itself isn't chosen
    sim.fill_diagonal_(-float("inf"))

    # 3. Get top‑k indices for each row
    _, knn_indices = torch.topk(sim, k, dim=1)  # (N, k)

    # 4. Fetch the target vectors for all k neighbours of all N queries
    # Shape: (N, k, num_terms)
    knn_targets = eval_targets[knn_indices]

    # Prepare query targets for broadcasting
    # Shape: (N, 1, num_terms)
    query_targets = eval_targets.unsqueeze(1)

    # 5. Calculate intersection size (number of shared terms)
    # Shape: (N, k)
    intersection_sizes = (query_targets * knn_targets).sum(dim=2)

    # 6. Calculate the number of terms the query actually has
    # Shape: (N, 1)
    query_sizes = query_targets.sum(dim=2)

    # Prevent division by zero if a query has absolutely no labels
    # (Though in taxonomic lineages, this shouldn't happen)
    query_sizes = torch.clamp(query_sizes, min=1.0)

    # 7. Calculate the overlap ratio (Intersection / Query Size)
    # Shape: (N, k)
    overlap_ratios = intersection_sizes / query_sizes

    # 8. Check which neighbours meet the minimum threshold
    # Shape: (N, k) -> boolean mask converted to float (1.0 or 0.0)
    matches = (overlap_ratios >= min_overlap_ratio).float()

    # 9. Calculate the fraction of the k neighbours that are matches for each query
    # Shape: (N,)
    agreement_fracs = matches.mean(dim=1)

    # Return the global mean score
    return agreement_fracs.mean().item()
