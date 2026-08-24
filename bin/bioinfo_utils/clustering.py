from collections import Counter

TERM_FREQUENCY_RANGES = [(30, 90), (91, 400), (401, 3000)]


def prepare_pca_clusters(eval_data, top_n=14):
    """
    Finds the top N terms in the evaluation set and assigns proteins
    to strictly mutually exclusive clusters for clean visualization.
    """
    print(f"Preparing PCA clusters for top {top_n} terms...")

    # 1. Find the top N most common terms in the eval set
    all_eval_terms = [term for protein in eval_data for term in protein]
    term_counts = Counter(all_eval_terms)
    top_terms = [term for term, _ in term_counts.most_common(top_n)]
    top_terms_set = set(top_terms)

    target_indices = []
    cluster_labels = []

    # 2. Filter proteins: must contain EXACTLY ONE of the top N terms
    for i, protein in enumerate(eval_data):
        if len(protein) >= 3:
            overlap = top_terms_set.intersection(protein)
            if len(overlap) <= 2 and len(overlap) > 0:
                usable_labels = list(overlap)
                usable_labels.sort(key=lambda l: term_counts[l])
                label = usable_labels[0]
                target_indices.append(i)
                cluster_labels.append(label)

    print(f"Selected {len(target_indices)} purely clustered proteins for PCA.")
    different_labels = len(set(cluster_labels))
    print(f"Different labels: {different_labels}")
    return target_indices, cluster_labels, top_terms


def prepare_pca_clusters_predef(eval_data, top_terms, include_unclustered=False):
    target_indices = []
    cluster_labels = []
    cluster_multilabels = []

    for i, protein in enumerate(eval_data):
        prot_label = "Other"
        prot_mlabel = "Other"
        clusters_in_labels = [label for label in top_terms if label in protein]
        if len(clusters_in_labels) > 0:
            prot_label = clusters_in_labels[0]
            if len(clusters_in_labels) > 1:
                prot_mlabel = ";".join(clusters_in_labels)
            else:
                prot_mlabel = prot_label

        if prot_label != "Other" or include_unclustered:
            target_indices.append(i)
            cluster_labels.append(prot_label)
            cluster_multilabels.append(prot_mlabel)

    print(f"Selected {len(target_indices)} clustered proteins for PCA.")
    different_labels = len(set(cluster_labels))
    print(f"Different labels: {different_labels}")
    return target_indices, cluster_labels, cluster_multilabels, top_terms


def create_mutually_exclusive_clustering(label_options, term_freqs, clean_data):
    # Creates a labeling of clean_data, where each annotation has a single representative label (or none)
    # Members of a cluster share the same label and members of other clusters dont have that label

    mutually_exclusive_labels = []  # clusters
    label_options.sort(key=lambda x: term_freqs[x])

    # inspect each label to try to make clusters from each, starting from least frequent
    new_cluster_added = True
    while new_cluster_added and len(label_options) > 0:
        label = label_options[0]
        indexes_with_label = [i for i, x in enumerate(clean_data) if label in x]
        intersecting_labels_set = set()
        for i in indexes_with_label:
            for term in clean_data[i]:
                if term in label_options and term != label:
                    intersecting_labels_set.add(term)
        mutually_exclusive_labels.append(label)
        to_remove = [label] + list(intersecting_labels_set)
        for term in to_remove:
            label_options.remove(term)
        new_cluster_added = True
    return mutually_exclusive_labels


def find_intepro_clusterings(term_freqs, clean_data):

    cluster_lists = {}
    for min_freq, max_freq in TERM_FREQUENCY_RANGES:
        range_name = f"freq_{min_freq}_to_{max_freq}"
        print(f"Creating clusters for range: {range_name}")

        label_options = [
            x for x, y in term_freqs.items() if y >= min_freq and y <= max_freq
        ]
        label_options_copy = label_options.copy()
        clusters = create_mutually_exclusive_clustering(
            label_options_copy, term_freqs, clean_data
        )

        print(f"Found {len(clusters)} clusters for range {range_name}")
        cluster_lists[range_name] = clusters

    return cluster_lists
