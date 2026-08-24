import os
from collections import Counter
import glob

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.manifold import TSNE


taxid_pca_labels = {
    # --- EUKARYOTA ---
    "40674": {"name": "Mammals", "color": "darkred"},
    "7742": {"name": "Other Vertebrates", "color": "hotpink"},
    "6656": {"name": "Arthropods & Insects", "color": "#FF9A00"},
    "6231": {"name": "Nematodes", "color": "#E9C46A"},  # C. elegans & worms
    "4751": {"name": "Fungi", "color": "magenta"},
    "33090": {"name": "Plants & Algae", "color": "forestgreen"},
    "5794": {"name": "Apicomplexans", "color": "#A8DADC"},  # Plasmodium / Parasites
    # --- BACTERIA ---
    "3379134": {"name": "Bacteria - Pseudomonadati", "color": "cornflowerblue"},
    "1783272": {"name": "Bacteria - Bacillati", "color": "dodgerblue"},
    "201174": {
        "name": "Bacteria - Actinomycetota",
        "color": "midnightblue",
    },  # Mycobacterium, etc.
    "2": {"name": "Bacteria - Others", "color": "#E0E1DD"},
    # --- VIRUSES ---
    "2732396": {"name": "Viruses - RNA", "color": "#457B9D"},  # Distinct cool blue
    "10239": {
        "name": "Viruses - DNA & Others",
        "color": "#1D3557",
    },  # Deep navy/graphite
    # --- NOISE ---
    "other": {"name": "Other Organisms", "color": "silver"},  # Very light gray
}


def prepare_pca_clusters_taxo(eval_data, top_n=16):
    """
    Finds the top N terms in the evaluation set and assigns proteins
    to strictly mutually exclusive clusters for clean visualization.
    """
    print(f"Preparing PCA clusters for top {top_n} terms...")

    # 1. Find the top N most common terms in the eval set
    all_eval_terms = [term for protein in eval_data for term in protein]

    target_indices = []
    cluster_labels = []

    other_label = "other"

    for i, protein in enumerate(eval_data):
        label = "other"
        for taxid, info in taxid_pca_labels.items():
            if taxid in protein:
                label = taxid
                break
        # pretty_label = taxid_pca_labels[label]["name"]
        cluster_labels.append(label)
        target_indices.append(i)

    undefined_indexes = [i for i in target_indices if cluster_labels[i] == other_label]

    print("Undefined organisms:", len(undefined_indexes) / len(eval_data) * 100, "%")

    print(eval_data[0])
    print(eval_data[-1])
    undefined_eval_data = [
        [x for x in eval_data[i] if x not in taxid_pca_labels]
        for i in undefined_indexes
    ]

    term_list_for_freq = []
    print(undefined_eval_data[0])
    print(undefined_eval_data[-1])
    for lineage in undefined_eval_data:
        for taxid in lineage:
            term_list_for_freq.append(taxid)

    term_counts = Counter(term_list_for_freq)
    top_terms = [term for term, _ in term_counts.most_common(32)]
    print("Top most common terms in undefined organisms:")
    for x in top_terms:
        print(x, term_counts[x])

    print(f"Selected {len(target_indices)} different lineages for PCA.")
    different_labels = len(set(cluster_labels))
    print(f"Different labels: {different_labels}")

    label_counts = Counter(cluster_labels)

    print("Frequencies of used labels:")
    for label, count in label_counts.items():
        print(taxid_pca_labels[label]["name"], count)

    labels_used = list(set(cluster_labels))

    return target_indices, cluster_labels, labels_used


def plot_epoch_pca(
    eval_embeddings, target_indices, cluster_labels, top_terms, epoch, directory
):
    """
    Calculates 2D t-SNE for the selected embeddings and saves a clean, ordered plot.
    """

    X = eval_embeddings[target_indices]

    tsne = TSNE(
        n_components=2, perplexity=30, learning_rate="auto", init="pca", random_state=42
    )
    X_pca = tsne.fit_transform(X)

    if len(top_terms) <= 20:
        fig, ax = plt.subplots(
            figsize=(4.8, 4.5), dpi=180
        )  # Slightly wider for the legend
    else:
        fig, ax = plt.subplots(figsize=(8, 7), dpi=180)

    is_taxid = all([t in taxid_pca_labels for t in cluster_labels])
    if is_taxid:
        colors = {term: taxid_pca_labels[term]["color"] for term in top_terms}
    else:
        if len(top_terms) <= 20:
            cmap = plt.get_cmap("tab20")
            colors = {term: cmap(i) for i, term in enumerate(top_terms)}
        else:
            n_colors_to_pick = len(top_terms)
            cmap = plt.get_cmap("hsv")

            # sample evenly from the cmap
            n_colors_in_colormap = cmap.N
            colors = {
                term: cmap(int(i * n_colors_in_colormap / n_colors_to_pick))
                for i, term in enumerate(top_terms)
            }

    # 1. Scatter plot by cluster
    for term in top_terms:
        term_mask = [label == term for label in cluster_labels]

        if any(term_mask):
            # We skip adding the 'label' argument here to prevent default legend creation
            ax.scatter(
                X_pca[term_mask, 0],
                X_pca[term_mask, 1],
                color=colors[term],
                alpha=0.95,
                s=12,  # Keep the actual plot points small and distinct
                edgecolors="none",
            )

    # 2. Build Custom Legend
    legend_handles = []

    if is_taxid:
        # Iterate through the dictionary to maintain exact taxonomic order
        for taxid, info in taxid_pca_labels.items():
            # Only add to the legend if the organism actually exists in this evaluation batch
            if taxid in top_terms or (taxid == "other" and "other" in cluster_labels):
                handle = mlines.Line2D(
                    [],
                    [],
                    color=info["color"],
                    marker="o",
                    linestyle="None",
                    markersize=9,  # LARGER CIRCLES in the legend
                    label=info["name"],
                )
                legend_handles.append(handle)
        legend_fontsize = 7
        labelspacing = 0.8
    else:
        for term in top_terms:
            handle = mlines.Line2D(
                [],
                [],
                color=colors[term],
                marker="o",
                linestyle="None",
                markersize=6,
                label=term,
            )
            legend_handles.append(handle)
        legend_fontsize = 6
        labelspacing = 0.4

    # 3. Formatting
    title = (
        f"Compact Representation of Taxonomy - {epoch}"
        if is_taxid
        else f"Compact Representation of Protein Families - {epoch}"
    )
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xticks([])
    ax.set_yticks([])

    # Apply the custom legend
    ax.legend(
        handles=legend_handles,
        fontsize=legend_fontsize,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),  # Push legend entirely outside the plot box
        frameon=False,
        title="Taxonomic Clades" if is_taxid else "Prot. Family IDs",
        title_fontsize=8,
        labelspacing=labelspacing,  # Add a little breathing room between items
    )

    # 4. Save the plot
    if directory is None:
        directory = "."
    os.makedirs(directory, exist_ok=True)

    filepath = os.path.join(directory, f"pca_epoch_{epoch:03d}.png")
    plt.savefig(filepath, format="png", bbox_inches="tight")
    plt.close(fig)


def plot_tsne_scatter_multilabel_pies(
    eval_embeddings,
    target_indices,
    cluster_labels,
    top_terms,
    epoch,
    directory,
    pretty_names={},
    max_other_points=800,
):
    focal_indices = []
    focal_labels = []
    other_indices = []
    other_labels = []

    # 1. Separate "Other" from the biologically distinct clusters
    for idx, label in zip(target_indices, cluster_labels):
        # Catch both capitalizations just in case (Interpro vs Taxid)
        if label == "Other" or label == "other":
            other_indices.append(idx)
            other_labels.append(label)
        else:
            focal_indices.append(idx)
            focal_labels.append(label)

    # 2. Subsample ONLY the "Other" category
    if len(other_indices) > max_other_points:
        np.random.seed(
            42
        )  # Fixed seed so the background scaffolding doesn't jitter wildly between epochs
        sampled_idxs = np.random.choice(
            len(other_indices), max_other_points, replace=False
        )
        other_indices = [other_indices[i] for i in sampled_idxs]
        other_labels = [other_labels[i] for i in sampled_idxs]

    # 3. Recombine for t-SNE
    target_indices = focal_indices + other_indices
    cluster_labels = focal_labels + other_labels

    X = eval_embeddings[target_indices]

    tsne = TSNE(
        n_components=2,
        perplexity=30,
        learning_rate="auto",
        init="pca",
        random_state=42,
        n_jobs=4,
    )
    X_pca = tsne.fit_transform(X)

    fig, ax = plt.subplots(figsize=(9, 6), dpi=180)

    n_colors_to_pick = len(top_terms)
    cmap = plt.get_cmap("nipy_spectral")
    # sample evenly from the cmap
    n_colors_in_colormap = cmap.N
    colors = {
        term: cmap(int(i * n_colors_in_colormap / n_colors_to_pick))
        for i, term in enumerate(top_terms)
    }

    # 1. Scatter plot by cluster
    for term in top_terms:
        term_mask = [label == term for label in cluster_labels]

        if any(term_mask):
            # We skip adding the 'label' argument here to prevent default legend creation
            ax.scatter(
                X_pca[term_mask, 0],
                X_pca[term_mask, 1],
                color=colors[term],
                alpha=0.95,
                s=12,  # Keep the actual plot points small and distinct
                edgecolors="none",
            )

    # 2. Build Custom Legend
    legend_handles = []
    for term in top_terms:
        term2 = pretty_names.get(term, term)
        handle = mlines.Line2D(
            [],
            [],
            color=colors[term],
            marker="o",
            linestyle="None",
            markersize=6,
            label=term2,
        )
        legend_handles.append(handle)
    legend_fontsize = 6
    labelspacing = 0.4

    # 3. Formatting
    title = f"Compact Representation of Protein Families - {epoch}"
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xticks([])
    ax.set_yticks([])

    # Apply the custom legend
    ax.legend(
        handles=legend_handles,
        fontsize=legend_fontsize,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),  # Push legend entirely outside the plot box
        frameon=False,
        title="Prot. Families",
        title_fontsize=8,
        labelspacing=labelspacing,  # Add a little breathing room between items
    )

    # 4. Save the plot
    if directory is None:
        directory = "."
    os.makedirs(directory, exist_ok=True)

    filepath = os.path.join(directory, f"pca_epoch_{epoch:03d}.png")
    plt.savefig(filepath, format="png", bbox_inches="tight")
    plt.close(fig)


def generate_pca_gif(directory, output_filename="latent_evolution.gif", duration=400):
    """
    Compiles all saved PCA scatter plots into an animated GIF.
    duration: milliseconds per frame (400ms = 2.5 frames per second)
    """
    if directory is None:
        return

    print(f"\nStitching PCA frames into GIF...")

    # Find all pca images and sort them alphabetically (which matches chronologically due to the 000 padding)
    search_pattern = os.path.join(directory, "pca_epoch_*.png")
    image_paths = sorted(glob.glob(search_pattern))

    if not image_paths:
        print("No PCA images found to make a GIF.")
        return

    # Open all images
    frames = [Image.open(img_path) for img_path in image_paths]

    # Save as GIF
    output_path = os.path.join(directory, output_filename)

    # Save the first frame, and append the rest
    frames[0].save(
        output_path,
        format="GIF",
        append_images=frames[1:],
        save_all=True,
        duration=duration,
        loop=0,  # 0 means infinite loop
    )
    print(f"GIF successfully updated at: {output_path}")
