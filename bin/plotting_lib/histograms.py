from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import LogLocator, LogFormatterMathtext


def plot_term_freq_histogram(values_for_freq: list, model_dir: str):
    n_bins = 32
    max_freq = 1000

    values_freq = Counter(values_for_freq)
    x1 = np.array(list(values_freq.values()))

    # Split the data into regular range and overflow
    ones = x1[x1 == 1]
    regular = x1[(x1 > 1) & (x1 <= max_freq)]
    overflow = x1[x1 > max_freq]
    global_max = max(x1) if len(x1) > 0 else 0

    if len(regular) == 0 and len(ones) == 0:
        raise ValueError(
            "No values <= 1000 were found. Cannot build the 16-bin histogram."
        )

    if len(regular) == 0:
        raise ValueError(
            "No values > 1 and <= 1000 were found. Cannot build the 16-bin histogram."
        )

    # 16 bins between 2 and 1000
    bin_edges_reg = np.linspace(2, max_freq, n_bins + 1)
    counts_reg, _ = np.histogram(regular, bins=bin_edges_reg)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(16, 6), constrained_layout=True)

    # 1. Plot the "ones" (Singletons)
    # Positioned at x=0
    bars_ones = ax.bar(0, [len(ones)], width=0.8, edgecolor="black", linewidth=0.8)

    # 2. Plot the "regular" bins
    # Positioned from x=2 to x=17. Width 1.0 makes the bars touch.
    x_reg = np.arange(2, 2 + n_bins)
    bars_reg = ax.bar(x_reg, counts_reg, width=1.0, edgecolor="black", linewidth=0.8)

    # 3. Plot the "overflow" bin
    # Positioned with a symmetric gap after the last regular bar
    x_over = x_reg[-1] + 2
    bar_over = ax.bar(
        x_over, [len(overflow)], width=0.8, edgecolor="black", linewidth=0.8
    )

    # Make the overflow bar visually distinct
    bar_over[0].set_hatch("//")
    bar_over[0].set_alpha(0.85)

    ax.set_title("Distribution of InterPro Term Frequencies in the Dataset", pad=12)
    ax.set_xlabel("Term frequency bin")
    ax.set_ylabel("Number of terms")

    # --- Tick and Label Formatting ---
    # Ticks for regular bins are placed at the edges of the bars (x +/- 0.5)
    reg_tick_positions = np.arange(1.5, 1.5 + n_bins + 1)
    tick_positions = [0] + list(reg_tick_positions) + [x_over]

    # Generate labels: Discrete "1", the actual edges for regular bins, and the overflow range
    tick_labels = (
        ["1"] + [str(int(e)) for e in bin_edges_reg] + [f"{max_freq}–{global_max}"]
    )

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")

    # Scientific-looking y-axis
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(LogFormatterMathtext())
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))

    # Clean styling
    ax.grid(axis="y", which="major", linestyle="--", linewidth=0.7, alpha=0.5)
    ax.grid(axis="y", which="minor", linestyle=":", linewidth=0.5, alpha=0.25)
    ax.grid(
        axis="x", visible=False
    )  # Turned off x-grid to prevent clutter around edges
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate bars
    all_rects = [bars_ones[0]] + list(bars_reg) + [bar_over[0]]
    all_counts = [len(ones)] + counts_reg.tolist() + [len(overflow)]

    for rect, c in zip(all_rects, all_counts):
        if c > 0:
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                c * 1.08,
                f"{c}",
                ha="center",
                va="bottom",
                fontsize=9,
                rotation=0,
            )

    os.makedirs(model_dir, exist_ok=True)
    fig.savefig(
        os.path.join(model_dir, "term_freq_histogram.png"), dpi=400, bbox_inches="tight"
    )
    plt.close(fig)
