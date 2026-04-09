import math
import heapq
import json

"""
## Creating Informative and Broad Compact Vocabularies

To construct a representative and information-dense vocabulary for the Protein Dimension DB autoencoders, we formulated the selection process as a greedy maximization problem. This approach iteratively identifies the most valuable terms to include in the vocabulary by evaluating a custom utility function designed to optimize both global coverage and specific functional depth. 

Global coverage was defined as the fraction of proteins that are annotated with at least one term in the vocabulary. Specificity as defined as Information Content (IC). Following the probabilistic framework established by Resnik (1995) and further applied in Gene Ontology (GO) assessments by Wang et al. (2014), we define IC from empirical annotation frequency as:

IC(term) = -log p(term).

To ensure computational feasibility across large datasets, this optimization is solved using a CELF-style (Cost-Effective Lazy Forward) selection procedure (Leskovec 2007). This procedure favors terms that increased protein coverage while still allowing highly specific terms to be retained for already covered proteins. This yields a reduced vocabulary that preserves broad functional and taxonomic diversity while maintaining high informational density. The overall objective is related to prior work on informative subset selection and ontology graph compression, although our implementation differs from classical information bottleneck formulations by directly optimizing a coverage-weighted IC utility.

## CELF-based Vocabulary Selection Algorithm

The implementation utilizes a lazy greedy approach to maximize the Average Total IC Retained while managing computational cost:

### 1. Precomputation:
- Compute IC(term) for all available classes.
- Calculate the total original IC sum for each protein (IC_total(p)) to allow O(1) tracking of the incremental fraction of information retained.
### 2. Heap Construction:
- Build a max-priority queue Q containing all possible terms.
- Initial scores are set as IC(term) * |Proteins(term)|, representing the maximum potential information gain.
### 3. Lazy Greedy Iteration:
- Pop the top candidate term from Q.
- Re-evaluation: Recalculate the Score(term) based on the current covered protein set, accounting for new coverage and the redundancy bonus alpha.
- Check Competition: Compare the updated score against the current top of Q.
- Selection: If Score(term) >= top(Q), add term to the vocabulary and update the global covered set.
- Lazy Step: Otherwise, re-insert term into Q with its updated score and repeat.
## 4. Termination:
- The process continues until the Average Total IC Retained target is met or the maximum vocabulary size is reached.

Resnik, P. (1995). Using Information Content to Evaluate Semantic Similarity in a Taxonomy. Proceedings of the 14th IJCAI.

Leskovec, J., Krause, A., Guestrin, C., Faloutsos, C., VanBriesen, J., & Glance, N. (2007). Cost-effective Outbreak Detection in Networks. KDD 2007.

Wang, H., Azuaje, F., & Zheng, H. (2014). An information theoretic approach to assessing Gene-Ontology-driven similarity and its application. Int. J. Data Mining and Bioinformatics, 9(2), 121–134.
"""

# Final class:


def calc_coverage(vocab, annots_by_class, all_uniprots):
    vocab_uniprots = set()
    for c in vocab:
        vocab_uniprots.update(annots_by_class[c])
    return len(vocab_uniprots) / len(all_uniprots)


def build_ic_map(annots_by_class, all_uniprots):
    N = len(all_uniprots)
    ic_map = {}
    for c, prots in annots_by_class.items():
        n = len(prots)
        if n == 0:
            continue
        p_c = n / N
        ic_map[c] = -math.log2(p_c)
    return ic_map


def evaluate_annotation_quality(all_uniprots, annots_by_class, vocab, ic_map):
    """
    Evaluates the quality of a reduced vocabulary by measuring
    how much Information Content (IC) is retained per protein.
    """
    vocab_set = set(vocab)

    # 1. Invert the mapping: Protein -> Original Classes
    protein_to_classes = {p: [] for p in all_uniprots}
    for c, prots in annots_by_class.items():
        for p in prots:
            protein_to_classes[p].append(c)

    total_ic_ratios = []
    max_ic_ratios = []

    for p, orig_classes in protein_to_classes.items():
        if not orig_classes:
            continue

        # Original IC metrics
        orig_ic_sum = sum(ic_map[c] for c in orig_classes)
        orig_ic_max = max(ic_map[c] for c in orig_classes)

        # Vocabulary-filtered metrics
        vocab_classes = [c for c in orig_classes if c in vocab_set]

        if not vocab_classes:
            total_ic_ratios.append(0.0)
            max_ic_ratios.append(0.0)
            continue

        vocab_ic_sum = sum(ic_map[c] for c in vocab_classes)
        vocab_ic_max = max(ic_map[c] for c in vocab_classes)

        total_ic_ratios.append(vocab_ic_sum / orig_ic_sum)
        max_ic_ratios.append(vocab_ic_max / orig_ic_max)

    # Calculate dataset-wide averages
    avg_total_ic = sum(total_ic_ratios) / len(total_ic_ratios)
    avg_max_ic = sum(max_ic_ratios) / len(max_ic_ratios)

    return avg_total_ic, avg_max_ic


def select_vocab_ic_weighted(
    annots_by_class,
    all_uniprots,
    target_ic_retained=0.90,  # Target metric: e.g., keep 90% of all Information Content
    alpha=0.1,  # How much to value adding IC depth to already-covered proteins
    max_vocab_size=24000,
    verbose=True,
    ic_map=None,
    total_ic_method="sum",
):
    """
    Greedy selection targeting 'Average Total IC Retained'.

    Score(c) = IA(c) * [new_proteins + (alpha * already_covered_proteins)]
    """

    print("IA Weighted Parameters:")
    print("target_ic_retained:", target_ic_retained)
    print("alpha:", alpha)
    print("max_vocab_size:", max_vocab_size)
    print("verbose:", verbose)

    N_total = len(all_uniprots)
    if ic_map is None:
        ic_map = build_ic_map(annots_by_class, all_uniprots)

    # 1. Precompute original IA sums per protein for O(1) tracking
    protein_orig_ic_sum = {}
    for c, prots in annots_by_class.items():
        for p in prots:
            if total_ic_method == "max":
                protein_orig_ic_sum[p] = max(protein_orig_ic_sum.get(p, 0.0), ic_map[c])
            else:
                protein_orig_ic_sum[p] = protein_orig_ic_sum.get(p, 0.0) + ic_map[c]

    # Count how many proteins actually have annotations
    N_annotated = sum(1 for v in protein_orig_ic_sum.values() if v > 0)

    # 2. Setup max-heap
    heap = []
    for c, prots in annots_by_class.items():
        if not prots:
            continue
        # Initially, all proteins are "new"
        initial_score = ic_map[c] * len(prots)
        heapq.heappush(heap, (-initial_score, c))

    covered = set()
    selected = []
    selected_set = set()

    # Incremental IC tracker
    current_retained_ratio_sum = 0.0
    avg_total_ic_retained = 0.0

    while heap and avg_total_ic_retained < target_ic_retained:
        if max_vocab_size is not None and len(selected) >= max_vocab_size:
            break

        neg_upper_score, c = heapq.heappop(heap)
        if c in selected_set:
            continue

        # Calculate exact score under current covered set
        total_prots_for_c = len(annots_by_class[c])
        remaining = annots_by_class[c] - covered
        new_count = len(remaining)
        already_covered_count = total_prots_for_c - new_count

        # Score now rewards BOTH new coverage AND adding depth to covered proteins
        exact_score = ic_map[c] * (new_count + (alpha * already_covered_count))

        # If it adds absolutely no value, skip it
        if exact_score <= 0:
            continue

        # Lazy greedy step
        if heap and exact_score < -heap[0][0]:
            heapq.heappush(heap, (-exact_score, c))
            continue

        # Accept this class
        selected.append(c)
        selected_set.add(c)
        covered.update(annots_by_class[c])

        """# Update incremental Average Total IC
        for p in annots_by_class[c]:
            if p in protein_orig_ic_sum and protein_orig_ic_sum[p] > 0:
                # The exact fractional increase this class provides to this protein
                delta = ic_map[c] / protein_orig_ic_sum[p]
                current_retained_ratio_sum += delta"""
        current_retained_ics = {}
        for c2 in selected:
            for p in annots_by_class[c2]:
                if p in covered:
                    if total_ic_method == "max":
                        current_retained_ics[p] = max(
                            current_retained_ics.get(p, 0.0), ic_map[c2]
                        )
                    else:
                        current_retained_ics[p] = (
                            current_retained_ics.get(p, 0.0) + ic_map[c2]
                        )
        current_retained_ratio_sum = sum(current_retained_ics.values())
        avg_total_ic_retained = current_retained_ratio_sum / N_annotated

        if verbose and (
            len(selected) <= 10
            or len(selected) % 500 == 0
            or avg_total_ic_retained >= target_ic_retained
        ):
            print(
                f"selected={len(selected):6d}  "
                f"ic_retained={avg_total_ic_retained:.4f}  "
                f"coverage={len(covered)/N_total:.4f}  "
                f"last={c}  "
                f"new_prots={new_count:5d}  "
                f"score={exact_score:.2f}"
            )

    return selected, covered, ic_map, avg_total_ic_retained


class ICRichVocabulary:
    def __init__(
        self,
        annots_by_class=None,
        input_path=None,
        target_ic_retained=0.95,
        alpha=0.2,
        max_vocab_size=24000,
        enrich_ic=True,
        vocab_method="ic_rich",
        ic_map=None,
        total_ic_method="sum",
    ):
        if annots_by_class is not None:
            self.target_ic_retained = target_ic_retained
            self.alpha = alpha
            self.max_vocab_size = max_vocab_size
            self.annots_by_class = annots_by_class

            self.instance_ids = set()
            self.classes = list(annots_by_class.keys())
            for c in self.classes:
                self.instance_ids.update(annots_by_class[c])

            if enrich_ic:
                self.vocab, self.covered, self.ic_map, self.retained_ic = (
                    select_vocab_ic_weighted(
                        annots_by_class=self.annots_by_class,
                        all_uniprots=self.instance_ids,
                        target_ic_retained=target_ic_retained,
                        alpha=alpha,
                        max_vocab_size=max_vocab_size,
                        verbose=True,
                        ic_map=ic_map,
                        total_ic_method=total_ic_method,
                    )
                )

            else:
                # Make vocab with max_k
                if len(self.classes) > max_vocab_size:
                    self.vocab = sorted(
                        self.classes,
                        key=lambda x: len(self.annots_by_class[x]),
                        reverse=True,
                    )[:max_vocab_size]
                else:
                    self.vocab = self.classes

                self.covered = set()
                for c in self.vocab:
                    self.covered.update(self.annots_by_class[c])
                if ic_map is None:
                    self.ic_map = build_ic_map(self.annots_by_class, self.instance_ids)
                else:
                    self.ic_map = ic_map
                self.retained_ic, _ = evaluate_annotation_quality(
                    self.instance_ids,
                    self.annots_by_class,
                    self.vocab,
                    self.ic_map,
                )
        elif input_path is not None:
            with open(input_path, "r") as f:
                data = json.load(f)
            self.target_ic_retained = data["target_ic_retained"]
            self.alpha = data["alpha"]
            self.max_vocab_size = data["max_vocab_size"]
            self.annots_by_class = data["annots_by_class"]
            self.instance_ids = set(data["instance_ids"])
            self.classes = data["classes"]
            self.vocab = data["vocab"]
            self.covered = set(data["covered"])
            self.ic_map = data["ic_map"]
            self.retained_ic = data["retained_ic"]
        else:
            raise ValueError("Either annots_by_class or input_path must be provided")

    def get_vocab(self, max_n=None):
        if max_n is None:
            return self.vocab
        else:
            vocab_sorted = sorted(
                self.vocab, key=lambda x: self.ic_map[x], reverse=True
            )
            return vocab_sorted[:max_n]

    def save(self, output_path):
        # Save class as json, but converting sets to lists before
        data = {
            "target_ic_retained": float(self.target_ic_retained),
            "alpha": float(self.alpha),
            "max_vocab_size": int(self.max_vocab_size),
            "vocab": self.vocab,
            "covered": list(self.covered),
            "ic_map": {k: float(v) for k, v in self.ic_map.items()},
            "retained_ic": float(self.retained_ic),
        }
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)
