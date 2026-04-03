#!/usr/bin/env python
import sys
import gzip
from collections import defaultdict
from bioinfo_utils.gene_ontology import expand_go_set, load_go_graph, gos_not_to_use


def format_cafa_terms(input_path, output_prefix, go_not_use_path, go_basic_path):
    # EntryID	term	aspect
    # Q5W0B1	GO:0000785	C

    # Store annotations: aspect -> entry_id -> list of terms
    annotations = {
        "F": defaultdict(set),  # MF
        "P": defaultdict(set),  # BP
        "C": defaultdict(set),  # CC
    }

    go_graph = load_go_graph(go_basic_path)
    not_use = gos_not_to_use(go_not_use_path)

    print(f"Reading {input_path}")
    opener = gzip.open if input_path.endswith(".gz") else open
    mode = "rt" if input_path.endswith(".gz") else "r"

    with opener(input_path, mode) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            entry_id, term, aspect = parts
            if aspect in annotations:
                annotations[aspect][entry_id].add(term)

    # Mapping aspect char to output aspect string and filename suffix
    aspect_map = {"F": "mf", "P": "bp", "C": "cc"}

    for aspect_char, aspect_name in aspect_map.items():
        output_file = f"{output_prefix}.{aspect_name}.tsv"
        print(f"Writing {output_file}")
        with open(output_file, "w") as f:
            for entry_id, terms in annotations[aspect_char].items():
                terms_updated = set(terms)
                for term in terms:
                    expanded_terms = expand_go_set(term, go_graph, not_use)
                    terms_updated.update(expanded_terms)
                terms_str = ",".join(sorted(list(terms_updated)))
                f.write(f"{entry_id}\t{terms_str}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python format_cafa_terms.py <input_train_terms> <output_prefix>")
        sys.exit(1)

    input_terms = sys.argv[1]
    output_prefix = sys.argv[2]  # e.g. "go.experimental"
    go_not_use_path = sys.argv[3]
    go_basic_path = sys.argv[4]

    format_cafa_terms(input_terms, output_prefix, go_not_use_path, go_basic_path)
