#!/usr/bin/env python3

import sys

import numpy as np
from tqdm import tqdm

from bioinfo_utils.vocab_builder import ICRichVocabulary

if __name__ == "__main__":
    """Usage:
    make_vocab.py <input_type> <output_path>  \
        <min_vocab_size> <max_vocab_size> <target_ic_retained> \
        <vocab_method> <input_file1> <input_file2...>

    Example:
    make_vocab.py interpro ./data/interpro_vocab \
        64 21000 0.90 ic_rich \
        ./data/interproscan.tsv
    """
    input_type = sys.argv[1]
    assert input_type in ["interpro", "taxallnomy"]
    output_prefix = sys.argv[2]
    min_vocab_size = int(sys.argv[3])
    max_vocab_size = int(sys.argv[4])
    target_ic_retained = float(sys.argv[5])
    vocab_method = sys.argv[6]
    assert vocab_method in ["ic_rich", "top_k", "dmu"]

    input_files = sys.argv[7:]

    if input_type == "interpro" or input_type == "dmu":
        interproscan_tsv_path = input_files[0]
        annots_by_class = {}
        for rawline in open(interproscan_tsv_path, "r"):
            cells = rawline.strip().split("\t")
            if len(cells) == 2:
                uniprot = cells[0]
                annots = cells[1].split(";")

                for a in annots:
                    if a not in annots_by_class:
                        annots_by_class[a] = set()
                    annots_by_class[a].add(uniprot)
    else:
        raise ValueError(f"Unknown input type: {input_type}")
    # elif input_type == 'taxallnomy':

    print("Total number of classes:", len(annots_by_class))
    print(
        "Total number of proteins:",
        len(set().union(*annots_by_class.values())),
    )
    print("vocab_method", vocab_method)

    n_test_vocabs = 16
    test_sizes = np.linspace(min_vocab_size, max_vocab_size, n_test_vocabs).astype(int)

    test_results = []
    rich_vocab = None
    test_target_ic = target_ic_retained
    bar = tqdm(test_sizes)
    for max_vocab in bar:
        attempt = ICRichVocabulary(
            annots_by_class=annots_by_class,
            target_ic_retained=test_target_ic,
            max_vocab_size=max_vocab,
            enrich_ic=vocab_method in ["ic_rich", "dmu"],
            vocab_method=vocab_method,
        )
        test_results.append(
            {
                "max_vocab_size": max_vocab,
                "vocab_len": len(attempt.vocab),
                "proteins_n": len(attempt.covered),
                "ic_retained": attempt.retained_ic,
                "coverage": len(attempt.covered) / len(attempt.instance_ids),
            }
        )
        if attempt.retained_ic >= target_ic_retained:
            if rich_vocab is None:
                rich_vocab = attempt
                test_target_ic = 0.9995

    cols = ["max_vocab_size", "ic_retained", "coverage", "vocab_len", "proteins_n"]
    with open(output_prefix + ".tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for result in test_results:
            f.write("\t".join([str(result[col]) for col in cols]) + "\n")
    if rich_vocab is not None:
        rich_vocab.save(output_prefix + ".json")
