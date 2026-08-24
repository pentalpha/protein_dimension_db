#!/usr/bin/env python3

import sys
import gzip
import os
from typing import Counter

from tqdm import tqdm
import numpy as np
import polars as pl

input_tsv = sys.argv[1]
final_parquet = sys.argv[2]


def load_taxallnomy_parquet():
    if os.path.exists(final_parquet):
        print("Loading parquet file")
        return pl.read_parquet(final_parquet)
    # rows are species, columns are taxonomy levels
    # last columns are more specific taxa levels, first columns are more general taxa levels
    # values are ncbi taxon ids
    p = input_tsv
    max_taxa = 2718183

    print("Loading tsv into numpy matrix")
    # load tsv into numpy matrix
    with gzip.open(p, "rt") as f:
        rows = []
        bar = tqdm(total=max_taxa)
        for line in f:
            # print(line)
            try:
                new_row = [float(x) for x in line.strip().split("\t")]
                rows.append(new_row)
                bar.update(1)
            except ValueError as err:
                print(f"Error parsing line: {line}")
                print(err)
                quit(1)
        bar.close()

    n_cols = len(rows[0])
    first_col_name = "taxid"
    other_col_names = [f"level_{i+1}" for i in range(n_cols - 1)]

    print("Parsing rows into columns")
    # separate by columns
    columns = {first_col_name: []}
    for i in range(n_cols - 1):
        col_name = other_col_names[i]
        columns[col_name] = []

    for row in tqdm(rows):
        columns[first_col_name].append(row[0])
        for i in range(n_cols - 1):
            col_name = other_col_names[i]
            columns[col_name].append(row[i + 1])

    columns["taxid"] = [int(x) for x in columns["taxid"]]
    taxallnomy = pl.DataFrame(columns)
    print(taxallnomy)
    print(taxallnomy.shape)

    print("Saving to file")
    taxallnomy.write_parquet(final_parquet)


df = load_taxallnomy_parquet()
"""# convert all columns to int
df = df.select(pl.all().cast(pl.Int64))

print(df.head())

vocab_size = 256

for col in df.columns:
    # Count unique values
    print(col, df[col].n_unique())

# First 12, except for first one
generic_vocab_cols = df.columns[1:13]
print(generic_vocab_cols)

# Last 12
specific_vocab_cols = df.columns[-12:]
print(specific_vocab_cols)"""

"""specific_counts = {}
for col in specific_vocab_cols:
    counts = Counter(df[col].to_list())
    for taxid, count in counts.items():
        if taxid not in specific_counts:
            specific_counts[taxid] = 0
        specific_counts[taxid] += count

specific_counts_vec = sorted(
    [(count, taxid) for taxid, count in specific_counts.items()], reverse=True
)

specific_vocab_size = vocab_size // 2
specific_vocab = [taxid for count, taxid in specific_counts_vec[:specific_vocab_size]]
print("Top in specific_counts_vec (count, taxid):")
print(specific_counts_vec[:10])
print("Bottom in specific_counts_vec (count, taxid):")
print(specific_counts_vec[-10:])"""
# TODO: count taxid frequencies using taxids of uniprot ids
