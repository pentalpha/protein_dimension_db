#!/usr/bin/env python
import polars as pl
import sys

from tqdm import tqdm

if __name__ == "__main__":
    go_by_protein_path = sys.argv[1]
    trembl_ids_path = "./trembl_ids.txt"

    print("Reading", go_by_protein_path)
    go_by_protein = pl.read_parquet(go_by_protein_path)

    trembl_proteins = set()
    used_terms = []
    for row in tqdm(go_by_protein.rows(named=True), total=go_by_protein.height):
        if row["ProteinSet"] == "TrEMBL":
            trembl_proteins.add(row["id"])
    trembl_proteins = sorted(trembl_proteins)
    print(len(trembl_proteins), "TrEMBL proteins")
    open(trembl_ids_path, "w").write("\n".join(trembl_proteins))
