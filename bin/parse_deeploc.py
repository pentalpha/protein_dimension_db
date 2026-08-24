#!/usr/bin/env python
import sys
import pandas as pd
import polars as pl

deeploc_path = sys.argv[1]
deeploc_membrane_path = sys.argv[2]
output_parquet_path = sys.argv[3]

by_uniprot = {}
classes_a = [
    "Membrane",
    "Cytoplasm",
    "Nucleus",
    "Extracellular",
    "Cell membrane",
    "Mitochondrion",
    "Plastid",
    "Endoplasmic reticulum",
    "Lysosome/Vacuole",
    "Golgi apparatus",
    "Peroxisome",
]
classes_membrane = ["Peripheral", "Transmembrane", "LipidAnchor", "Soluble"]
for _, row in pd.read_csv(deeploc_path).iterrows():
    uniprot = row["ACC"]
    by_uniprot[uniprot] = {k: True if row[k] == 1 else False for k in classes_a}
    for k in classes_membrane:
        by_uniprot[uniprot][k] = False

for _, row in pd.read_csv(deeploc_membrane_path).iterrows():
    uniprot = row["ACC"]
    if not uniprot in by_uniprot:
        by_uniprot[uniprot] = {k: False for k in classes_a}
    for k in classes_membrane:
        by_uniprot[uniprot][k] = True if row[k] == 1 else False

rows = []
for uniprot_id, annots in by_uniprot.items():
    new_row = {"id": uniprot_id}
    for k, v in annots.items():
        new_row[k] = v
    rows.append(new_row)

df = pl.DataFrame(rows)
df.write_parquet(output_parquet_path)

for col in classes_a + classes_membrane:
    counts = df[col].sum()
    print(f"class {col}: {counts} positives")
