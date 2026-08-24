#!/usr/bin/env python
import gzip
import sys

import polars as pl

# --
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def consult_ncbi(tax_id=3372387):
    # Get a list of all ancestor TaxIDs
    lineage_ids = ncbi.get_lineage(tax_id)

    # Fetch the names and ranks for those specific IDs
    """names = ncbi.get_taxid_translator(lineage_ids)
    ranks = ncbi.get_rank(lineage_ids)

    print("Detailed Lineage:")
    for tid in lineage_ids:
        rank = ranks.get(tid, "no rank").capitalize()
        name = names.get(tid, "Unknown")
        print(f"{rank}: {name} (ID: {tid})")"""

    return lineage_ids


def read_uniprot_fasta(fasta_path, protein_taxa=None):
    print("Loading", fasta_path)
    if protein_taxa is None:
        protein_taxa = {}

    opener = gzip.open if fasta_path.endswith(".gz") else open

    ids = []

    for line in opener(fasta_path, "rt"):
        if line.startswith(">"):
            header_parts = line.rstrip("\n").lstrip(">").split("|")
            if len(header_parts) == 1:
                uniprot_id = header_parts[0].split()[0]
            else:
                uniprot_id = header_parts[1]
            ids.append(uniprot_id)

            if "OX=" in line:
                taxid = line.split("OX=")[-1].split()[0].rstrip("\n")
            else:
                parts = line.rstrip("\n").split()
                if len(parts) > 1:
                    taxid = parts[-1]
                else:
                    taxid = "0"
                    raise Exception("Taxid not found in line:", line)
            protein_taxa[uniprot_id] = taxid
    return protein_taxa, ids


if __name__ == "__main__":
    # $swissprot_fasta $trembl_fasta
    fasta_path1 = sys.argv[1]
    ids_sorted_path = sys.argv[2]
    taxallnomy_parquet_path = sys.argv[3]
    output_path = sys.argv[4]
    '''fasta_path1 = "uniprot_sprot.fasta.gz"
    ids_sorted_path = "ids.swissprot.txt"
    taxallnomy_parquet_path = "taxallnomy.parquet"
    output_path = "taxid.tsv"'''

    print("Loading fasta", fasta_path1)
    protein_taxa, ids = read_uniprot_fasta(fasta_path1)

    ids_sorted = open(ids_sorted_path, "r").read().strip().split("\n")

    taxa_found = list(int(v) for v in protein_taxa.values())
    print("Loading taxallnomy", taxallnomy_parquet_path)
    taxallnomy_parquet_all = pl.read_parquet(taxallnomy_parquet_path)
    taxallnomy_parquet = taxallnomy_parquet_all.filter(
        pl.col("taxid").is_in(taxa_found)
    )
    lineages = {}
    lineage_columns = [c for c in taxallnomy_parquet.columns if c.startswith("level_")]
    lineage_columns.sort(key=lambda x: int(x.split("_")[1]))

    print("Building lineages")
    for row in taxallnomy_parquet.rows(named=True):
        lineage_int = [int(row[c]) for c in lineage_columns]
        no_repeat = []
        for taxid in lineage_int:
            if taxid not in no_repeat:
                no_repeat.append(taxid)
        lineages[row["taxid"]] = no_repeat

    print("Building tsv")

    output_file = open(output_path, "w")
    columns = ["uniprot_id", "taxid", "lineage"]
    output_file.write("\t".join(columns) + "\n")

    not_found = 0
    for uniprot in ids_sorted:
        taxid = int(protein_taxa[uniprot])
        if taxid in lineages:
            lineage = [str(x) for x in lineages[taxid]]
        else:
            print(taxid, "not found in local taxallnomy, consulting api...")
            lineage = consult_ncbi(tax_id=taxid)
            print("Api result:", lineage)
            lineage = [str(x) for x in lineage]
            not_found += 1

        newline = "\t".join([uniprot, str(taxid), ";".join(lineage)])
        output_file.write(newline + "\n")

    print("not_found", not_found)
    output_file.close()
    print("Done")
