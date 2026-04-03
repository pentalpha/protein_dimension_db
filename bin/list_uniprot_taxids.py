#!/usr/bin/env python
import gzip
import sys


def read_uniprot_fasta(fasta_path, protein_taxa=None):
    print("Loading", fasta_path)
    if protein_taxa is None:
        protein_taxa = {}

    opener = gzip.open if fasta_path.endswith(".gz") else open

    for line in opener(fasta_path, "rt"):
        if line.startswith(">"):
            header_parts = line.rstrip("\n").lstrip(">").split("|")
            if len(header_parts) == 1:
                uniprot_id = header_parts[0].split()[0]
            else:
                uniprot_id = header_parts[1]

            if "OX=" in line:
                taxid = line.split("OX=")[-1].split()[0].rstrip("\n")
            else:
                parts = line.rstrip("\n").split()
                if len(parts) > 1:
                    taxid = parts[-1]
                else:
                    taxid = "0"
            protein_taxa[uniprot_id] = taxid
    return protein_taxa


if __name__ == "__main__":
    # $swissprot_fasta $trembl_fasta
    fasta_path1 = sys.argv[1]
    fasta_path2 = sys.argv[2]
    ids_path = sys.argv[3]
    output_path = sys.argv[4]

    uniprots_list = open(ids_path, "r").read().split("\n")
    protein_taxa = read_uniprot_fasta(fasta_path1)
    protein_taxa = read_uniprot_fasta(fasta_path2, protein_taxa)

    protein2taxa = []
    for protein in uniprots_list:
        # print(protein)
        taxon = protein_taxa[protein]
        protein2taxa.append(protein + "\t" + taxon)

    open(output_path, "w").write("\n".join(protein2taxa))
