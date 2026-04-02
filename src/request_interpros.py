import sys
from data.interpro_api.consult_interpro import consult_uniprot_ids
from bioinfo_utils.fasta import read_fasta
from os import path, remove

if __name__ == "__main__":
    # input fasta or txt file with uniprot ids
    input_file = sys.argv[1]
    # tsv file to save annotations
    output_file = sys.argv[2]

    results_cache = output_file + ".cache"

    if path.exists(results_cache + ".lock"):
        remove(results_cache + ".lock")

    fasta_seqs = read_fasta(input_file, type="uniprot")
    ids = [seq[0] for seq in fasta_seqs]
    consult_uniprot_ids(ids, results_cache, info_type="family")

    # Copy cache to final output
    with open(results_cache, "r") as cache_file, open(output_file, "w") as out_file:
        for line in cache_file:
            out_file.write(line)
