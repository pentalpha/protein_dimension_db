import sys

from bioinfo_utils.fasta import split_fasta_simple, read_fasta, write_fasta

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    max_tokens = int(sys.argv[2])

    output_prefix = sys.argv[3]

    seq_tuples = read_fasta(fasta_path)
    split_seqs = split_fasta_simple(seq_tuples, max_tokens)

    for i, seq_list in enumerate(split_seqs):
        output_path = output_prefix + "." + str(i + 1) + ".fasta"
        write_fasta(output_path, seq_list)
