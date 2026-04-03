#!/usr/bin/env python
import gzip
import sys


def filter_by_len(input_fastas, fasta_out_path, ids_path, maxlen):
    half = int(maxlen / 2)
    print("Loading", input_fastas, file=sys.stderr)
    print("Filtering to", fasta_out_path, file=sys.stderr)

    header = None
    content = ""
    sequences = []
    for input_path in input_fastas:
        input = (
            gzip.open(input_path, "rt")
            if input_path.endswith(".gz")
            else open(input_path, "r")
        )
        for line in input:
            if line.startswith(">"):
                if header != None:
                    if len(content) > maxlen:
                        protein_seq = content[:half] + content[-half:]
                    else:
                        protein_seq = content
                    sequences.append((len(content), ">" + header, protein_seq))
                    # output.write('>'+header+'\n')
                    # output.write(content+'\n')
                    content = ""
                header_parts = line.rstrip("\n").lstrip(">").split("|")
                if len(header_parts) > 1:
                    header = header_parts[1]
                else:
                    header = header_parts[0]
            else:
                content += line.rstrip("\n").strip()

        if content != "":
            if len(content) > maxlen:
                protein_seq = content[:half] + content[-half:]
            else:
                protein_seq = content
            sequences.append((len(content), ">" + header, protein_seq))
        input.close()

    sequences.sort()
    output = (
        gzip.open(fasta_out_path, "wt")
        if fasta_out_path.endswith(".gz")
        else open(fasta_out_path, "w")
    )
    for length, header, content in sequences:
        output.write(header + "\n")
        output.write(content + "\n")
    output.close()

    id_list = [header.lstrip(">").split()[0] for length, header, content in sequences]
    with open(ids_path, "w") as id_out:
        id_out.write("\n".join(id_list))


if __name__ == "__main__":
    # Usage: python filter_fasta_by_len.py <input_fasta1> ... <input_fastaN> <output_fasta> <ids_path> <maxlen>
    inputs = sys.argv[1:-3]
    output_fasta = sys.argv[-3]
    ids_path = sys.argv[-2]
    maxlen = int(sys.argv[-1])

    filter_by_len(inputs, output_fasta, ids_path, maxlen)
