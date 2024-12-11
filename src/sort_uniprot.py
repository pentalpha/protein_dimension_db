import gzip
import sys

def read_uniprot_fasta(fasta_path):
    print('Loading', fasta_path)
    proteins = []
    header = None
    content = ""
    for line in gzip.open(fasta_path, 'rt'):
        if line.startswith('>'):
            if header != None:
                proteins.append((header, content))
                content = ""
            header = line.rstrip('\n').lstrip('>').split('|')[1]
        else:
            content += line.rstrip('\n').strip()
    
    if content != "":
        proteins.append((header, content))

    return proteins

def sort_fasta(protein_tuples):
    print('Sorting fasta')
    protein_tuples.sort(key = lambda prot: (len(prot[1]), prot[0]))

    ids = [a for a, b in protein_tuples]
    fasta_lines = []
    for header, seq in protein_tuples:
        fasta_lines.append('>'+header)
        fasta_lines.append(seq)
    
    return ids, fasta_lines

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    ids_path = sys.argv[2]
    new_fasta_path = sys.argv[3]
    proteins = read_uniprot_fasta(fasta_path)
    ids, fasta_lines = sort_fasta(proteins)

    open(ids_path, 'w').write('\n'.join(ids))
    gzip.open(new_fasta_path, 'wt').write('\n'.join(fasta_lines))