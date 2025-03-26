import gzip
import sys

def read_uniprot_fasta(fasta_path):
    print('Loading', fasta_path)
    protein_taxa = {}
    for line in gzip.open(fasta_path, 'rt'):
        if line.startswith('>'):
            header_parts = line.rstrip('\n').lstrip('>').split('|')
            if len(header_parts) == 1:
                uniprot_id = header_parts[0]
            else:
                uniprot_id = header_parts[1]
            taxid = line.split('OX=')[-1].split()[0].rstrip('\n')
            protein_taxa[uniprot_id] = taxid
    return protein_taxa

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    ids_path = sys.argv[2]
    output_path = sys.argv[3]

    uniprots_list = open(ids_path, 'r').read().split('\n')
    protein_taxa = read_uniprot_fasta(fasta_path)

    protein2taxa = []
    for protein in uniprots_list:
        #print(protein)
        taxon = protein_taxa[protein]
        protein2taxa.append(protein+'\t'+taxon)

    open(output_path, 'w').write('\n'.join(protein2taxa))