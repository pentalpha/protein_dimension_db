from sort_uniprot import read_uniprot_fasta
import sys

if __name__ == '__main__':
    fastas = sys.argv[1:]

    for f in fastas:
        print(f)
        fasta_content = read_uniprot_fasta(f)
        names = [n for n, s in fasta_content]
        seqs = [s for n, s in fasta_content]
        for varname, var in [('sequence names', names), ('sequences', seqs)]:
            print('\t'+varname)
            print('\tTotal length:', len(var))
            print('\tUniq values length:', len(set(var)))

        already_seen_seqs = {}
        for n, s in fasta_content:
            if s in already_seen_seqs and len(s) > 30:
                print('Repeated detected:')
                print(already_seen_seqs[s], '\t', s)
                print(n, '\t', s)
                break
            else:
                already_seen_seqs[s] = n