import gzip
import sys

def filter_by_len(fasta_path, fasta_out_path, maxlen):
    print('Loading', fasta_path, file=sys.stderr)
    print('Filtering to', fasta_out_path, file=sys.stderr)
    
    header = None
    content = ""
    input = gzip.open(fasta_path, 'rt') if fasta_path.endswith('.gz') else open(fasta_path, 'r')
    output = gzip.open(fasta_out_path, 'wt') if fasta_out_path.endswith('.gz') else open(fasta_out_path, 'w')
    for line in input:
        if line.startswith('>'):
            if header != None:
                if len(content) <= maxlen:
                    output.write('>'+header+'\n')
                    output.write(content+'\n')
                content = ""
            header_parts = line.rstrip('\n').lstrip('>').split('|')
            if len(header_parts) > 1:
                header = header_parts[1]
            else:
                header = header_parts[0]
        else:
            content += line.rstrip('\n').strip()
    
    if content != "":
        if len(content) <= maxlen:
            output.write('>'+header+'\n')
            output.write(content+'\n')
    
    input.close()
    output.close()

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    fasta_path2 = sys.argv[2]
    maxlen = int(sys.argv[3])
    filter_by_len(fasta_path, fasta_path2, maxlen)