from util_base import open_file, write_file

def filter_fasta(fasta_path, allowed_ids, output_path, id_pos = 0):
    keep = []

    last_title = None
    current_seq = ''
    for rawline in open_file(fasta_path):
        if rawline.startswith('>'):
            if last_title:
                keep.append((last_title, current_seq))
                current_seq = ""
            title_parts = rawline.lstrip('>').rstrip('\n').split('|')
            current_id = title_parts[id_pos]
            if current_id in allowed_ids:
                last_title = current_id
            else:
                last_title = None
        else:
            if last_title:
                current_seq += rawline.rstrip('\n')
    
    output = write_file(output_path)
    for name, seq in keep:
        output.write('>'+name+'\n')
        output.write(seq+'\n')

def remove_from_fasta(fasta_path, to_remove, output_path, id_pos = 0):
    keep = []

    last_title = None
    current_seq = ''
    for rawline in open_file(fasta_path):
        if rawline.startswith('>'):
            if last_title:
                keep.append((last_title, current_seq))
                current_seq = ""
            title_parts = rawline.lstrip('>').rstrip('\n').split('|')
            current_id = title_parts[id_pos]
            if not current_id in to_remove:
                last_title = current_id
            else:
                last_title = None
        else:
            if last_title:
                current_seq += rawline.rstrip('\n')
    
    output = write_file(output_path)
    for name, seq in keep:
        output.write('>'+name+'\n')
        output.write(seq+'\n')

    return [n for n, seq in keep]

def ids_from_fasta(fasta_path):
    ids = []
    taxons = []
    for rawline in open_file(fasta_path):
        if rawline.startswith('>'):
            title_parts = rawline.lstrip('>').rstrip('\n').split('|')
            ids.append(title_parts[0])
            taxids = [x for x in title_parts if 'taxon:' in x]
            taxid = taxons[0] if len(taxids) > 0 else None
            taxons.append(taxid)
    return ids, taxons

def fasta_n_seqs(fasta_path):
    n = 0
    taxons = []
    for rawline in open_file(fasta_path):
        if rawline.startswith('>'):
            n += 1
    return n

def fasta_equal_split(fasta_path, n_fastas: int):
    n_seqs = fasta_n_seqs(fasta_path)
    seqs_per_fasta = int(n_seqs / n_fastas)
    print('Split', fasta_path, 'with', n_seqs, 'sequences',
          'into', n_fastas, 'fastas with', seqs_per_fasta, 'each')
    
    fastas = []
    current_total = 0
    seq = ""
    last_title = None
    current_fasta_path = fasta_path + '.'+str(len(fastas)+1)+'.fasta'
    current_fasta = open(current_fasta_path, 'w')
    for rawline in open_file(fasta_path):
        if rawline.startswith('>'):
            if last_title:
                current_fasta.write('>'+last_title+'\n')
                current_fasta.write(seq+'\n')
                current_total += 1
                if current_total == seqs_per_fasta:
                    current_fasta.close()
                    fastas.append(current_fasta_path)
                    current_fasta_path = fasta_path + '.'+str(len(fastas)+1)+'.fasta'
                    current_fasta = open(current_fasta_path, 'w')
                    current_total = 0
            last_title = rawline.lstrip('>').rstrip('\n')
            seq = ""
        else:
            seq += rawline.rstrip('\n')
    current_fasta.write('>'+last_title+'\n')
    current_fasta.write(seq+'\n')
    current_fasta.close()
    fastas.append(current_fasta_path)

    return fastas