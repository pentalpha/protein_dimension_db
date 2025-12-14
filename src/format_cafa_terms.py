import sys
import gzip
from collections import defaultdict

def format_cafa_terms(input_path, output_prefix):
    # EntryID	term	aspect
    # Q5W0B1	GO:0000785	C
    
    # Store annotations: aspect -> entry_id -> list of terms
    annotations = {
        'F': defaultdict(set), # MF
        'P': defaultdict(set), # BP
        'C': defaultdict(set)  # CC
    }
    
    print(f"Reading {input_path}")
    opener = gzip.open if input_path.endswith('.gz') else open
    mode = 'rt' if input_path.endswith('.gz') else 'r'
    
    with opener(input_path, mode) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            entry_id, term, aspect = parts
            if aspect in annotations:
                annotations[aspect][entry_id].add(term)
    
    # Mapping aspect char to output aspect string and filename suffix
    aspect_map = {
        'F': 'mf',
        'P': 'bp',
        'C': 'cc'
    }
    
    for aspect_char, aspect_name in aspect_map.items():
        output_file = f"{output_prefix}.{aspect_name}.tsv"
        print(f"Writing {output_file}")
        with open(output_file, 'w') as f:
            for entry_id, terms in annotations[aspect_char].items():
                terms_str = ','.join(sorted(list(terms)))
                f.write(f"{entry_id}\t{terms_str}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python format_cafa_terms.py <input_train_terms> <output_prefix>")
        sys.exit(1)
        
    input_terms = sys.argv[1]
    output_prefix = sys.argv[2] # e.g. "go.experimental"
    
    format_cafa_terms(input_terms, output_prefix)
