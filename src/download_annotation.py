from collections import Counter
import gzip
from os import path, mkdir
import subprocess
import pandas as pd
from tqdm import tqdm

from gene_ontology import expand_go_set, gos_not_to_use, load_go_graph
from util import (load_parsed_goa, open_file, run_command, write_file, 
    goa_parsed, goa_parsed_expanded, goa_parsed_frequent, 
    config, count_lines, write_parsed_goa)

#manual urls:
#quickgo:
#https://www.ebi.ac.uk/QuickGO/annotations?aspect=molecular_function&evidenceCode=ECO:0000352,ECO:0000269,ECO:0000314,ECO:0000315,ECO:0000316,ECO:0000353,ECO:0000270,ECO:0007005,ECO:0007001,ECO:0007003,ECO:0007007,ECO:0006056,ECO:0000318,ECO:0000320,ECO:0000321,ECO:0000304,ECO:0000305&evidenceCodeUsage=descendants&withFrom=UniProtKB
#geneontology.org
#http://current.geneontology.org/annotations/filtered_goa_uniprot_all_noiea.gaf.gz


    
if __name__ == '__main__':
    dbs_dir = 'databases'
    if not path.exists(config['go_annotation_raw']):
        download_cmd = ['wget', 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz']
        run_command(download_cmd)
        run_command(['mv', 'goa_uniprot_all.gaf.gz', config['go_annotation_raw']])

    '''print('Counting length of ', config['go_annotation_raw'])
    count_ann_total_lines = count_lines(config['go_annotation_raw'])
    
    bar = tqdm(total = count_ann_total_lines)
    goa_gaf = open_file(config['go_annotation_raw'])
    evi_df = pd.read_csv('evi_not_to_use.txt',sep=',')
    evi_not_use = set(evi_df['code'].tolist())

    parsed = ['prot_id\tgoid\tevi\ttaxonid']
    droped = 0
    other_ontos = 0
    quickgolines = 0
    other_dbs = 0
    incorrect_line_number = 0
    other_dbs_names = set()
    try:
        for line in goa_gaf:
            quickgolines += 1
            cells = line.rstrip('\n').split('\t')
            if len(cells) >= 13:
                if cells[0] == 'UniProtKB':
                    if cells[8] == 'F':
                        protid = cells[1]
                        goid = cells[4]
                        evi = cells[6]
                        taxid = cells[12]
                        if len(evi) < 2 or len(evi) > 3:
                            print('strange evidence:', evi)
                        else:
                            if not evi in evi_not_use:
                                parsed.append('\t'.join([protid,goid,evi,taxid]))
                            else:
                                droped += 1
                    else:
                        other_ontos += 1
                else:
                    other_dbs += 1
                    other_dbs_names.add(line.split('\t')[0])
            else:
                incorrect_line_number += 1
            
            bar.update(1)
    except EOFError as err:
        print('GAF file download incomplete')
        print(err)
    bar.close()
    print(quickgolines, 'lines in quickgo original')
    print(other_dbs, 'from other dbs:', other_dbs_names)
    print(other_ontos, 'from other ontologies')
    print(droped, 'with evi codes we cant use')
    print(incorrect_line_number, 'incorrect_line_number')
    print(quickgolines - other_dbs - other_ontos - droped - incorrect_line_number, len(parsed))
    output_path = goa_parsed
    write_file(output_path).write('\n'.join(parsed))'''

    parsed = open_file(goa_parsed).read().split('\n')

    print('Loading GO')
    goes_to_not_use = gos_not_to_use()
    go_graph = load_go_graph()

    already_added = set()

    new_parsed = []

    print('Expanding GO')
    for rawline in tqdm(parsed[1:]):
        protid, goid, evi, taxids = rawline.split('\t')
        expanded_set = expand_go_set(goid, go_graph, goes_to_not_use)
        for goid2 in expanded_set:
            if not (protid, goid2) in already_added:
                for taxid in taxids.split('|'):
                    new_parsed.append('\t'.join([protid, goid2, evi, taxid]))
                    already_added.add((protid, goid2))
    
    print(len(parsed), 'annotations from goa_uniprot_all')
    print(len(new_parsed), 'annotations with expansion')
    run_command(['mkdir', 'input'])
    write_file(goa_parsed_expanded).write('\n'.join(new_parsed))

    #new_parsed = open_file(goa_parsed_expanded).read().split('\n')
    all_annotations = load_parsed_goa(file_path=goa_parsed_expanded)
    
    all_goids = [goid for protid, goid, evi, taxid in all_annotations]
    print(len(set(all_goids)), 'GO IDs from uniprot proteins in', goa_parsed_expanded)
    go_counts = Counter(all_goids)
    frequent_go_ids = set()
    for goid, freq in go_counts.items():
        if freq >= config['min_annotations']:
            frequent_go_ids.add(goid)
    print(len(frequent_go_ids), 'frequent GO IDs in', goa_parsed_expanded)

    all_annotations = [x for x in all_annotations if x[1] in frequent_go_ids]
    write_parsed_goa(all_annotations, goa_parsed_frequent)