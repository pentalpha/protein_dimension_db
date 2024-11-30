from collections import Counter
import gzip
from os import path, mkdir
import subprocess
import pandas as pd
from tqdm import tqdm

from gene_ontology import expand_go_set, gos_not_to_use, load_go_graph
from util import (count_lines_large, open_file, run_command, write_file, 
    goa_parsed, goa_parsed_expanded_mf, goa_parsed_expanded_bp, goa_parsed_expanded_cc, 
    config, uniprot_fasta)

#manual urls:
#quickgo:
#https://www.ebi.ac.uk/QuickGO/annotations?aspect=molecular_function&evidenceCode=ECO:0000352,ECO:0000269,ECO:0000314,ECO:0000315,ECO:0000316,ECO:0000353,ECO:0000270,ECO:0007005,ECO:0007001,ECO:0007003,ECO:0007007,ECO:0006056,ECO:0000318,ECO:0000320,ECO:0000321,ECO:0000304,ECO:0000305&evidenceCodeUsage=descendants&withFrom=UniProtKB
#geneontology.org
#http://current.geneontology.org/annotations/filtered_goa_uniprot_all_noiea.gaf.gz


    
if __name__ == '__main__':
    dbs_dir = 'databases'
    if not path.exists(config['go_annotation_raw']):
        download_cmd = ['wget', 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz']
        if not path.exists('goa_uniprot_all.gaf.gz'):
            run_command(download_cmd)
        run_command(['mv', 'goa_uniprot_all.gaf.gz', config['go_annotation_raw']])

    print('Reading swissprot ids')
    uniprots = set()
    for line in open_file(uniprot_fasta):
        if line.startswith('>'):
            uniprots.add(line.rstrip('\n').lstrip('>').split('|')[1])
    if not path.exists(goa_parsed):
        print('Counting length of ', config['go_annotation_raw'])
        count_ann_total_lines = count_lines_large(config['go_annotation_raw'])
        print(count_ann_total_lines, 'lines of annotation')
        
        print('opening', config['go_annotation_raw'])
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
        not_swissprot = 1
        try:
            for line in goa_gaf:
                quickgolines += 1
                if line.startswith('UniProtKB'):
                    cells = line.rstrip('\n').split('\t')
                    if len(cells) >= 13:
                        if cells[1] in uniprots:
                            protid = cells[1]
                            goid = cells[4]
                            evi = cells[6]
                            taxid = cells[12]
                            if len(evi) < 2 or len(evi) > 3:
                                print('strange evidence:', evi)
                            else:
                                if not evi in evi_not_use:
                                    parsed.append('\t'.join([protid,goid,evi,taxid,cells[8]]))
                                else:
                                    droped += 1
                        else:
                            not_swissprot += 1
                    else:
                        incorrect_line_number += 1
                else:
                    other_dbs += 1
                bar.update(1)
        except EOFError as err:
            print('GAF file download incomplete')
            print(err)
        bar.close()
        print(quickgolines, 'lines in goa_uniprot_all original')
        print(other_dbs, 'from other dbs:', other_dbs_names)
        print(other_ontos, 'from other ontologies')
        print(droped, 'with evi codes we cant use')
        print(not_swissprot, 'not swissprot')
        print(incorrect_line_number, 'incorrect_line_number')
        print(quickgolines - other_dbs - other_ontos - droped - incorrect_line_number, len(parsed))
        output_path = goa_parsed
        write_file(output_path).write('\n'.join(parsed))

    parsed = open_file(goa_parsed).read().split('\n')

    print('Loading GO')
    goes_to_not_use = gos_not_to_use()
    go_graph = load_go_graph()
    new_parsed = set()

    print('Expanding GO')
    for rawline in tqdm(parsed[1:]):
        protid, goid, evi, taxids, onto = rawline.split('\t')
        expanded_set = expand_go_set(goid, go_graph, goes_to_not_use)
        for goid2 in expanded_set:
            for taxid in taxids.split('|'):
                new_parsed.add((protid, goid2, evi, taxid, onto))
    
    print(len(parsed), 'annotations from goa_uniprot_all')
    print(len(new_parsed), 'annotations with expansion')
    mf_parsed = ['\t'.join(x).rstrip('\tF') for x in new_parsed if x[-1] == 'F']
    bp_parsed = ['\t'.join(x).rstrip('\tP') for x in new_parsed if x[-1] == 'P']
    cc_parsed = ['\t'.join(x).rstrip('\tC') for x in new_parsed if x[-1] == 'C']
    run_command(['mkdir', path.dirname(goa_parsed_expanded_mf)])
    write_file(goa_parsed_expanded_mf).write('\n'.join(sorted(mf_parsed)))
    write_file(goa_parsed_expanded_bp).write('\n'.join(sorted(bp_parsed)))
    write_file(goa_parsed_expanded_cc).write('\n'.join(sorted(cc_parsed)))