import gzip
from os import path
import sys
import pandas as pd
from tqdm import tqdm

from gene_ontology import expand_go_set, gos_not_to_use, load_go_graph
from util_base import count_lines_large, open_file, write_file

#manual urls:
#quickgo:
#https://www.ebi.ac.uk/QuickGO/annotations?aspect=molecular_function&evidenceCode=ECO:0000352,ECO:0000269,ECO:0000314,ECO:0000315,ECO:0000316,ECO:0000353,ECO:0000270,ECO:0007005,ECO:0007001,ECO:0007003,ECO:0007007,ECO:0006056,ECO:0000318,ECO:0000320,ECO:0000321,ECO:0000304,ECO:0000305&evidenceCodeUsage=descendants&withFrom=UniProtKB
#geneontology.org
#http://current.geneontology.org/annotations/filtered_goa_uniprot_all_noiea.gaf.gz

def gafparsed_to_id2go(gaf_info, uniprot_ids, output_file):
    go_lists = {x: set() for x in uniprot_ids}
    for protid, goid2, evi, taxid, onto in gaf_info:
        if protid in go_lists:
            go_lists[protid].add(goid2)
    
    id2go = [x+'\t'+','.join(go_lists[x]) for x in uniprot_ids]
    gzip.open(output_file, 'wt').write('\n'.join(id2go))

def parse_gaf(goa_gaf, evi_not_use, bar):
    parsed = ['prot_id\tgoid\tevi\ttaxonid']
    droped = 0
    other_ontos = 0
    quickgolines = 0
    other_dbs = 0
    incorrect_line_number = 0
    other_dbs_names = set()
    not_swissprot = 1
    #cols used: 1, protein_id, go_id, evi_type, taxid, 8
    try:
        for line in goa_gaf:
            quickgolines += 1
            if line.startswith('UniProtKB'):
                cells = line.rstrip('\n').split('\t')
                if len(cells) >= 13:
                    if cells[1] in uniprots_set:
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

    return parsed



if __name__ == '__main__':
    evi_not_use_path = sys.argv[1]
    go_not_use_path = sys.argv[2]
    go_basic_path = sys.argv[3]
    go_annotation_raw = sys.argv[4]
    ids_path = sys.argv[5]
    output_dir = './'
    calc_go_expanded = sys.argv[6] == "rerun"

    #evi_not_use_path = proj_dir+'/evi_not_to_use.txt'
    #go_not_use_path = proj_dir+"/databases/gocheck_do_not_annotate.json"
    #go_basic_path = proj_dir+"/databases/go-basic.obo"
    
    goa_parsed = output_dir+'/go.experimental.tsv.gz'
    goa_expanded = output_dir+'/go.expanded.tsv.gz'
    
    goa_parsed_expanded_mf = output_dir+'/go.experimental.mf.tsv.gz'
    goa_parsed_expanded_bp = output_dir+'/go.experimental.bp.tsv.gz'
    goa_parsed_expanded_cc = output_dir+'/go.experimental.cc.tsv.gz'

    print('Reading swissprot ids')
    uniprots_list = open(ids_path).read().split('\n')

    uniprots_set = set(uniprots_list)
    
    if calc_go_expanded or not path.exists(goa_expanded):
        print('Counting length of ', go_annotation_raw)
        count_ann_total_lines = count_lines_large(go_annotation_raw)
        print(count_ann_total_lines, 'lines of annotation')
        
        print('opening', go_annotation_raw)
        bar = tqdm(total = count_ann_total_lines)
        goa_gaf = open_file(go_annotation_raw)
        evi_df = pd.read_csv(evi_not_use_path,sep=',')
        evi_not_use = set(evi_df['code'].tolist())

        parsed = parse_gaf(goa_gaf, evi_not_use, bar)
        output_path = goa_parsed
        write_file(output_path).write('\n'.join(parsed))

        parsed = open_file(goa_parsed).read().split('\n')

        print('Loading GO')
        goes_to_not_use = gos_not_to_use(go_not_use_path)
        go_graph = load_go_graph(go_basic_path)
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
        parsed_all = ['\t'.join(x) for x in new_parsed]
        write_file(goa_expanded).write('\n'.join(sorted(parsed_all)))
    
    print('Reading', goa_expanded)
    go_expanded_ann = [l.split('\t') for l in gzip.open(goa_expanded, 'rt').read().split('\n')]

    mf_lines = [x for x in go_expanded_ann if x[-1] == 'F']
    bp_lines = [x for x in go_expanded_ann if x[-1] == 'P']
    cc_lines = [x for x in go_expanded_ann if x[-1] == 'C']

    print('Writing id2go files')
    gafparsed_to_id2go(mf_lines, uniprots_list, goa_parsed_expanded_mf)
    gafparsed_to_id2go(bp_lines, uniprots_list, goa_parsed_expanded_bp)
    gafparsed_to_id2go(cc_lines, uniprots_list, goa_parsed_expanded_cc)
    
    '''mf_parsed = ['\t'.join(x).rstrip('\tF') for x in new_parsed if x[-1] == 'F']
    bp_parsed = ['\t'.join(x).rstrip('\tP') for x in new_parsed if x[-1] == 'P']
    cc_parsed = ['\t'.join(x).rstrip('\tC') for x in new_parsed if x[-1] == 'C']
    run_command(['mkdir', path.dirname(goa_parsed_expanded_mf)])
    write_file(goa_parsed_expanded_mf).write('\n'.join(sorted(mf_parsed)))
    write_file(goa_parsed_expanded_bp).write('\n'.join(sorted(bp_parsed)))
    write_file(goa_parsed_expanded_cc).write('\n'.join(sorted(cc_parsed)))'''