import sys
from os import path
import gzip
from collections import Counter
from typing import List
import numpy as np
from tqdm import tqdm
from glob import glob
import polars as pl
from util_base import run_command

#prj_dir = path.dirname(path.dirname(__file__))

class TaxaProfileModel():
    
    def __init__(self, taxids_to_use: List[str], taxallnomy_df_path: str = None) -> None:
        print('Starting TaxaProfileModel graph')
        self.taxid_parent_lists = {}

        assert taxallnomy_df_path != None and path.exists(taxallnomy_df_path)

        for rawline in gzip.open(taxallnomy_df_path, 'rt'):
            cells = rawline.rstrip("\n").split('\t')
            taxids = [float(x) for x in cells]
            parent_taxa = list(reversed(taxids[1:]))
            self.taxid_parent_lists[taxids[0]] = parent_taxa
        
        for vec in self.taxid_parent_lists.values():
            self.taxonomy_depth = len(vec)
            break

        self.results_cache = {}
        self.onehot_cache = {}
        self.taxids_to_use = taxids_to_use
    
    '''
    Find frequency of taxids missing from taxid_parent_lists
    '''
    def find_missing_taxids(self, taxids):
        taxid_frequency = {k: c for k, c in Counter(taxids).items()}
        missing_taxids = [k for k, c in taxid_frequency.items() if k not in self.taxid_parent_lists]
        total_missing = sum(taxid_frequency[k] for k in missing_taxids)
        freq_relative = total_missing / len(taxids)
        print('Taxids missing from taxalnomy:', len(missing_taxids), 'total missing:', total_missing, file=sys.stderr)
        print('Total missing percentage:', freq_relative, file=sys.stderr)

        for missing in missing_taxids:
            self.results_cache[missing] = np.zeros(len(self.taxids_to_use))

    def calc(self, taxid):
        if taxid in self.results_cache:
            return self.results_cache[taxid]
        
        if not taxid in self.taxid_parent_lists:
            print(taxid, 'not found in parent lists', file=sys.stderr)
            print(len(self.taxid_parent_lists), 'length of parent lists', file=sys.stderr)
            print('some valid keys:', file=sys.stderr)
            m = 3
            for key in self.taxid_parent_lists.keys():
                print(key, type(key))
                m -= 1
                if m == 0:
                    break
            quit(1)
        
        taxid_parents = self.taxid_parent_lists[taxid]
        similarities = []
        for common_taxid in self.taxids_to_use:
            if common_taxid == taxid:
                similarities.append(1.0)
            else:
                n_equals = sum([t1 == t2 for t1, t2 in zip(taxid_parents, self.taxid_parent_lists[common_taxid])])
                similarities.append((n_equals/self.taxonomy_depth)*0.95)
        result = np.array(similarities)
        if not taxid in self.results_cache:
            self.results_cache[taxid] = result

        return result
    
    def calc_onehot(self, taxid):
        if taxid in self.onehot_cache:
            return self.onehot_cache[taxid]
        similarities = []
        for common_taxid in self.taxids_to_use:
            similarities.append(1.0 if common_taxid == taxid else 0.0)
        result = np.array(similarities)
        if not taxid in self.onehot_cache:
            self.onehot_cache[taxid] = result

        return result

if __name__ == "__main__":
    release_dir = '.'
    taxallnomy_df_path = sys.argv[1]
    #taxallnomy_df_path = prj_dir + '/databases/taxallnomy.tsv.gz'
    profile_lengths = [128, 256]
    #go_experimental_mf_path = release_dir + '/go.experimental.mf.tsv.gz'
    go_experimental_mf_path = sys.argv[2]
    #taxid_path = release_dir + '/taxid.tsv'
    taxid_path = sys.argv[3]

    print(sys.argv)
    print([path.exists(a) for a in sys.argv])
    
    # Check if files exist before proceeding
    if not path.exists(go_experimental_mf_path):
        print(f"Error: File {go_experimental_mf_path} does not exist", file=sys.stderr)
        sys.exit(1)
    if not path.exists(taxid_path):
        print(f"Error: File {taxid_path} does not exist", file=sys.stderr)
        sys.exit(1)
    
    print('Loading taxon of uniprotids')
    uniprot2taxid = {}
    uniprots = []
    taxids = []
    for rawline in open(taxid_path, 'r'):
        uniprotid, taxid = rawline.rstrip("\n").split('\t')
        uniprot2taxid[uniprotid] = float(taxid)
        uniprots.append(uniprotid)
        taxids.append(uniprot2taxid[uniprotid])
    
    print('Counting n annotations of taxa')
    gos_by_species = {t: 0 for t in uniprot2taxid.values()}
    for rawline in gzip.open(go_experimental_mf_path, 'rt'):
        cells = rawline.rstrip("\n").split('\t')
        if len(cells) == 2:
            uniprot_id = cells[0]
            taxid = uniprot2taxid[uniprot_id]
            gos_by_species[taxid] += len(cells[1].split(','))

    print('Sorting according to number of annotations')
    most_annotated_taxids = [(t, n_gos) for t, n_gos in gos_by_species.items()]
    most_annotated_taxids.sort(key = lambda tp: -tp[1])
    most_annotated_taxids = most_annotated_taxids[:max(profile_lengths)]
    for x, y in most_annotated_taxids:
        print(x, y)
    
    print('Creating models')
    taxa_graph_nx = None
    for profile_len in profile_lengths:
        if len(most_annotated_taxids) > profile_len:
            taxids_to_use = [t for t, n_gos in most_annotated_taxids][:profile_len]
        else:
            taxids_to_use = [t for t, n_gos in most_annotated_taxids]

        tops_savepath = release_dir + '/top_taxa_'+str(profile_len)+'.txt'
        open(tops_savepath, 'w').write('\n'.join([str(x) for x  in taxids_to_use]))

        profile_model = TaxaProfileModel(taxids_to_use, taxallnomy_df_path=taxallnomy_df_path)
        profile_model.find_missing_taxids(taxids)
        #quit(1)
        profiled = np.asarray([profile_model.calc(taxid) for taxid in tqdm(taxids)])
        save_path = release_dir + '/emb.taxa_profile_'+str(profile_len)+'.parquet'
        pl.DataFrame({'id': uniprots, 'emb': profiled}).write_parquet(save_path)
        
        onehot = np.asarray([profile_model.calc_onehot(taxid) for taxid in tqdm(taxids)])
        save_path2 = release_dir + '/onehot.taxa_'+str(profile_len)+'.parquet'
        pl.DataFrame({'id': uniprots, 'emb': onehot}).write_parquet(save_path2)
        
        del profile_model