import sys
from os import path
import gzip
from typing import List
import numpy as np
from tqdm import tqdm
from util_base import run_command

prj_dir = path.dirname(path.dirname(__file__))

class TaxaProfileModel():
    
    def __init__(self, taxids_to_use: List[str], taxallnomy_df_path: str = None) -> None:
        print('Starting TaxaProfileModel graph')
        self.taxid_parent_lists = {}
        if taxallnomy_df_path != None and path.exists(taxallnomy_df_path):
            for rawline in gzip.open(taxallnomy_df_path, 'rt'):
                cells = rawline.rstrip("\n").split('\t')
                taxids = [float(x) for x in cells]
                self.taxid_parent_lists[taxids[0]] = taxids
                self.taxonomy_depth = len(taxids)

        self.results_cache = {}
        self.onehot_cache = {}
        self.taxids_to_use = taxids_to_use

    def calc(self, taxid):
        if taxid in self.results_cache:
            return self.results_cache[taxid]
        taxid_parents = self.taxid_parent_lists[taxid]
        similarities = []
        for common_taxid in self.taxids_to_use:
            n_equals = sum([t1 == t2 for t1, t2 in zip(taxid_parents, self.taxid_parent_lists[common_taxid])])
            similarities.append(((n_equals)/self.taxonomy_depth))
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
    
    print('Creating models')
    taxa_graph_nx = None
    for profile_len in profile_lengths:
        taxids_to_use = [t for t, n_gos in most_annotated_taxids][:profile_len]
        profile_model = TaxaProfileModel(taxids_to_use, taxallnomy_df_path=taxallnomy_df_path)
        
        profiled = [profile_model.calc(taxid) for taxid in tqdm(taxids)]
        tops_savepath = release_dir + '/top_taxa_'+str(profile_len)+'.txt'
        open(tops_savepath, 'w').write('\n'.join([str(x) for x  in taxids_to_use]))
        save_path = release_dir + '/emb.taxa_profile_'+str(profile_len)+'.npy.gz'
        np.save(save_path.rstrip('.gz'), np.asarray(profiled))
        if path.exists(save_path):
            run_command(['rm', save_path])
        run_command(['gzip', save_path.rstrip('.gz')])
        
        onehot = [profile_model.calc_onehot(taxid) for taxid in tqdm(taxids)]
        save_path2 = release_dir + '/onehot.taxa_'+str(profile_len)+'.npy.gz'
        np.save(save_path2.rstrip('.gz'), np.asarray(onehot))
        if path.exists(save_path2):
            run_command(['rm', save_path2])
        run_command(['gzip', save_path2.rstrip('.gz')])
        del profile_model