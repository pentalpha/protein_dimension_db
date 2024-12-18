import sys
from os import path
import gzip
from typing import List
import networkx as nx
import numpy as np
from tqdm import tqdm
from util_base import run_command

prj_dir = path.dirname(path.dirname(__file__))

class TaxaProfileModel():
    def taxallnomy_df_to_nx(taxallnomy_df_path: str) -> nx.DiGraph:
        edges = set()
        nodes = set()
        print('Reading taxallnomy edges')
        edges_file_path = 'taxon_edges.tsv.tmp.gz'
        '''edges_file = gzip.open(edges_file_path, 'wt')
        for rawline in gzip.open(taxallnomy_df_path, 'rt'):
            cells = rawline.rstrip("\n").split('\t')
            taxids = [float(x) for x in cells]
            edges_file.write(str(taxids[0])+'\t'+str(taxids[-1])+'\n')
            #edges.add((taxids[0], taxids[-1]))
            nodes = set()
            for taxid_i in range(2, len(taxids)):
                if taxids[taxid_i] in nodes:
                    break
                edges_file.write(str(taxids[taxid_i])+'\t'+str(taxids[taxid_i-1])+'\n')
                #edges.add((taxids[taxid_i], taxids[taxid_i-1]))
                nodes.add(taxids[taxid_i])
            if not taxids[1] in nodes:
                edges_file.write(str(taxids[1])+'\t'+str(0.0)+'\n')
                #edges.add((taxids[1], 0.0))
                nodes.add(taxids[1])
        edges_file.close()'''

        print('Reading edges temp file')
        edges = []
        for rawline in gzip.open(edges_file_path, 'rt'):
            cells = rawline.rstrip("\n").split('\t')
            if len(cells) == 2:
                edges.append((float(cells[0]), float(cells[1])))

        print('Creating graph')
        G = nx.DiGraph(edges)

        return G
    
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
        self.taxids_to_use = taxids_to_use

        '''print('Preparing taxonomy lists of main taxa')
        self.used_taxids_parents = {}
        self.taxonomy_depth = None
        for taxid in self.taxids_to_use:
            parents_list = list(reversed(nx.shortest_path(self.taxallnomy_graph, taxid, 0.0)))
            self.used_taxids_parents[taxid] = parents_list
            self.taxonomy_depth = len(parents_list)-1'''

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
        if taxid in self.results_cache:
            return self.results_cache[taxid]
        similarities = []
        for common_taxid in self.taxids_to_use:
            similarities.append(1.0 if common_taxid == taxid else 0.0)
        result = np.array(similarities)
        if not taxid in self.results_cache:
            self.results_cache[taxid] = result

        return result

if __name__ == "__main__":
    release_dir = sys.argv[1]

    taxallnomy_df_path = prj_dir + '/databases/taxallnomy.tsv.gz'
    profile_lengths = [128, 256]
    go_experimental_mf_path = release_dir + '/go.experimental.mf.tsv.gz'
    taxid_path = release_dir + '/taxid.tsv'
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
        #taxa_graph_nx = profile_model.taxallnomy_graph
        profiled = [profile_model.calc(taxid) for taxid in tqdm(taxids)]