from glob import glob
import gzip
from os import path
import subprocess
from tqdm import tqdm
import yaml
import numpy as np

from util_base import open_file, proj_dir, write_file

config = yaml.safe_load(open(proj_dir+"/config.yml", "r"))

release_dir = config['release_dir']

uniprot_fasta = release_dir+"/databases/uniprot_sprot.fasta.gz"
goa_parsed = proj_dir+"/databases/goa_parsed.tsv.gz"
goa_parsed_expanded_mf = release_dir+"/databases/goa_parsed_expanded.mf.tsv.gz"
goa_parsed_expanded_bp = release_dir+"/databases/goa_parsed_expanded.bp.tsv.gz"
goa_parsed_expanded_cc = release_dir+"/databases/goa_parsed_expanded.cc.tsv.gz"
annotation_path = release_dir+'/annotation.tsv'
taxon_profile_path = release_dir+'/taxa_profile.tsv.gz'
esm_features_prefix = release_dir+'/esm2_t'
esm_features_pattern = release_dir+'/esm2_t*.tsv.gz'
esm_model_ids = [str(x) for x in config['esm_models_to_use']]
esm_features_paths = [esm_features_prefix+x+'.tsv.gz' for x in esm_model_ids]
labels_path = release_dir+'/go_labels.tsv'
ids_path = release_dir+'/ids.txt'

go_not_use_path = proj_dir+"/databases/gocheck_do_not_annotate.json"
go_basic_path = proj_dir+"/databases/go-basic.obo"

def load_parsed_goa(file_path=goa_parsed):
    anns = []
    for rawline in open_file(file_path):
        protid, goid, evi, taxid = rawline.rstrip('\n').split('\t')
        if goid.startswith('GO') and taxid.startswith('tax'):
            anns.append([protid, goid, evi, taxid])

    return anns

def load_label_lists(uniprotids: list):
    label_lists = {}
    for rawline in open_file(labels_path):
        uniprotid, taxid, gos = rawline.rstrip('\n').split('\t')
        if uniprotid in uniprotids:
            label_lists[uniprotid] = gos.split(',')
    label_lists2 = [label_lists[uniprotid] for uniprotid in uniprotids]
    return label_lists2

def load_features_from_dir(dirname: str, ids_allowed: list = []):
    features_path = dirname+'/features.npy'
    ids_path = dirname+'/ids.txt'

    features = np.load(features_path)
    ids = open(ids_path, 'r').read().split('\n')

    protein_indexes = {ids[i]: i for i in range(len(ids))}

    new_index = {}
    local_features = []
    i = 0
    print('Selecting features of proteins')
    for protein in tqdm(ids_allowed):
        feature_index = protein_indexes[protein]
        feature_vec = features[feature_index]
        local_features.append(feature_vec)
        new_index[protein] = i
        i+= 1
    print('Converting')
    local_features = np.asarray(local_features)

    return local_features, new_index

features_cache = {}

def load_features(feature_file_path: str, subset: list, converter):
    ids_path = path.dirname(feature_file_path)+'/ids.txt'

    ids = open(ids_path, 'r').read().split('\n')
    features = []
    id_to_line = {}
    line_i = 0
    for protid in subset:
        features.append(None)
        id_to_line[protid] = line_i
        line_i += 1

    bar = tqdm(total=len(subset))
    '''if not feature_file_path in features_cache:'''
    print('Loading', feature_file_path)
    feature_file = open_file(feature_file_path)
    lines = []
    for rawline in feature_file:
        cols = rawline.rstrip('\n').split('\t')
        vals = [converter(x) for x in cols]
        lines.append(vals)
    #features_cache[feature_file_path] = lines
    print('Loaded', feature_file_path)
    feature_file.close()

    #lines = features_cache[feature_file_path]
    line_i = 0
    for vals in lines:
        current_id = ids[line_i]
        if current_id in subset:
            correct_line = id_to_line[current_id]
            features[correct_line] = vals
            bar.update(1)
        line_i += 1
    bar.close()
    assert not (None in features)
    return features


def filter_features(feature_file_path: str, subset: list, new_feature_file_path: str):
    ids_path = path.dirname(feature_file_path)+'/ids.txt'

    ids = open(ids_path, 'r').read().split('\n')
    feature_file = open_file(feature_file_path)
    
    features = []
    id_to_line = {}
    line_i = 0
    for protid in subset:
        features.append(None)
        id_to_line[protid] = line_i
        line_i += 1

    bar = tqdm(total=len(subset))
    line_i = 0
    while line_i < len(ids):
        current_line = feature_file.readline()
        current_id = ids[line_i]
        if current_id in subset:
            correct_line = id_to_line[current_id]
            features[correct_line] = current_line.rstrip('\n')
            bar.update(1)
        line_i += 1
    bar.close()
    assert not (None in features)

    feature_file2 = write_file(new_feature_file_path)
    feature_file2.write('\n'.join(features))
    return features

def load_labels_from_dir(dirname: str, ids_allowed: list = []):
    labels_path = dirname+'/go_labels.tsv'

    labels = {}
    anns = {}
    bar = tqdm(total=len(ids_allowed))
    for rawline in open(labels_path, 'r'):
        uniprotid, taxonid, gos = rawline.rstrip('\n').split('\t')
        protid = uniprotid+'\t'+taxonid
        if protid in ids_allowed:
            go_list = gos.split(',')
            labels[protid] = go_list
            for go in go_list:
                if not go in anns:
                    anns[go] = set()
                anns[go].add(protid)
            bar.update(1)
    bar.close()

    return labels, anns

def load_dataset_from_dir(subset: list = [], to_load = ['taxa_profile', 'esm']):
    if len(subset) == 0:
        subset = open(ids_path, 'r').read().split('\n')
    labels, annotations = load_labels_from_dir(release_dir, ids_allowed=subset)

    features = []

    if 'taxa_profile' in to_load:
        taxa_profile_path = taxon_profile_path
        taxa_profile_features = load_features(taxa_profile_path, subset, float)
        features.append(('taxa_profile', taxa_profile_features))

    if 'esm' in to_load:
        esm_paths = glob(esm_features_pattern)
        for esm_path in esm_paths:
            esm_len_str = esm_path.split('.')[-3].split('_')[-1]
            esm_len = int(esm_len_str)
            print('loading', esm_path)
            esm_features = load_features(esm_path, subset, float)
            features.append(('esm_'+str(esm_len), esm_features))
    
    return features, labels, annotations