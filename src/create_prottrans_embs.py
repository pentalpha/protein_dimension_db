import gzip
import sys
import h5py
import numpy as np
from tqdm import tqdm
import pandas as pd

if __name__ == '__main__':
    proj_dir = sys.argv[1]
    output_dir = sys.argv[2]

    databases_dir = proj_dir + "/databases"
    original_embs_path = databases_dir + "/per-protein.h5"
    ids_path = output_dir + "/databases/ids.txt"

    ids_list = open(ids_path, 'r').read().split('\n')
    print('Ids from swiss prot:', len(ids_list))
    id_to_index = {x: i for i, x in enumerate(ids_list)}
    sorted_embs = [None for _ in range(len(ids_list))]

    file = h5py.File(original_embs_path, "r")
    print('ProtTrans-T5 embeddings from uniprot:', len(file))
    indexes_used = set()
    for uniprotid, rawemb in tqdm(file.items()):
        if uniprotid in id_to_index:
            index = id_to_index[uniprotid]
            sorted_embs[index] = np.array(rawemb)
            indexes_used.add(index)
    
    print('Replacing none with nan')
    emb_shape = None
    for i in range(len(sorted_embs)):
        if i in indexes_used:
            emb_shape = sorted_embs[i].shape
            break

    for i in tqdm(range(len(sorted_embs))):
        if not i in indexes_used:
            sorted_embs[i] = np.empty(emb_shape)
            sorted_embs[i].fill(np.nan)
    
    print('Creating numpy matrix')
    sorted_embs = np.asarray(sorted_embs)
    print('Saving in npy.gz')
    np.save(gzip.open(output_dir+'/prottrans.npy.gz', 'w'), sorted_embs)
    #embs_df = pd.DataFrame()
    #embs_df.index = ids_list
    #embs_df['uniprot_id'] = ids_list
    #embs_df['embedding'] = sorted_embs
    #print('Dataframe created:', embs_df.shape)

    #print(embs_df)

    #embs_df.to_parquet(output_dir+'/prottrans.parquet', compression='gzip')