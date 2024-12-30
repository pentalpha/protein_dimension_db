from glob import glob
import gzip
import json
from multiprocessing import Pool
import multiprocessing
from os import mkdir, path
import sys
from time import time
import numpy as np
from tqdm import tqdm
import ankh
import torch

from sort_uniprot import read_uniprot_fasta
from util_base import chunks, run_command

class Embedder():
    def __init__(self, is_large, caches_dir) -> None:
        if is_large:
            self.model, self.tokenizer = ankh.load_large_model()
            self.model_name = 'ankh-large'
            self.emb_len = 1536
        else:
            self.model, self.tokenizer = ankh.load_base_model()
            self.model_name = 'ankh-base'
            self.emb_len = 768
        self.emb_shape = (self.emb_len,)
        self.model.eval()
        self.cache_dir  = caches_dir + '/'+self.model_name

        self.base_cache_name = self.cache_dir+'/'+self.model_name+'_CACHEN_.json.gz'
        
        self.load_cache()

    def load_cache(self):
        if not path.exists(self.cache_dir):
            mkdir(self.cache_dir)
        cache_paths = self.base_cache_name.replace('_CACHEN_', '_*')
        caches_available = glob(cache_paths)
        self.cached_seqs = {}
        print('Loading sequence embedding caches:', cache_paths, file=sys.stderr)
        for c_path in tqdm(caches_available):
            c_dict = json.load(gzip.open(c_path, 'rt'))
            for seq, emb in c_dict.items():
                self.cached_seqs[seq] = np.array(emb)
        print('Sequences in cache:', len(self.cached_seqs))
    
    def aminos_to_embeddings(self, protein_sequences):
        #torch.set_num_threads(int(multiprocessing.cpu_count()*0.7))
        protein_sequences = [list(seq) for seq in protein_sequences]
        outputs = self.tokenizer.batch_encode_plus(protein_sequences, 
            add_special_tokens=True, 
            padding=True, 
            is_split_into_words=True, 
            return_tensors="pt")
        
        with torch.no_grad():
            embeddings = self.model(input_ids=outputs['input_ids'], attention_mask=outputs['attention_mask'])
            n_seqs = embeddings.last_hidden_state.shape[0]
            embeddings_np = [embeddings.last_hidden_state[i].mean(dim=0).numpy() for i in range(n_seqs)]

            return embeddings_np
    
    def cache_results(self, seqs, embeddings):
        not_cached = {s: e.tolist() for s, e in zip(seqs, embeddings) if not s in self.cached_seqs.keys()}
        if len(not_cached) > 0:
            cache_i = 1
            cache_name = self.base_cache_name.replace('_CACHEN_', '_'+str(cache_i))
            while path.exists(cache_name):
                cache_i += 1
                cache_name = self.base_cache_name.replace('_CACHEN_', '_'+str(cache_i))
            output_stream = gzip.open(cache_name, 'wt')
            json.dump(not_cached, output_stream)
    
    def calc_embeddings_batched(self, seqs, use_cache=True):
        parts = 500
        part_size = int(len(seqs)/parts)
        seq_chunks = chunks(seqs, part_size)
        iterator = tqdm(total=parts)
        embeddings_list = []
        started = time()
        for seq_chunk in seq_chunks:
            not_cached = [(i, seq_chunk[i]) 
                          for i in range(len(seq_chunk)) 
                          if not seq_chunk[i] in self.cached_seqs]
            not_cached_indexes = [i for i, s in not_cached]
            not_cached_seqs = [s for i, s in not_cached]

            new_embs = [None if not s in self.cached_seqs else self.cached_seqs[s] 
                        for s in range(len(seq_chunk))]
            
            if len(not_cached_seqs) > 0:
                print('Chunk is', (1.0-(len(not_cached_seqs)/len(seq_chunk)))*100, '% cached')
                calculated_embs = self.aminos_to_embeddings(not_cached_seqs)
                for i, emb in zip(not_cached_indexes, calculated_embs):
                    new_embs[i] = emb
            else:
                print('Cached chunk')
            
            if use_cache:
                self.cache_results(seq_chunk, new_embs)
            iterator.update(1)
            embeddings_list += new_embs
        #embeddings_list = [self.amino_to_embedding(s) for s in iterator]
        elapsed = time() - started
        iterator.close()
        return embeddings_list, elapsed

def embed_sequences(is_large, fasta_path, caches_path):
    fasta = read_uniprot_fasta(fasta_path)
    seq_names = [h for h, x in fasta]
    seqs = [s for h, s in fasta]

    print('embedding on', len(seqs), 'sequences with', file=sys.stderr)
    model = Embedder(is_large, caches_path)
    print(model.model_name, file=sys.stderr)
    
    embeddings, total_time = model.calc_embeddings_batched(seqs, use_cache=True)
    
    print(model.model_name, total_time)
    return seq_names, embeddings

if __name__ == "__main__":
    fasta_path = sys.argv[1]
    caches_path = sys.argv[2]
    all_ids_path = sys.argv[3]
    all_ids = open(all_ids_path, 'r').read().split('\n')

    nonetype = type(None)

    for is_large in [False, True]:
        if is_large:
            output_np = 'emb.ankh_large.npy'
        else:
            output_np = 'emb.ankh_base.npy'
        
        seq_names, embeddings = embed_sequences(is_large, fasta_path, caches_path)

        print("Filling missing proteins with np.nan")
        emb_shape = embeddings[0].shape

        uniprot_to_index = {id: i for i, id in enumerate(all_ids)}
        all_embeddings = [None for _ in all_ids]
        for id, emb in zip(seq_names, embeddings):
            index = uniprot_to_index[id]
            all_embeddings[index] = emb
        
        for i in range(len(all_embeddings)):
            if type(all_embeddings[i]) == nonetype:
                all_embeddings[i] = np.empty(emb_shape)
                all_embeddings[i].fill(np.nan)
        
        print('Converting to numpy matrix')
        all_embeddings = np.asarray(all_embeddings)
        print('Saving to file', output_np)
        np.save(output_np, all_embeddings, allow_pickle=False)
        print('compressing', output_np)
        run_command(['gzip', output_np])