from glob import glob
import gzip
import json
from multiprocessing import Pool
#import multiprocessing
from os import mkdir, path
import sys
from time import time
import numpy as np
from tqdm import tqdm
import ankh
import torch
import polars as pl

from sort_uniprot import read_uniprot_fasta
from util_base import chunks, run_command, split_list_by_maxtokens

class Embedder():
    def __init__(self, is_large, caches_dir) -> None:
        print('Loading model', file=sys.stderr)
        started = time()
        if is_large:
            self.model, self.tokenizer = ankh.load_large_model()
            self.model_name = 'ankh-large'
            self.emb_len = 1536
        else:
            self.model, self.tokenizer = ankh.load_base_model()
            self.model_name = 'ankh-base'
            self.emb_len = 768
        print('Model loaded in', time() - started, file=sys.stderr)
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
        self.cached_seqs = set()
        self.seq_to_path = {}
        print('Loading sequence embedding caches:', cache_paths, file=sys.stderr)
        for c_path in tqdm(caches_available):
            try:
                c_dict = json.load(gzip.open(c_path, 'rt'))
                for seq, emb in c_dict.items():
                    #self.cached_seqs[seq] = np.array(emb)
                    self.seq_to_path[seq] = c_path
                self.cached_seqs.update(c_dict.keys())

            except:
                continue
            
        print('Sequences in cache:', len(self.cached_seqs))

    def load_existing_embeddings(self, sequences):
        embeddings = [None for _ in sequences]
        seq_to_indexes = {}
        for i, s in enumerate(sequences):
            #sequences can repeat
            if s not in seq_to_indexes:
                seq_to_indexes[s] = []
            seq_to_indexes[s].append(i)

        loaded_cache = {}
        corrupted_caches = set()
        n_corrupted = 0
        n_not_really_in_cache = 0
        for seq in tqdm(sequences):
            if seq in self.cached_seqs:
                if seq not in loaded_cache:
                    c_path = self.seq_to_path[seq]
                    try:
                        for loaded_seq, emb in json.load(gzip.open(c_path, 'rt')).items():
                            loaded_cache[loaded_seq] = np.array(emb)
                    except:
                        corrupted_caches.add(c_path)
                        n_corrupted += 1
                emb = loaded_cache[seq]
                for i in seq_to_indexes[seq]:
                    embeddings[i] = emb
            else:
                n_not_really_in_cache += 1
        
        n_nans = sum([1 if a is None else 0 for a in embeddings])
        total_loaded = len(embeddings) - n_nans
        print(f'Loaded {total_loaded} embeddings from cache ({total_loaded/len(embeddings)*100:.2f}%)', 
            file=sys.stderr)
        if n_corrupted > 0:
            print(f'Seqs in corrupted caches: {n_corrupted}', file=sys.stderr)
            print('Corrupted caches:', file=sys.stderr)
            for c in corrupted_caches:
                print(c, file=sys.stderr)
        print(f'Not really in cache: {n_not_really_in_cache}', file=sys.stderr)
        for i, e in enumerate(embeddings):
            if e is None:
                embeddings[i] = np.nan
        return embeddings
    
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
        not_cached = {s: e.tolist() for s, e in zip(seqs, embeddings) if not s in self.cached_seqs}
        if len(not_cached) > 0:
            cache_i = 1
            cache_name = self.base_cache_name.replace('_CACHEN_', '_'+str(cache_i))
            while path.exists(cache_name):
                cache_i += 1
                cache_name = self.base_cache_name.replace('_CACHEN_', '_'+str(cache_i))
            output_stream = gzip.open(cache_name, 'wt')
            json.dump(not_cached, output_stream)
            for seq, emb in not_cached.items():
                self.cached_seqs.add(seq)
                self.seq_to_path[seq] = cache_name
    
    def calc_embeddings_batched(self, original_seqs, use_cache=True):
        started = time()
        print('Total sequences:', len(original_seqs))
        print('Cache size:', len(self.cached_seqs))
        n_cached = len([1 for s in original_seqs if s in self.cached_seqs])
        print('Sequences to embed already in cache:', n_cached)
        seqs = [s for s in original_seqs if not s in self.cached_seqs]
        print('Not cached:', len(seqs))
        seqs = list(set(seqs))
        seqs.sort(key=lambda s: len(s))
        print('Not cached (non redundant):', len(seqs))
        if len(seqs) > 0:
            seq_chunks = split_list_by_maxtokens(seqs, 5000)
            total_seqs = sum([len(s) for s in seq_chunks])
            iterator = tqdm(total=total_seqs)
            for not_cached_seqs in seq_chunks:
                calculated_embs = self.aminos_to_embeddings(not_cached_seqs)

                if use_cache:
                    self.cache_results(not_cached_seqs, calculated_embs)
                if iterator:
                    iterator.update(len(not_cached_seqs))
                for s, e in zip(not_cached_seqs, calculated_embs):
                    self.cached_seqs.add(s)
            if iterator:
                iterator.close()
        embeddings_list = self.load_existing_embeddings(original_seqs)
        #embeddings_list = [np.array(self.cached_seqs[s]) for s in original_seqs]
        #embeddings_list = [np.nan for _ in original_seqs]
        #seq_to_index = {s: i for i, s in enumerate(original_seqs)}
        #for c_path in self.cached_seqs:
        #    embeddings_list[seq_to_index[s]] = self.cached_seqs[s]

        elapsed = time() - started
        return embeddings_list, elapsed

def embed_sequences(is_large, fasta_path, caches_path):
    fasta = read_uniprot_fasta(fasta_path)
    seq_names = [h.replace('|', ' ').split(' ')[0] for h, x in fasta]
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
    
    output_suffix = ""
    if len(sys.argv) > 4:
        output_suffix = sys.argv[4]

    nonetype = type(None)

    for is_large in [False]:
        size_str = 'large' if is_large else 'base'
        if output_suffix == '':
            output_pq = 'emb.ankh_'+size_str+'.parquet'
        else:
            output_pq = 'emb.ankh_'+size_str+'.'+output_suffix+'.parquet'
        seq_names, embeddings = embed_sequences(is_large, fasta_path, caches_path)

        emb_shape = embeddings[0].shape

        print('Embeddings shapes:')
        n_shapes = {'float': 0}
        for e in embeddings:
            try:
                n_shapes[e.shape] = n_shapes.get(e.shape, 0) + 1
            except AttributeError as err:
                n_shapes['float'] += 1
        for shape, count in n_shapes.items():
            print(f'{shape}: {count}')

        print('Sorting embeddings')
        uniprot_to_index = {id: i for i, id in enumerate(all_ids)}
        all_embeddings = [None for _ in all_ids]
        for id, emb in zip(seq_names, embeddings):
            index = uniprot_to_index[id]
            all_embeddings[index] = emb
        
        print("Filling missing proteins with np.nan")
        for i in range(len(all_embeddings)):
            if type(all_embeddings[i]) == nonetype:
                all_embeddings[i] = np.empty(emb_shape)
                all_embeddings[i].fill(np.nan)
        
        print('Embeddings shapes after filling missing proteins with np.nan:')
        n_shapes = {'float': 0}
        for e in all_embeddings:
            try:
                n_shapes[e.shape] = n_shapes.get(e.shape, 0) + 1
            except AttributeError as err:
                n_shapes['float'] += 1
        for shape, count in n_shapes.items():
            print(f'{shape}: {count}')

        print('Converting to numpy matrix')
        all_embeddings = np.asarray(all_embeddings)
        print('Creating polars df')
        df = pl.DataFrame({
            'id': all_ids,
            'emb': all_embeddings
        })
        print('Saving to file', output_pq)
        df.write_parquet(output_pq)
        #np.save(output_np, all_embeddings, allow_pickle=False)
        #print('compressing', output_np)
        #run_command(['gzip', output_np])#