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
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

from sort_uniprot import read_uniprot_fasta
from util_base import chunks, run_command, split_list_by_maxtokens

class ESMCEmbedder():
    def __init__(self, model_name, caches_dir) -> None:
        client = ESMC.from_pretrained(model_name).to("cpu") # or "cpu"
        if model_name == 'esmc_600m':
            self.emb_len = 1152
        elif model_name == 'esmc_300m':
            self.emb_len = 960
        self.model_name = model_name
        self.emb_shape = (self.emb_len,)
        
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
        proteins = [ESMProtein(sequence=s) for s in protein_sequences]
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
    
    def calc_embeddings_batched(self, original_seqs, use_cache=True):
        started = time()
        print('Total sequences:', len(original_seqs))
        print('Cache size:', len(self.cached_seqs.keys()), len(self.cached_seqs.values()))
        n_cached = len([1 for s in original_seqs if s in self.cached_seqs])
        print('Sequences to embed already in cache:', n_cached)
        seqs = [s for s in original_seqs if not s in self.cached_seqs]
        print('Not cached:', len(seqs))
        seqs = list(set(seqs))
        seqs.sort(key=lambda s: len(s))
        print('Not cached (non redundant):', len(seqs))
        if len(seqs) > 0:
            seq_chunks = split_list_by_maxtokens(seqs, 100000)
            total_iter = len(seq_chunks)
            iterator = tqdm(total=total_iter)
            for not_cached_seqs in seq_chunks:
                calculated_embs = self.aminos_to_embeddings(not_cached_seqs)

                if use_cache:
                    self.cache_results(not_cached_seqs, calculated_embs)
                if iterator:
                    iterator.update(1)
                for s, e in zip(not_cached_seqs, calculated_embs):
                    self.cached_seqs[s] = e.tolist()
            if iterator:
                iterator.close()
        embeddings_list = [np.array(self.cached_seqs[s]) for s in original_seqs]
        elapsed = time() - started
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