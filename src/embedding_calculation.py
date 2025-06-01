from glob import glob
import gzip
import json
from multiprocessing import Pool
from os import mkdir, path
import sys
from time import time
import numpy as np
import torch
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from transformers import EsmTokenizer, EsmModel

from sort_uniprot import read_uniprot_fasta
from util_base import chunks, run_command

class Embedder():
    def __init__(self, model_name, caches_dir, max_prot_length=2000) -> None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Using device: {device}", file=sys.stderr)
        self.model_library = "ankh"
        self.cache_dir = caches_dir
        self.model_name = model_name
        self.model_name_simple = self.model_name.replace('/', '_')
        if 'esm' in model_name:
            self.model_library = 'esm'
        
        if self.model_library == 'ankh':
            tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
            model = AutoModelForSeq2SeqLM.from_pretrained(model_name, trust_remote_code=True).to(device)
        elif self.model_library == 'esm':
            tokenizer = EsmTokenizer.from_pretrained(model_name, trust_remote_code=True)
            model = EsmModel.from_pretrained(model_name, trust_remote_code=True).to(device)
        print('Loaded models', file=sys.stderr)
        decoder_input_ids = tokenizer("<s>", return_tensors="pt").input_ids.to(device)

        self.device = device
        self.tokenizer = tokenizer
        self.model = model
        self.decoder_input_ids = decoder_input_ids

        self.base_cache_name = self.cache_dir+'/'+self.model_name_simple+'_CACHEN_.json.gz'
        
        self.max_prot_length = max_prot_length
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
    
    def amino_to_embedding(self, seq, use_cache=True):
        if seq in self.cached_seqs:
            return self.cached_seqs[seq]
        else:
            inputs = self.tokenizer(seq, return_tensors="pt", add_special_tokens=True).to(self.device)
            #print('inputs')
            if self.model_library == 'esm':
                outputs = self.model(input_ids=inputs["input_ids"], output_hidden_states=True)
                emb = outputs.last_hidden_state
            else:
                outputs = self.model(input_ids=inputs["input_ids"], decoder_input_ids=self.decoder_input_ids)
                emb = outputs.encoder_last_hidden_state
            #print('outputs')
            #print('Processed', seq)
            #emb = outputs.encoder_last_hidden_state
            protein_emb = emb[0].mean(dim=0)
            sequence_embedding = protein_emb.detach().cpu().numpy().squeeze()
            #print(sequence_embedding)
            return sequence_embedding
    
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
    
    def calc_embeddings_batched(self, seqs, seq_chunk_length, use_cache=True, detailed=False):
        seq_chunks = chunks(seqs, seq_chunk_length)
        iterator = tqdm(seqs)
        embeddings_list = []
        last_elapsed = time()
        started = None
        emb_shape = None
        for seq_chunk in seq_chunks:
            new_embs = []
            for s in seq_chunk:
                if detailed:
                    print('input seq has', len(s), 'aminoacids', file=sys.stderr)
                if len(s) > self.max_prot_length:
                    emb = np.empty(emb_shape)
                    emb = emb.fill(np.nan)
                else:
                    emb = self.amino_to_embedding(s, use_cache=use_cache)
                    emb_shape = emb.shape
                new_embs.append(emb)
                iterator.update(1)
                if detailed:
                    print('Done in', time() - last_elapsed, file=sys.stderr)
                if not started:
                    started = time()
                last_elapsed = time()
            if use_cache:
                self.cache_results(seq_chunk, new_embs)
            embeddings_list += new_embs
        #embeddings_list = [self.amino_to_embedding(s) for s in iterator]
        elapsed = time() - started
        return np.asarray(embeddings_list), elapsed

def run_embedding_time_benchmark(model_name, fasta_path, caches_path):
    #model_names = ['ElnaggarLab/ankh-base', 'ElnaggarLab/ankh-large', 'ElnaggarLab/ankh2-ext2']
    '''seqs = ['HSDAIFTDNYSRFRKQMAVKKYLNSVLT',
        'VGCEECPMHCKGKHAVPTCDDGVCNCNV',
        'HSDAVFTDNYTRLRKQMAVKKYLNSILN',
        'ALFSILRGLKKLGNMGQAFVNCKIYKKC',
        'GIMDSVKGLAKNLAGKLLDSLKCKITGC',
        'SYSMEHFRWGKPVGRKRRPVKVYPNGVQEETSEGFPLEF',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL',
        'SYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPELSYSMEHFRWGKPMGRKRRPIKVYPNSFEDESVENMGPEL']'''

    fasta = read_uniprot_fasta(fasta_path)
    seqs = [s for h, s in fasta]

    seq_lengths = [50, 100, 200, 400, 500, 
                   600, 800, 900, 1000,
                   1100, 1200, 1300, 1400,
                   1500, 1600, 1700, 1800]

    n_proteins = [(l, len([s for s in seqs if len(s) >= l])) for l in seq_lengths]
    with open(fasta_path+'.n_proteins.txt', 'w') as nprots_out:
        for l, n in n_proteins:
            print(n, '>=', l, file=nprots_out)
    
    seq_len_i = 0
    seq_examples = []
    seq_i = 0
    while seq_len_i < len(seq_lengths) and seq_i < len(seqs):
        if len(seqs[seq_i]) > seq_lengths[seq_len_i]:
            seq_examples.append(seqs[seq_i])
            seq_len_i += 1
        seq_i += 1
    print('benchmarking on', len(seq_examples), 'sequences', file=sys.stderr)
    print(model_name, file=sys.stderr)
    model = Embedder(model_name, caches_path)
    total_time = 0.0
    final_shape = None
    #for _ in range(4):
    embeddings, elapsed = model.calc_embeddings_batched(seq_examples, 100, use_cache=False, detailed=True)
    total_time += elapsed
    final_shape = embeddings.shape
        
    print(model_name, total_time, final_shape)

def embed_sequences(model_name, fasta_path, caches_path, max_prot_length):
    fasta = read_uniprot_fasta(fasta_path)
    seqs = [s for h, s in fasta]

    print('embedding on', len(seqs), 'sequences with', file=sys.stderr)
    print(model_name, file=sys.stderr)
    model = Embedder(model_name, caches_path, max_prot_length=max_prot_length)
    
    embeddings, total_time = model.calc_embeddings_batched(seqs, 10000, use_cache=True, detailed=False)
    final_shape = embeddings.shape
        
    print(model_name, total_time, final_shape)
    return embeddings

if __name__ == "__main__":
    model_name = sys.argv[1]
    fasta_path = sys.argv[2]
    caches_path = sys.argv[3]
    if len(sys.argv) >= 6:
        output_npzgz = sys.argv[4]
        max_prot_length = int(sys.argv[5])
        embeddings = embed_sequences(model_name, fasta_path, caches_path, max_prot_length)
        print('saving', model_name)
        np.save(output_npzgz.rstrip(".gz"), embeddings, allow_pickle=False)
        print('compressing', output_npzgz)
        run_command(['gzip', output_npzgz.rstrip(".gz")])
    else:
        run_embedding_time_benchmark(model_name, fasta_path, caches_path)