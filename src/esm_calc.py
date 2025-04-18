from os import mkdir, path
from glob import glob
from pickle import dump
import sys
import torch
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import polars as pl

from util_base import proj_dir, run_command
from fasta import fasta_equal_split, fasta_equal_split_by_len, remove_from_fasta

#esm_script_path = proj_dir + '/libs/esm/scripts/extract.py'
esm_script_path = 'esm/scripts/extract.py'

def run_esm_extract(args: dict):
    fasta_input = args['fasta_input']
    esm_model_name = args['esm_model_name']
    esm_output_dir = args['esm_output_dir']
    batch_size = args['batch_size']
    #tmp_output_dir = fasta_input.rstrip('.fasta')+'_output'
    #run_command(['mkdir', tmp_output_dir])
    esm_cmd = ['python', esm_script_path, esm_model_name, 
        fasta_input, esm_output_dir, '--include mean --toks_per_batch', str(batch_size)]
    run_command(esm_cmd)

def download_esm_model(model_name):
    fake_fasta_content = '>seqnameA\nMATPGASSARDEFVYMAKLAEQAERYEEMVTHPIRLGLALNFSVFYYEI\n'
    fake_fasta_path = 'fakefasta.fasta'
    open(fake_fasta_path, 'w').write(fake_fasta_content)

    download_cmd = ['python', esm_script_path, model_name, 
        fake_fasta_path, 'output_tmp', '--include mean --toks_per_batch 10000']
    run_command(download_cmd)

    run_command(['rm -Rf', fake_fasta_path, 'output_tmp'])

class ESM_Embedder():

    models_names = {
        'esm2_t36_3B_UR50D',
        'esm2_t33_650M_UR50D',
        'esm2_t30_150M_UR50D',
        'esm2_t12_35M_UR50D',
        'esm2_t6_8M_UR50D'
    }

    models_dims = {
        'esm2_t36_3B_UR50D': 2560,
        'esm2_t33_650M_UR50D': 1280,
        'esm2_t30_150M_UR50D': 640,
        'esm2_t12_35M_UR50D': 480,
        'esm2_t6_8M_UR50D': 320
    }

    models_n_layers = {
        'esm2_t36_3B_UR50D': 36,
        'esm2_t33_650M_UR50D': 33,
        'esm2_t30_150M_UR50D': 30,
        'esm2_t12_35M_UR50D': 12,
        'esm2_t6_8M_UR50D': 6
    }

    models_params = {
        'esm2_t36_3B_UR50D': {'batch_max': 10000, 'processes': 1},
        'esm2_t33_650M_UR50D': {'batch_max': 10000, 'processes': 1},
        'esm2_t30_150M_UR50D': {'batch_max': 100000, 'processes': 1},
        'esm2_t12_35M_UR50D': {'batch_max': 200000, 'processes': 1},
        'esm2_t6_8M_UR50D': {'batch_max': 600000, 'processes': 1}
    }

    def __init__(self, cache_dir, model_full_name) -> None:
        self.model_name = model_full_name.split('/')[-1]
        assert self.model_name in ESM_Embedder.models_names, (model_full_name 
            + ' / '+ self.model_name + ' not in' + str(ESM_Embedder.models_names))
        self.processes = ESM_Embedder.models_params[self.model_name]['processes']
        self.batch_max = ESM_Embedder.models_params[self.model_name]['batch_max']
        self.model_dim = ESM_Embedder.models_dims[self.model_name]
        self.n_layers = ESM_Embedder.models_n_layers[self.model_name]
        self.shape = (self.model_dim,)

        print('Loading fair-esm embedding utility')
        self.embs_dir = cache_dir + '/' + model_full_name.replace('/', '_') + '_cache'
        if not path.exists(self.embs_dir):
            mkdir(self.embs_dir)
        self.calculated = set()
        self.find_calculated()
    
    def find_calculated(self):
        emb_files = glob(self.embs_dir+'/*.pt')
        uniprotids = [path.basename(f.rstrip('.pt')) for f in emb_files]
        self.calculated.update(uniprotids)
        print('ESM', self.model_name, ':', len(self.calculated), 'proteins calculated')

    def load_embedding(self, uniprotid: str):
        if uniprotid in self.calculated:
            emb_path = path.join(self.embs_dir, uniprotid+'.pt')
            x = torch.load(emb_path)
            emb = x['mean_representations'][self.n_layers].tolist()
            return emb
        else:
            print(uniprotid, 'not found at dictionary with', 
                len(self.calculated), 'keys', file=sys.stderr)
            empty_vec = np.empty(self.shape)
            empty_vec.fill(np.nan)
            return empty_vec
    
    def get_embeddings(self,uniprotids: list):
        print('Loading embeddings for', len(uniprotids), 'proteins')
        embs_list = [self.load_embedding(ID) for ID in tqdm(uniprotids)]
        return embs_list
    
    def export_embeddings(self, uniprotids: list, output_path: str):
        embs_list = self.get_embeddings(uniprotids)

        all_embeddings = np.asarray(embs_list)
        print('Saving to file', output_path)
        df = pl.DataFrame({
            'id': uniprotids,
            'emb': all_embeddings
        })
        print('Saving to file', output_path)
        df.write_parquet(output_path)
        #np.save(output_path, all_embeddings, allow_pickle=False)
        #print('compressing', output_path)
        #run_command(['gzip', output_path])

        return output_path

    def calc_embeddings(self, input_fasta: str):
        calculated_proteins = self.calculated
        print('Q6UY62 is in calculated_proteins', 'Q6UY62' in calculated_proteins)
        if len(calculated_proteins) > 0:
            print('Some embeddings have already been calculated')
            print('Removing them from the input fasta')
            to_process_fasta = self.embs_dir+'/to_calc.fasta'
            kept = remove_from_fasta(input_fasta, calculated_proteins, to_process_fasta)
            print('Q6UY62 is in kept', 'Q6UY62' in kept)
            if len(kept) == 0:
                run_command(['rm', to_process_fasta])
                print('All proteins calculated')
                return
            input_fasta = to_process_fasta
        
        fasta_parts = fasta_equal_split_by_len(input_fasta, self.processes)
        for f in fasta_parts:
            print(f)
        
        download_esm_model(self.model_name)
        
        '''if path.exists(esm_output_dir):
            run_command(['rm -rf', esm_output_dir])'''
        esm_run_params = [{'fasta_input': fasta_part, 'esm_model_name': self.model_name,
                        'esm_output_dir': self.embs_dir, 'batch_size': self.batch_max}
            for fasta_part in fasta_parts]
        
        with Pool(self.processes) as pool:
            print("Parallel processing with", self.processes, 'processes')
            pool.map(run_esm_extract, esm_run_params)

        pt_files = glob(self.embs_dir+'/*.pt')
        print(len(pt_files), 'proteins with embedding')
        uniprotids = [path.basename(f.rstrip('.pt')) for f in pt_files]
        self.calculated.update(uniprotids)

        for f in fasta_parts:
            run_command(['rm', f])

if __name__ == "__main__":
    fasta_input_path = sys.argv[1]
    cache_dir = sys.argv[2]
    all_uniprot_ids_path = sys.argv[3]
    all_ids = open(all_uniprot_ids_path, 'r').read().split('\n')
    models_meta_info_csv_path = 'others/model_sizes.csv'
    facebook_models = []

    for rawline in open(models_meta_info_csv_path, 'r'):
        cells = rawline.rstrip('\n').split(',')
        model_full_name = cells[0].strip('"')
        short_name = cells[-1].strip('"')
        if 'facebook' in model_full_name:
            print(cells)
            facebook_models.append((model_full_name, short_name))

    for model_full_name, short_name in facebook_models:
        embedder = ESM_Embedder(cache_dir, model_full_name)
        embedder.calc_embeddings(fasta_input_path)
        output_path = 'emb.'+short_name+'.parquet'
        if not 'esm2_t' in output_path:
            output_path = output_path.replace('esm2_', 'esm2_t')
        embedder.export_embeddings(all_ids, output_path)