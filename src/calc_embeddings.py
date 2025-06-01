import sys
from os import path

from util_base import run_command

prj_dir = path.dirname(path.dirname(__file__))

MAX_PROTEIN_LENGTH = 1800

def calc_embedding_process(params: str):
    short_name = params['short_name']
    model_name = params['model_name']
    release_dir = params['release_dir']
    fasta_path = release_dir + '/uniprot_sorted.fasta.gz'
    embedding_path = release_dir + '/emb.'+short_name+'.npy.gz'

    cache_dir = path.dirname(release_dir)+'/embs_cache_tmp'
    if not path.exists(cache_dir):
        run_command(['mkdir -p', cache_dir])
    if path.exists(embedding_path):
        run_command(['rm', embedding_path])
    cmd = ['python', prj_dir+'/src/embedding_calculation.py', model_name, fasta_path, cache_dir, embedding_path, str(MAX_PROTEIN_LENGTH)]
    run_command(cmd)

if __name__ == "__main__":
    release_dir = sys.argv[1]

    models_df_path = prj_dir + '/others/model_sizes.csv'
    model_names = []
    models_input_df = open(models_df_path, 'r')
    models_input_df.readline()
    for rawline in models_input_df:
        cells = rawline.rstrip("\n").split(',')
        if len(cells[-1]) > 1:
            model_names.append((cells[0], cells[-1]))
    models_input_df.close()

    print('Models listed to use:')
    params_list = []
    for a, b in model_names:
        print(a, b)
        params_list.append({
            'short_name': b.strip('"'),
            'model_name': a.strip('"'),
            'release_dir': release_dir
        })
    
    for p in params_list:
        print(p)
        calc_embedding_process(p)