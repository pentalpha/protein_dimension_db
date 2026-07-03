#!/bin/bash
#SBATCH --mail-user=alves.pitagoras@gmail.com
#SBATCH --partition=intel-128
#SBATCH --mail-type=ALL
#SBATCH --mem 5G
#SBATCH --cpus-per-task=4        
#SBATCH --time=01:00:00

eval "$(conda shell.bash hook)" # inclua esta linha antes executar o conda
conda activate /home/pdaasobrinho/.conda/envs/pytorch
module load singularity/3.7.1

nextflow -C nextflow-slurm.config run infer_embeddings.nf -resume --release_dir /home/pdaasobrinho/data/dimension_db/release_2