#!/bin/bash
#SBATCH --mail-user=alves.pitagoras@gmail.com
#SBATCH --partition=intel-128
#SBATCH --mail-type=ALL
#SBATCH --mem 5G
#SBATCH --cpus-per-task=4        
#SBATCH --time=01:00:00

CONFIG=$1

eval "$(conda shell.bash hook)" # inclua esta linha antes executar o conda
conda activate /home/pdaasobrinho/.conda/envs/nextflow

nextflow -C ${CONFIG} run main.nf -resume --mode test --release_dir /home/pdaasobrinho/data/pddb_test