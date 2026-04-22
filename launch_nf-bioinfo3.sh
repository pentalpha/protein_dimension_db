# 1. Create the cache and temp directories inside your current working directory
mkdir -p ${PWD}/.cache/tmp
mkdir -p ${PWD}/.cache/singularity

# 2. Tell Singularity where to unpack images and store cache
export SINGULARITY_TMPDIR="${PWD}/.cache/tmp"
export SINGULARITY_CACHEDIR="${PWD}/.cache/singularity"

# 3. Tell Nextflow and standard Linux tools where to write temp files
export TMPDIR="${PWD}/.cache/tmp"
export NXF_TEMP="${PWD}/.cache/tmp"

nextflow -C nextflow-bioinfo3.config run main-uniprot2026.nf -resume --mode full --release_dir /data/home/pitagoras/data/dimension_db/release_2 --old_release_paths_str /data/home/pitagoras/data/dimension_db/cafa6,/data/home/pitagoras/data/dimension_db/release_1