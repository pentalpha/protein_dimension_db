#!/bin/bash
mkdir -p singularity/sif
sudo singularity build singularity/sif/ankh.sif singularity/def/ankh.def > singularity/sif/ankh.stdout 2>&1
sudo singularity build singularity/sif/env2.sif singularity/def/env2.def > singularity/sif/env2.stdout 2>&1
sudo singularity build singularity/sif/basic_env.sif singularity/def/basic_env.def > singularity/sif/basic_env.stdout 2>&1
