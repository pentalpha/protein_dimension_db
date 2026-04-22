#!/bin/bash
mkdir -p singularity/sif
#sudo singularity build singularity/sif/ankh_24.sif singularity/def/ankh_24.def > singularity/sif/ankh_24.stdout 2>&1
sudo singularity build singularity/sif/ankh.sif singularity/def/ankh.def > singularity/sif/ankh.stdout 2>&1
sudo singularity build singularity/sif/env2.sif singularity/def/env2.def > singularity/sif/env2.stdout 2>&1
sudo singularity build singularity/sif/basic_env.sif singularity/def/basic_env.def > singularity/sif/basic_env.stdout 2>&1
sudo singularity build singularity/sif/obonet_scipy.sif singularity/def/obonet_scipy.def > singularity/sif/obonet_scipy.stdout 2>&1
sudo singularity build singularity/sif/python_legacy.sif singularity/def/python_legacy.def > singularity/sif/python_legacy.stdout 2>&1