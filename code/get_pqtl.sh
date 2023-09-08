#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=16GB # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o pqtl.out # STDOUT
#SBATCH -e pqtl.err # STDERR

# Using R v4.2.1
R_CONTAINER="/home/saulm/singularity/rocker_rstudio_4.2.1.sif"
module load singularity

# Changing to Kaczorowski Lab folder
cd /projects/kaczorowski-lab/USERS/saulm/adbxd_proteomics/

# Running pQTL analysis
singularity exec $R_CONTAINER Rscript -e 'workflowr::wflow_publish("./analysis/104_pqtl_mapping.Rmd", "Initial commit of pQTL analysis.")'
