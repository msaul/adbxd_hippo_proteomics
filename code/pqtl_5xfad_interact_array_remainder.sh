#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=16GB # memory pool for all cores
#SBATCH -t 0-02:00 # time (D-HH:MM)
#SBATCH -o pqtl_%A-%a.out # STDOUT
#SBATCH -e pqtl_%A-%a.err # STDERR
#SBATCH --array=85,88,90,96,97,99,100

# Using R v4.2.1
R_CONTAINER="/home/saulm/singularity/rocker_rstudio_4.2.1.sif"
module load singularity

singularity exec $R_CONTAINER Rscript /projects/kaczorowski-lab/USERS/saulm/adbxd_proteomics/code/adbxd_pqtl_5xfad_interact_analysis.R
