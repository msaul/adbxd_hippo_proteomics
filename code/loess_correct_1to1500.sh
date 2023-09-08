#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=16GB # memory pool for all cores
#SBATCH -t 1-00:00 # time (D-HH:MM)
#SBATCH -o ploess_%A-%a.out # STDOUT
#SBATCH -e ploess_%A-%a.err # STDERR
#SBATCH --array=1-1500%100

# Using R v4.2.1
R_CONTAINER="/home/saulm/singularity/rocker_rstudio_4.2.1.sif"
module load singularity

singularity exec $R_CONTAINER Rscript /projects/kaczorowski-lab/USERS/saulm/adbxd_proteomics/code/loess_correct_1to1500.R
