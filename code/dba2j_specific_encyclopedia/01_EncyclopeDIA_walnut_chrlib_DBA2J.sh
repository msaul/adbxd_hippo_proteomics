#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=16GB # memory pool for all cores
#SBATCH -t 0-08:00 # time (D-HH:MM)
#SBATCH -o chrlib_EncyclopeDIA_walnut_%A-%a.out # STDOUT
#SBATCH -e chrlib_EncyclopeDIA_walnut_%A-%a.err # STDERR
#SBATCH --array=1-186%40

module load singularity
PROT_CONTAINER="/home/saulm/singularity/EncyclopeDIA.sif"

MZML_FILE=`head -n $SLURM_ARRAY_TASK_ID ../adbxd_proteomics_chrlib_mzML.txt | tail -n 1 -`
cd /fastscratch/saulm/mzML/

singularity exec $PROT_CONTAINER \
    java -Djava.awt.headless=true -Xmx14G -jar \
    /code/encyclopedia-1.4.10-executable.jar \
    -walnut -i $MZML_FILE -f ADBXD_qc_proteins.fasta
