#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=16GB # memory pool for all cores
#SBATCH -t 0-08:00 # time (D-HH:MM)
#SBATCH -o EncyclopeDIA_quant_%A-%a.out # STDOUT
#SBATCH -e EncyclopeDIA_quant_%A-%a.errc # STDERR
#SBATCH --array=1-441%40

module load singularity
PROT_CONTAINER="/home/saulm/singularity/EncyclopeDIA.sif"

MZML_FILE=`head -n $SLURM_ARRAY_TASK_ID adbxd_hippo_proteomics_mzML.txt | tail -n 1 -`
cd /fastscratch/saulm/

singularity exec $PROT_CONTAINER \
    java -Djava.awt.headless=true -Xmx14G -jar \
    /code/encyclopedia-1.4.10-executable.jar \
    -i "./quant/"$MZML_FILE \
    -l big_chromatogram_library.elib \
    -f mm_uniprot.fasta
