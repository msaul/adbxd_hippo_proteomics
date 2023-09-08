#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 36 # number of cores
#SBATCH --mem=150GB # memory pool for all cores
#SBATCH -t 0-08:00 # time (D-HH:MM)
#SBATCH -o EncyclopeDIA_libprep_walnut.out # STDOUT
#SBATCH -e EncyclopeDIA_libprep_walnut.err # STDERR
#SBATCH --array=1

module load singularity
PROT_CONTAINER="/home/saulm/singularity/EncyclopeDIA.sif"

cd /fastscratch/saulm/

singularity exec $PROT_CONTAINER \
    java -Djava.awt.headless=true -Xmx148G -jar \
    /code/encyclopedia-1.4.10-executable.jar \
    -walnut -libexport -walnut -o DBA2J_chromatogram_library.elib \
    -i mzML -f ADBXD_qc_proteins.fasta -a false
