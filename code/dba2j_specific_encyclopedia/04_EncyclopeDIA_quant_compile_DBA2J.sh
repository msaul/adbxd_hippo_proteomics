#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 36 # number of cores
#SBATCH --mem=500GB # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o EncyclopeDIA_quant_compile.out # STDOUT
#SBATCH -e EncyclopeDIA_quant_compile.err # STDERR
#SBATCH --array=1

# Loading singularity
module load singularity
PROTEOMICS_SIF="/home/saulm/singularity/EncyclopeDIA.sif"

# Setting temporary directory location
mkdir -p /fastscratch/saulm/tmp
TMPDIR="/fastscratch/saulm/tmp"

# Getting to directory
cd /fastscratch/saulm/

singularity exec $PROTEOMICS_SIF \
    java -Xmx498G -Djava.io.tmpdir=/fastscratch/saulm/tmp -Djava.awt.headless=true -jar /code/encyclopedia-1.4.10-executable.jar -libexport -o ADBXD_out_report.elib -l DBA2J_chromatogram_library.elib -i /fastscratch/saulm/quant/ -f ADBXD_qc_proteins.fasta -a true
