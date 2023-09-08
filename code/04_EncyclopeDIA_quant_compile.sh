#!/bin/bash
#SBATCH -p high_mem
#SBATCH -N 1 # number of nodes
#SBATCH -n 142 # number of cores
#SBATCH --mem=3020GB # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o EncyclopeDIA_quant_compile.out # STDOUT
#SBATCH -e EncyclopeDIA_quant_compile.err # STDERR
#SBATCH --array=1

# Loading singularity
module load singularity
PROTEOMICS_SIF="/home/saulm/singularity/EncyclopeDIA.sif"

# Setting temporary directory location
TMPDIR="/fastscratch/saulm/tmp"

# Getting to directory
cd /fastscratch/saulm/

singularity exec $PROTEOMICS_SIF \
    java -Xmx3018G -Djava.io.tmpdir=/fastscratch/saulm/tmp \
    -Djava.awt.headless=true -jar /code/encyclopedia-1.4.10-executable.jar \
    -libexport -o big_report.elib -l big_chromatogram_library.elib \ 
    -i quant -f mm_uniprot.fasta -a true
