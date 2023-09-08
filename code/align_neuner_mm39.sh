#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem=64GB # memory pool for all cores
#SBATCH -t 0-04:00 # time (D-HH:MM)
#SBATCH -o mm10_align_%A-%a.out # STDOUT
#SBATCH -e mm10_align_%A-%a.err # STDERR
#SBATCH --array=1-12

# Loading Singularity
# Singularity allows for containerized workflows, which facilitates long-term reproducibility
module load singularity

# Bringing in environment variables
# Note: GENOME, GTF, and TEGTF come from David Anderson's work
# Getting genome and annotation files
GENOMEDIR="/projects/kaczorowski-lab/david/genome"
GTF="/projects/kaczorowski-lab/david/genome/Mus_musculus.GRCm38.99.gtf"
TEGTF="/projects/kaczorowski-lab/david/genome/GRCm38_Ensembl_rmsk_TE.gtf"

# Getting FASTQ files
FASTQ_FILE_DIR="/projects/kaczorowski-lab/USERS/saulm/adbxd_te_coexp/data/GSE162526/"
FASTQ_FILE_KEY=$FASTQ_FILE_DIR"GSE162526_fastq_files.txt"
FASTQ_NAME=`cat $FASTQ_FILE_KEY | head -n $SLURM_ARRAY_TASK_ID | tail -n -1`
PREFIX=`echo $FASTQ_NAME | cut -d. -f1`
INFILE1=$FASTQ_FILE_DIR""$FASTQ_NAME

# Getting trim files
TRIMDIR="/fastscratch/saulm/goate_te_trimfiles/"
TRIMFILE1=$TRIMDIR""$PREFIX"_trimmed.fq"
ADAPTERFILE="/projects/kaczorowski-lab/USERS/saulm/adbxd_te_coexp/code/adapter/TruSeq3-SE.fa:2:30:10"

# Getting alignment files
BAMDIR="/fastscratch/saulm/goate_te_bamfiles/"$PREFIX"/"
BAMFILE=$BAMDIR"Aligned.sortedByCoord.out.bam"

# Getting out files
TEOUTFILE=$FASTQ_FILE_DIR""$PREFIX"_TE_counts"

echo "Started processing of $PREFIX on: `date`"

# Making directory for the alignment data
echo "Making $TRIMDIR on /fastscratch"
mkdir -p $TRIMDIR
echo "Making $BAMDIR on /fastscratch"
mkdir -p $BAMDIR

echo "Trimming FASTQ for $PREFIX on: `date`"

# Trimming RNAseq files
singularity exec docker://staphb/trimmomatic:0.38 \
    java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar \
    SE -threads 10 -phred33 \
    $INFILE1 $TRIMFILE1 \
    ILLUMINACLIP:$ADAPTERFILE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Aligning trimmed FASTQ for $PREFIX to GRCm38 with $GTF on: `date`"

# Running STAR for alignment to GRCm38/mm10
singularity exec docker://biocontainers/rna-star:v2.7.0adfsg-1-deb_cv1 \
	STAR --runThreadN 10 \
	--genomeDir $GENOMEDIR \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM \
	--quantTranscriptomeBAMcompression -1 \
	--quantTranscriptomeBan IndelSoftclipSingleend \
	--outFileNamePrefix $BAMDIR \
	--outSAMattributes NH HI AS nM \
	--readFilesIn $TRIMFILE1 \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100 \
	--sjdbGTFfile $GTF

echo "Running TETranscripts TEcount on $PREFIX for $TEGTF on: `date`"

# Running TETranscripts
singularity exec docker://mhammelllab/tetranscripts:2.2.3 \
    TEcount -b $BAMFILE --GTF $GTF \
    --TE $TEGTF \
    --sortByPos --stranded reverse \
    --project $TEOUTFILE

echo "Finished processing of $PREFIX on: `date`"
echo "Counts file saved as $TEOUTIFLE"
