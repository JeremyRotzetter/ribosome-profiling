#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=03:00:00
#SBATCH --job-name=mapping
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/03_output_mapping_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/03_error_mapping_%j.e
#SBATCH --partition=pall

# This script will use bowtie to map the trimmed reads from the previous step to the genome, sort them. These sorted BAM files
# would be used for QC purpose and differential expression analysis.

# How to use: User needs to input path to the directory where the input data is (SOURCE), the path to the directory where the indices for alignment are
# and the path to the directory where the output data should go (DESTINATION) when running the script with sbatch
# Example: sbatch 03_mapping.sh </path/to/source/> </path/to/indices/> </path/to/destination/>

SOURCE=$1 # the trimmed reads
INDICES=$2
DESTINATION=$3

# Load all the modules required for the alignment
module load UHTS/Aligner/bowtie/1.2.0
module load UHTS/Analysis/samtools/1.10


# Mapping the trimmed reads to undesired RNAs and filtering them out
# --------------------
mkdir -p ${DESTINATION}/filtered_data
cd ${DESTINATION}/filtered_data

for x in ${SOURCE}/*tr.fastq.gz; do gunzip ${x}; done # piping the gunzip -c ${x} output directly into bowtie does not work
for x in ${SOURCE}/*tr.fastq; \
do bowtie -t -p 4 ${INDICES}/GRCh38_p14_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb \
${x} --un $(basename ${x} .fastq)_no_r-t-sno-sn-RNA.fastq 2> $(basename ${x} .fastq.gz)_no_r-t-sno-sn-RNA_log.txt > /dev/null; done

for FILE in ${SOURCE}/*tr.fastq; do gzip ${FILE}; done

# Mapping the filtered reads to the genome
# --------------------
mkdir -p ${DESTINATION}/mapped_bamfiles
cd ${DESTINATION}/mapped_bamfiles

for x in ${DESTINATION}/filtered_data/*no_r-t-sno-sn-RNA.fastq; \
do bowtie -S -t -p 4 -v 1 -m 1 --best --strata ${INDICES}/GRCh38.p14.genome -q ${x} 2> $(basename ${x} .fastq)_GRCh38_p14_log.txt | \
samtools view -h -F 4 -b > $(basename ${x} .fastq)_GRCh38_p14.bam; done

# Sort the BAM file
for x in ${DESTINATION}/mapped_bamfiles/*.bam; do samtools sort -@ 4 ${x} -o $(basename ${x} .bam)_sorted.bam; done


# Remove the unsorted BAM file
rm *GRCh38_p14.bam

# Mapping filtered reads to transcriptome 
# --------------------
cd ${DESTINATION}/

mkdir -p codon_occupancy
cd codon_occupancy

for x in ${DESTINATION}/filtered_data/*no_r-t-sno-sn-RNA.fastq; \
do bowtie -t -p 4 -v 1 -m 1 --best --strata --norc ${INDICES}/GRCh38_p14_APPRIS_CDS_plus18 -q ${x} -S $(basename ${x} .fastq)_GRCh38_p14_APPRIS_CDS.sam 2> $(basename ${x} .fastq)_GRCh38_p14_APPRIS_CDS_log.txt; done