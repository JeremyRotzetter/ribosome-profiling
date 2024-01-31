#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --job-name=featureCounts
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/05_01_output_featureCounts_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/05_01_error_featureCounts_%j.e
#SBATCH --partition=pall

# This script will generate a table containing raw read counts for each gene (in rows) for each sample (in column)

module load UHTS/Analysis/subread/2.0.1

GENOME_ANNOT=$1 # directory where the gene transfer file is located
SORTED_BAM=$2 # directory where sorted bam files are located
DESTINATION=$3 # directory where output should go

# change directory to source where the gene transfer file is located
cd ${GENOME_ANNOT}

# unzip the annotation file for genome
gunzip Homo_sapiens.GRCh38.108.gtf.gz

# change directory to DESTINATION
cd ${DESTINATION}

# make a directory to store counts table and change to that directory
mkdir -p counts_table
cd counts_table

# use featureCounts to map the bam files to the coding sequence in the genome annotation and collapse to yield one value per gene 
featureCounts -T 8 -t CDS -g gene_id -a ${GENOME_ANNOT}/Homo_sapiens.GRCh38.108.gtf -o CDS_counts_rawfile.txt ${SORTED_BAM}/*_GRCh38_p14_sorted.bam

# use featureCounts to map the bam files to the exonic sequence in the genome annotation and collapse to yield one value per biotype 
featureCounts -T 8 -t exon -g gene_biotype -a ${GENOME_ANNOT}/Homo_sapiens.GRCh38.108.gtf -o biotype_counts_rawfile.txt ${SORTED_BAM}/*_GRCh38_p14_sorted.bam

# extract columns 1,7,8,9,10 from the raw counts table for further analysis
cut -f 1,7-10 CDS_counts_rawfile.txt > CDS_counts_processed.txt
cut -f 1,7-10 biotype_counts_rawfile.txt > biotype_counts_processed.txt

# change to source DESTINATION and compress the gene transfer file again
cd ${GENOME_ANNOT}
gzip Homo_sapiens.GRCh38.108.gtf