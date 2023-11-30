#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --job-name=fastqc
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/01_output_fastqc_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/01_error_fastqc_%j.e
#SBATCH --partition=pall

# load the needed module(s)
module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Analysis/MultiQC/1.8

# define working, data and output directories
WORKDIR=/data/users/jrotzetter/ribosome-profiling
READS_DIR=/data/users/jrotzetter/ribosome-profiling/raw_reads
OUTPUT_DIR=/data/users/jrotzetter/ribosome-profiling/qc

mkdir -p "$OUTPUT_DIR"

cd ${READS_DIR}

# perform fastqc on all fastq.gz present in the $READS_DIR
fastqc -o ${OUTPUT_DIR} *fastq.gz

cd ${OUTPUT_DIR}

# MultiQC searches a given directory for analysis logs and compiles a HTML report.
# It's a general use tool, perfect for summarising the output from numerous bioinformatics tools. https://multiqc.info/
multiqc -o ${OUTPUT_DIR} .