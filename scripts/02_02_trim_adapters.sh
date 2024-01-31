#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --job-name=trim_adapters
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/02_02_output_trim_adapters_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/02_02_error_trim_adapters_%j.e
#SBATCH --partition=pall

# Data Pre-processing
# Script will trim the raw reads to remove the adatper sequences and perform fastqc on the trimmed reads afterwards.

# How to use: User needs to input path to directory where the input data is (SOURCE) and the path to the directory where output data should go (DESTINATION) when running the script with sbatch
# Example: sbatch 01_02_trim_adapters.sh </path/to/source/> </path/to/destination/> 

# </path/to/source/> where the input files are present (in this case the raw reads).
SOURCE=$1 

# </path/to/destination/> where the processed data should be stored.
DESTINATION=$2


# Load all required modules

module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Quality_control/cutadapt/2.5
module add UHTS/Analysis/MultiQC/1.8


mkdir -p ${DESTINATION}/trimmed_reads
cd ${DESTINATION}/trimmed_reads

# Command for clipping 3' adapter
# 1> is used to send the standard output of the code to the filename that follows

for FILE in ${SOURCE}/*.fastq.gz;
do cutadapt -j 4 -q 25 -a CTGTAGGCACCATCAAT -m 25 --discard-untrimmed -o $(basename ${FILE} .fastq.gz)_clpd.fastq.gz \
${FILE} 1> $(basename ${FILE} .fastq.gz)_clpd_log.txt;
done

# for trimming 4 nt (randomized bases) from the 3' end

for i in *_clpd.fastq.gz;
do cutadapt -j 4 -q 25 --cut -4 -m 25 -o $(basename ${i} .fastq.gz)_tr.fastq.gz ${i} 1> $(basename ${i} .fastq.gz)_tr_log.txt;
done

# QC analysis of the processed reads

mkdir -p ${DESTINATION}/QC_trimmed_reads
cd ${DESTINATION}/QC_trimmed_reads

fastqc -t 32 -o ${DESTINATION}/QC_trimmed_reads ${DESTINATION}/trimmed_reads/*_tr.fastq.gz

multiqc .

# The fastqc and multiqc reports can now be downloaded to a local computer to assess the quality of the trimmed reads.