#!/usr/bin/env bash

#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=6
#SBATCH --job-name=ex1
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/output_ex1_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/error_ex1_%j.e
#SBATCH --partition=pall

WORKDIR=/data/users/jrotzetter/ribosome-profiling/raw_reads

cd ${WORKDIR}

module load UHTS/Analysis/sratoolkit/2.10.7

# fasterq-dump will automatically download the files from the specified accession numbers and convert them to fastq format.
fasterq-dump --progress SRR9596295
fasterq-dump --progress SRR9596296
fasterq-dump --progress SRR9596300
fasterq-dump --progress SRR9596303
fasterq-dump --progress SRR9596304
fasterq-dump --progress SRR9596310

# create a compressed file for any file ending in .fastq
gzip *.fastq

# rename the files
mv ./SRR9596295.fastq.gz ./Biever_Somata_Poly_1.fastq.gz
mv ./SRR9596296.fastq.gz ./Biever_Somata_Poly_2.fastq.gz
mv ./SRR9596300.fastq.gz ./Biever_Somata_Poly_3.fastq.gz
mv ./SRR9596310.fastq.gz ./Biever_Neuropil_Poly_1.fastq.gz
mv ./SRR9596303.fastq.gz ./Biever_Neuropil_Poly_2.fastq.gz
mv ./SRR9596304.fastq.gz ./Biever_Neuropil_Poly_3.fastq.gz