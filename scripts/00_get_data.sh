#!/usr/bin/env bash

#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=wget
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/00_output_wget_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/00_error_wget_%j.e
#SBATCH --partition=pall

WORKDIR=/data/users/jrotzetter/ribosome-profiling/raw_reads

cd ${WORKDIR}

# will download the files from the given URL and name them based on the given <FILE-NAME.FORMAT>. The links below expired on 22nd November 2023.
wget -O RPF_KO_Rep1.fastq.gz 'https://filesender.switch.ch/filesender2/download.php?token=7a5bdf5e-2b35-432a-a90d-0e982cc2f442&files_ids=525885'
wget -O RPF_KO_Rep2.fastq.gz 'https://filesender.switch.ch/filesender2/download.php?token=7a5bdf5e-2b35-432a-a90d-0e982cc2f442&files_ids=525886'
wget -O RPF_WT_Rep1.fastq.gz 'https://filesender.switch.ch/filesender2/download.php?token=7a5bdf5e-2b35-432a-a90d-0e982cc2f442&files_ids=525887'
wget -O RPF_WT_Rep2.fastq.gz 'https://filesender.switch.ch/filesender2/download.php?token=7a5bdf5e-2b35-432a-a90d-0e982cc2f442&files_ids=525888'