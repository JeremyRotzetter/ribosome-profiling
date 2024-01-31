#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=03:00:00
#SBATCH --job-name=prepare_indices
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/02_01_output_prepare_indices_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/02_01_error_prepare_indices_%j.e
#SBATCH --partition=pall

# Script will prepare indices using bowtie aligner for every annotation that was downloaded for the Ribo-seq project and can be found under annotation_...

# How to use: User needs to input path to the parent directory of where the genomic files are (SOURCE) (<parent/annotation_...>). The prepared indices will
# be stored in a new folder called annotation_indices, whose location must be specified.
# Example: sbatch 02_01_prepare_indices.sh </path/to/source/> 

# load bowtie 
module load UHTS/Aligner/bowtie/1.2.0
module load SequenceAnalysis/blat/36

SOURCE=$1
DESTINATION=$2/annotation_indices

# Use bowtie-build to make indices

mkdir -p ${DESTINATION}

# For the genome
# --------------------
cd ${SOURCE}/annotation_genome

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bowtie-build Homo_sapiens.GRCh38.dna.primary_assembly.fa ${DESTINATION}/GRCh38.p14.genome
faToTwoBit Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.p14.genome.2bit # done once for RiboseqQC analysis (after mapping to the genome)
gzip Homo_sapiens.GRCh38.dna.primary_assembly.fa


# For the "undesired" RNAs
# --------------------
cd ${SOURCE}/annotation_undesired_RNA

# Combine all the different files containing the different types of undesired RNA. This needs to be done only once.
cat *.txt > GRCh38_p14_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb.fa 

bowtie-build GRCh38_p14_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb.fa ${DESTINATION}/GRCh38_p14_r-t-sno-sn-RNA_ENSEMBL_NCBI_GtRNAdb


# For the transcriptome
# --------------------
cd ${SOURCE}/annotation_transcriptome

bowtie-build GRCh38_p14_APPRIS_CDS_plus18.fa ${DESTINATION}/GRCh38_p14_APPRIS_CDS_plus18


# Converting into single-line format for generating codon occupancy plots. This needs to be done only once.

awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < GRCh38_p14_APPRIS_CDS_plus18.fa > GRCh38_p14_APPRIS_CDS_plus18_SingleLine.fa

mv GRCh38_p14_APPRIS_CDS_plus18_SingleLine.fa ${DESTINATION}