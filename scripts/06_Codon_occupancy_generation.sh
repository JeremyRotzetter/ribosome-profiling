#!/usr/bin/env bash

#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=codon_occupancy_generation
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/06_output_codon_occupancy_generation_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/06_error_codon_occupancy_generation_%j.e
#SBATCH --partition=pall

# This script will perform codon occupancy calculations by calling Codon_occupancy_cal.sh and using SAM files obtained
# with 03_mapping.sh as input to determine the A-site ribosome occupancy that is normalized by positions +5, +6 and +7.

CDS_FA=$1
MAP_SAM=$2
SCRIPT=$3
OUTPUT=$4

cd ${OUTPUT}

${SCRIPT}/Codon_occupancy_cal.sh \
${CDS_FA}/GRCh38_p14_APPRIS_CDS_plus18_SingleLine.fa \
${MAP_SAM}/RPF_KO_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep1_Codon_occupancy.txt

${SCRIPT}/Codon_occupancy_cal.sh \
${CDS_FA}/GRCh38_p14_APPRIS_CDS_plus18_SingleLine.fa \
${MAP_SAM}/RPF_KO_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep2_Codon_occupancy.txt

${SCRIPT}/Codon_occupancy_cal.sh \
${CDS_FA}/GRCh38_p14_APPRIS_CDS_plus18_SingleLine.fa \
${MAP_SAM}/RPF_WT_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep1_Codon_occupancy.txt

${SCRIPT}/Codon_occupancy_cal.sh \
${CDS_FA}/GRCh38_p14_APPRIS_CDS_plus18_SingleLine.fa \
${MAP_SAM}/RPF_WT_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep2_Codon_occupancy.txt