################################################################################
######################### Ribo-seQC ############################################
################################################################################

# Ribo-seQC is a package that performs quality control analysis of small RNA-seq
# data (in .bam format), with a focus on Ribo-seq and related techniques.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# will set working directory to the script file path in Rstudio
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Installation of the package (to be done only once!)
# Despite having bugs that need fixing, should be done to install the
# dependencies.

# install.packages("devtools")
#
# library("devtools")
#
# install_github(repo = "ohlerlab/RiboseQC")

###### Debugging part ######

# 1. Download the fork of RiboseQC either by git clone or direct download of zip
# file from here: https://github.com/Tim-Yu/RiboseQC/tree/master

# 2. Open RiboseQC/R/riboseqc.R in RStudio (or any other IDE)

# 3. Comment out line 2283 i.e. it should look like
# # genome_sequence<-get(library(GTF_annotation$genome,character.only = TRUE))
# i.e. put a # at the starting of this line.

# 4. Edit line 2694 from n_transcripts = length(unique(gtfdata$transcript_id)) to
# n_transcripts = length(unique(na.omit(gtfdata$transcript_id)))

# 5. Save your changes

# 6. Press Cmd/Ctrl + Shift + 0 or Session > Restart R

# 7. Type remove.packages("RiboseQC") to remove your current installation of
# RiboseQC that has the bug.

# 8. In the terminal, navigate to the directory which contains correct RiboseQC code

# 9. Issue the command (in terminal) R CMD build RiboseQC/ which should make
# RiboseQC_0.99.0.tar.gz

# 10. Open RStudio and install.packages("/path/to/RiboseQC_0.99.0.tar.gz", repos = NULL, type="source")

# 11. Try to do the QC analysis again.

###### Analysis part ######


# Load the package
library(RiboseQC)


# Prepare genome file (to be done only once!!!)
prepare_annotation_files(
  annotation_directory = "./",
  twobit_file = "GRCh38.p14.genome.2bit",
  gtf_file = "Homo_sapiens.GRCh38.108.gtf",
  scientific_name = "Homo.sapiens",
  annotation_name = "GRCh38",
  export_bed_tables_TxDb = F,
  forge_BSgenome = T,
  create_TxDb = T
)


# Genome mapped sorted-BAM files

genome_bam <- c(
  "RPF_KO_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_sorted.bam",
  "RPF_KO_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_sorted.bam",
  "RPF_WT_Rep1_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_sorted.bam",
  "RPF_WT_Rep2_clpd_tr_no_r-t-sno-sn-RNA_GRCh38_p14_sorted.bam"
)

load_annotation("Homo_sapiens.GRCh38.108.gtf_Rannot")

###### QC plots ######

RiboseQC_analysis(
  annotation_file = "Homo_sapiens.GRCh38.108.gtf_Rannot",
  bam_files = genome_bam,
  fast_mode = T,
  report_file = "RPF_samples_QC.html",
  sample_names = c(
    "WT_Rep1", "WT_Rep2",
    "KO_Rep1", "KO_Rep2"
  ),
  dest_names = c(
    "WT_Rep1", "WT_Rep2",
    "KO_Rep1", "KO_Rep2"
  ),
  write_tmp_files = F
)
