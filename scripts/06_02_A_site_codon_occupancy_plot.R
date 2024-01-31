######################### A-site codon occupancy plot ##########################
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This script will combine the data obtained with
# 06_01_Codon_occupancy_generation.sh and generate the A-site occupancy plot.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Open explorer window to select the working directory (assumes script is run
## in Rstudio, requires graphical interactive desktop environment)
## setwd(rstudioapi::selectDirectory())

## Alternatively:
# setwd("/some/path/")

## Load package
library(tidyverse)

## Datasets
df_1 <- read.table("./codon_occupancy/RPF_WT_Rep1_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")

df_2 <- read.table("./codon_occupancy/RPF_WT_Rep2_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")

df_3 <- read.table("./codon_occupancy/RPF_KO_Rep1_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")

df_4 <- read.table("./codon_occupancy/RPF_KO_Rep2_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")

## Function to add column names
process_data <- function(input_df, sample) {
  colnames(input_df) <- c("Codon", "Occupancy")
  input_df$Sample <- sample
  input_df
}

df_1 <- process_data(df_1, "WT_Rep1")
df_2 <- process_data(df_2, "WT_Rep2")
df_3 <- process_data(df_3, "KO_Rep1")
df_4 <- process_data(df_4, "KO_Rep2")

## Merging datasets
data_wt <- inner_join(df_1, df_2, by = "Codon", suffix = c("_1", "_2"))
data_ko <- inner_join(df_3, df_4, by = "Codon", suffix = c("_1", "_2"))

## Take average of the two samples
data_wt$WT <- (data_wt$Occupancy_1 + data_wt$Occupancy_2) / 2
data_ko$KO <- (data_ko$Occupancy_1 + data_ko$Occupancy_2) / 2

## Extract the relevant columns
data_wt <- data_wt[ , c(1, 6)]
data_ko <- data_ko[ , c(1, 6)]

## Merge KO and WT datasets
data_plot <- inner_join(data_wt, data_ko, 
                        by = "Codon", 
                        suffix = c("", ""))

## Calculate normalized occupancy
data_plot$Norm_Occupancy_KO <- data_plot$KO / data_plot$WT

## Extract the relevant columns
data_plot <- data_plot[ , c(1, 4)]

## Plot
pdf("./codon_occupancy/A_site_codon_occupancy.pdf",width = 8, height = 10)
ggplot(data_plot,
       aes(x = Norm_Occupancy_KO,
           y = Codon)) +
  geom_point(color = "#dd1c77") +
  xlab("Codon occupancy (KO/WT)") +
  ylab("Codons") +
  theme_bw()
dev.off()