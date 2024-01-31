######################### codon decoding speed plot ############################
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This script will combine the data obtained with
# 06_01_Codon_occupancy_generation.sh and generate the codon decoding speed plot.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Open explorer window to select the working directory (assumes script is run
## in Rstudio, requires graphical interactive desktop environment)
# setwd(rstudioapi::selectDirectory())

# setwd("/some/path/")

## Load package
library(data.table)
library(tidyr)
library(reshape2)
library(devtools)
library(ggplot2)

## Datasets
df_1 <- read.table("./codon_occupancy/RPF_WT_Rep1_Codon_occupancy.txt",
                   header = F,
                   sep = "\t") %>% as.data.table()

df_2 <- read.table("./codon_occupancy/RPF_WT_Rep2_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")%>% as.data.table()

df_3 <- read.table("./codon_occupancy/RPF_KO_Rep1_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")%>% as.data.table()

df_4 <- read.table("./codon_occupancy/RPF_KO_Rep2_Codon_occupancy.txt",
                   header = F,
                   sep = "\t")%>% as.data.table()

## Function to add column names, condition and replicate
process_data <- function(input_df, Condition, Replicate) {
  colnames(input_df) <- c("Codon", "Occupancy")
  input_df$Condition <- Condition
  input_df$Replicate <- Replicate
  input_df
}

df_1 <- process_data(df_1, "WT","Rep_1")
df_2 <- process_data(df_2, "WT","Rep_2")
df_3 <- process_data(df_3, "KO","Rep_1")
df_4 <- process_data(df_4, "KO","Rep_2")

## Combine all the datasets by row
df_all = rbind(df_1, df_2) %>% rbind(.,df_3) %>% rbind(.,df_4)


## Separate the samples by condition and replicate for each codon
df_all_melt = dcast(df_all, Codon + Condition ~ Replicate, 
                    value.var = 'Occupancy') %>% as.data.table()
colnames(df_all_melt)[3:4] = c('Rep_1','Rep_2')

## Calculate mean of occupancy values for replicates in WT condition for all codons
df_all_melt[,WT_mean:=mean(c(Rep_1[Condition=='WT'], 
                             Rep_2[Condition=='WT'])), by = 'Codon']

## For each codon calculate 1.) values relative to mean WT values for all
## replicates in all conditions (min and max value) and 2.) values relative to
## mean WT values for mean values of both conditions (average)
df_all_melt[,":="(min = min(c(Rep_1,Rep_2))/WT_mean,
                  average = mean(c(Rep_1,Rep_2))/WT_mean,
                  max = max(c(Rep_1,Rep_2))/WT_mean),
            by=c('Codon','Condition')]

## Order by mean value for KO condition
df_all_melt = df_all_melt[order(average)]

## Designate codon as factor 
df_all_melt$Codon = factor(df_all_melt$Codon, level = df_all_melt[Condition=='KO']$Codon)

## take the 25 % and 75 % quantile of (mean) occupancy seen in WT as min and max respectively.
xmin <- quantile(df_all_melt$WT_mean[df_all_melt$Condition == "WT"], 0.25)
xmax <- quantile(df_all_melt$WT_mean[df_all_melt$Condition == "WT"], 0.75)

## Plot
ggplot(df_all_melt, aes(y=average, x=Codon, colour=Condition, fill=Condition)) + 
  geom_point(position=position_dodge(width=0.3), size=1) +
  geom_errorbar(aes(ymin=min, ymax=max), position = position_dodge(0.3)) + 
  # geom_hline(yintercept= c(0.9,1.1),linetype=3) + # arbitrary cutoff lines
  theme_bw(base_size=16) + coord_flip() +
  ylab("A-site codon occupancy (normalized to average WT)")
  # geom_rect(aes(xmin = Inf, # something isn't right, commented for now
  #               xmax = -Inf,
  #               ymin = xmin,
  #               ymax = xmax),
  #           inherit.aes = FALSE,
  #           fill = "gray", alpha = 0.01)
ggsave("codon_decoding_speed.svg", device = "svg", path = "./codon_occupancy/",
       width = 310, height = 240, dpi = 300, units = "mm")
