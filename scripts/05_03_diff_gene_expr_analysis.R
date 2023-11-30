# Open explorer window to select the working directory (assumes script is run
# in Rstudio, requires graphical interactive desktop environment)
setwd(rstudioapi::selectDirectory())

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18") # to upgrade to the at the time of
# writing most recent version

# BiocManager::install(c("org.Hs.eg.db", "pathview", "clusterProfiler", "topGO"))


library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
library(org.Hs.eg.db) # Human AnnotationDbi

library(grid)
library(gridExtra)
library(pathview)
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)

library(tidyverse)

# register(MulticoreParam(8)) # Change this based on your computer core count
# NOTE: MulticoreParam() not supported on Windows

# Setting up color profiles from colorbrewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_1 <- "WT"
sample_2 <- "KO"
no_of_reps <- 2

project_name <- paste("RPF", sample_2, "vs", sample_1, sep = "_")

## Create Samples dataframe
condition_column <- c(rep(sample_1, no_of_reps),
                      rep(sample_2, no_of_reps))

run_column <- c(paste(sample_1, "1", sep = "_"),
                paste(sample_1, "2", sep = "_"),
                paste(sample_2, "1", sep = "_"),
                paste(sample_2, "2", sep = "_"))

rep_column <- c("A", "B",
                "A", "B")

samples_df <- data.frame(condition_column,
                         run_column,
                         rep_column)

colnames(samples_df) <- c("condition", "run", "rep")

rownames(samples_df) <- samples_df$run

samples_df$condition <- factor(rep(c(rep(sample_1, no_of_reps),
                                     rep(sample_2, no_of_reps))))

featurecount_data <- read.table("CDS_counts_processed.txt", header = TRUE, row.names = 1)

# Reorder to make the order consistent with samples$run
featurecount_data <- featurecount_data[ , c(3, 4, 1, 2)]

# Change colnames
colnames(featurecount_data) <- rownames(samples_df)

# Import as DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = featurecount_data,
                              colData = samples_df,
                              design = ~ condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Factor levels
# WT columns should be first followed by KO/treatment
dds$condition <- factor(c(rep(sample_1, no_of_reps),
                          rep(sample_2, no_of_reps)),
                        levels = c(sample_1,
                                   sample_2))

# Differential expression analysis
dds <- DESeq(dds)


# colData(dds) # to check whether names are correct

################################################################################
################################################################################
# QC
################################################################################
################################################################################

# Log transformation for data quality assessment
rld <- rlog(dds, blind = FALSE)

# Sample distance matrix
sampleDists <- as.matrix(dist(t(assay(rld))))
pdf(paste(project_name, "QC_sample_distance_matrix_CDS.pdf", sep = "_"))
heatmap.2(as.matrix(sampleDists),
          key = T,
          trace = "none",
          col = colorpanel(100, "#2b8cbe", "#7bccc4", "white"),
          ColSideColors = mycols[dds$condition],
          RowSideColors = mycols[dds$condition],
          margin = c(10, 10), main = "Sample Distance Matrix")
dev.off()

# Count matrix heatmap
select <- order(rowMeans(counts(dds, normalized = TRUE)))
df <- as.data.frame(colData(dds)[ , c("condition","rep")])

pdf(paste(project_name, "QC_count_matrix_CDS.pdf", sep = "_"))
pheatmap(assay(rld)[select, ],
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = df)
dev.off()


# PCA plot

pdf(paste(project_name, "QC_PCA_CDS.pdf", sep = "_"))
pcaData <- plotPCA(rld, intgroup = c("condition", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData,
       aes(PC1,
           PC2,
           color = condition,
           shape = rep)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

################################################################################
################################################################################
## Resume the analysis
################################################################################
################################################################################

res_2_vs_1 <- results(dds, contrast = c("condition", sample_2, sample_1), alpha = 0.05)

# Adding gene names using org.Hs.eg.db
# Source: http://bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html
# Also: https://support.bioconductor.org/p/66288/
# This function takes a list of IDs as first argument and their key type as the second argument.
# The third argument is the key type we want to convert to, the fourth is the AnnotationDb object to use.
# Finally, the last argument specifies what to do if one source ID maps to several target IDs:
# should the function return an NA or simply the first of the multiple IDs

convertIDs <- function( ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(
    db, keys = ids, keytype = from, columns = c(from, to)))
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[duplicated(selRes[ , 1] ), 1]
    selRes <- selRes[!selRes[ , 1] %in% duplicatedIds, ]
  }
  return(selRes[match(ids, selRes[ , 1]), 2])
}

# # Check columns in the database:
columns(org.Hs.eg.db)

# Actual adding of the column

res_2_vs_1$GeneID <- row.names(res_2_vs_1)
res_2_vs_1$gene_symbol <- convertIDs(row.names(res_2_vs_1), "ENSEMBL", "SYMBOL", org.Hs.eg.db)


## Data summary
summary(res_2_vs_1)

## Convert DESeq object to dataframe
res_2_vs_1_df <- as.data.frame(res_2_vs_1)

## Define differentially expressed genes as:
## Upregulated (LFC > 0.5 & padj < 0.05)
## Downregulated (LFC < - 0.5 & padj < 0.05)

res_2_vs_1_df$regulation_level <- ifelse((res_2_vs_1_df$log2FoldChange > 0.5 & res_2_vs_1_df$padj < 0.05), "Upregulated",
                                         ifelse((res_2_vs_1_df$log2FoldChange < - 0.5 & res_2_vs_1_df$padj < 0.05), "Downregulated",
                                                "Unchanged"))

write.table (res_2_vs_1_df,
             file = paste(project_name, "DESeq2_res.csv", sep = "_"),
             sep = ",",
             row.names = F,
             col.names = T,
             quote = F)

res_2_vs_1_df$regulation_level <- factor(res_2_vs_1_df$regulation_level, levels = c("Upregulated", "Downregulated", "Unchanged"))

## Remove rows with NA in padj column before plotting
res_2_vs_1_df <- res_2_vs_1_df[!is.na(res_2_vs_1_df$padj), ]

pdf(paste(project_name, "Volcano_plot1.pdf", sep = "_"), width = 4, height = 5)
ggplot(res_2_vs_1_df,
       aes(x = -log10(padj),
           y = log2FoldChange,
           color = regulation_level)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
  xlab("-log10(adjusted p-value)") +
  ylab("Log2 fold change") +
  labs(color = "Regulation level") +
  theme_bw()
dev.off()


pdf(paste(project_name, "Volcano_plot2.pdf", sep = "_", paper = "a4"), width = 8, height = 4)
ggplot(res_2_vs_1_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = regulation_level)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
  xlab("Log2 fold change") +
  ylab("-log10(adjusted p-value)") +
  labs(color = "Regulation level") +
  theme_bw()
dev.off()

########################################################################################
######################### Gene Ontology Analysis #######################################
########################################################################################

## This part of the script is designed to take the DESeq2 result file as input.
## This file was generated in the previous step (separate script) and contains
## GeneID column containing Ensembl gene ids.
## If you are using your own script for DGE analysis then make sure you have
## the following columns:
## GeneID, log2FoldChange, padj

sample_name = "RPF_KO_vs_WT"

## Load DESeq2 output file.

df <- read.csv("RPF_KO_vs_WT_DESeq2_res.csv", sep = ",", header = T)

## Set rownames as GeneID column

rownames(df) <- df$GeneID
df <- df[order(df$padj), ]

## Define Up and Down regulated genes
## Here we are considering padj values only so that we do not lose important terms

genes_up <- which(df$padj < 0.05 & df$log2FoldChange > 0)
genes_down <- which(df$padj < 0.05 & df$log2FoldChange < 0)

all_genes_names <- rownames(df)

genes_up <- rownames(df)[genes_up]
genes_down <- rownames(df)[genes_down]

genelist_up <- factor(as.integer(all_genes_names %in% genes_up))
names(genelist_up) <- all_genes_names

genelist_down <- factor(as.integer(all_genes_names %in% genes_down))
names(genelist_down) <- all_genes_names


## Map IDs
allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Hs.eg.db",
                          ID = "ensembl")

## Function

topgo_object <- function(ontology, gene_list) {
  
  #' Create TopGO object
  
  object <- new("topGOdata",
                ontology = ontology,
                allGenes = gene_list,
                annot = annFUN.GO2genes,
                GO2genes = allGO2genes,
                nodeSize = 10)
  object
  
}

## Upregulates genes
GOdata_up_bp <- topgo_object("BP", genelist_up)
GOdata_up_mf <- topgo_object("MF", genelist_up)
GOdata_up_cc <- topgo_object("CC", genelist_up)

## Downregulated genes
GOdata_down_bp <- topgo_object("BP", genelist_down)
GOdata_down_mf <- topgo_object("MF", genelist_down)
GOdata_down_cc <- topgo_object("CC", genelist_down)

# Fischer exact test

resultFis_up_bp <- runTest(GOdata_up_bp, algorithm = "elim", statistic = "fisher")
resultFis_up_mf <- runTest(GOdata_up_mf, algorithm = "elim", statistic = "fisher")
resultFis_up_cc <- runTest(GOdata_up_cc, algorithm = "elim", statistic = "fisher")
resultFis_down_bp <- runTest(GOdata_down_bp, algorithm = "elim", statistic = "fisher")
resultFis_down_mf <- runTest(GOdata_down_mf, algorithm = "elim", statistic = "fisher")
resultFis_down_cc <- runTest(GOdata_down_cc, algorithm = "elim", statistic = "fisher")

## Function to extract data from S4 objects

parse_tables <- function(GO_data, statistics) {
  
  #' Parse TopGO data
  
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
  
}

GOres_up_bp <- parse_tables(GOdata_up_bp, resultFis_up_bp)
GOres_up_mf <- parse_tables(GOdata_up_mf, resultFis_up_mf)
GOres_up_cc <- parse_tables(GOdata_up_cc, resultFis_up_cc)

GOres_down_bp <- parse_tables(GOdata_down_bp, resultFis_down_bp)
GOres_down_mf <- parse_tables(GOdata_down_mf, resultFis_down_mf)
GOres_down_cc <- parse_tables(GOdata_down_cc, resultFis_down_cc)

## Function
plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  
  #' Make lollipop plots
  
  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  
  ggplot(GO_data, 
         aes(x = log_weightFisher,
             y = Term)) +
    geom_segment(aes(x = 0,
                     xend = log_weightFisher,
                     y = Term,
                     yend = Term),
                 colour = use_color)  +
    geom_point(aes(size = Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"))
}

## Make Plots

plot_up_BP <- plot_GO(GOres_up_bp, "Biological Proccess", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_MF <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_CC <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "#404788FF")


plot_down_BP <- plot_GO(GOres_down_bp, "Biological Proccess", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_MF <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_CC <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "#73D055FF")

# Prepare data for writing
# retrieve genes2GO list from the "expanded" annotation in GOdata

extract_genes <- function(GOdata, regulation_level) {
  
  #' Extract our genes in all identified GO terms
  
  all_GO_dereg <- genesInTerm(GOdata)
  data_to_extract <- lapply(all_GO_dereg, function(x) x[x %in% regulation_level])
  no_of_observations <- sapply(data_to_extract, length)
  seq_max <- seq_len(max(no_of_observations))
  GO_df <- t(sapply(data_to_extract, "[", i = seq_max))
  GO_df <- data.frame(GO_df[complete.cases(GO_df[ , 1:2]), ])
  GO_df$GO_ID <- row.names(GO_df)
  GO_df <- pivot_longer(GO_df,
                        cols = - GO_ID,
                        names_to = "Column_no",
                        values_to = "GeneID")
  GO_df <- GO_df[ , c(1, 3)]
  GO_df <- na.omit(GO_df)
  GO_df
  
}


GO_df_up_bp <- extract_genes(GOdata_up_bp, genes_up)
GO_df_up_mf <- extract_genes(GOdata_up_mf, genes_up)
GO_df_up_cc <- extract_genes(GOdata_up_cc, genes_up)

GO_df_down_bp <- extract_genes(GOdata_down_bp, genes_down)
GO_df_down_mf <- extract_genes(GOdata_down_mf, genes_down)
GO_df_down_cc <- extract_genes(GOdata_down_cc, genes_down)

## Export plots as PDF

pdf(paste(sample_name, "_Biological_Proccess_TopGO_Up_fisher.pdf", sep = ""))
plot_up_BP
dev.off()

pdf(paste(sample_name, "_Molecular_Function_TopGO_Up_fisher.pdf", sep = ""))
plot_up_MF
dev.off()

pdf(paste(sample_name, "_Cellular_Component_TopGO_Up_fisher.pdf", sep = ""))
plot_up_CC
dev.off()

pdf(paste(sample_name, "_Biological_Proccess_TopGO_down_fisher.pdf", sep = ""))
plot_down_BP
dev.off()

pdf(paste(sample_name, "_Molecular_Function_TopGO_down_fisher.pdf", sep = ""))
plot_down_MF
dev.off()

pdf(paste(sample_name, "_Cellular_Component_TopGO_down_fisher.pdf", sep = ""))
plot_down_CC
dev.off()

# Tables
write.table(GO_df_up_bp,
            file = paste(sample_name, "_BP_Up_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_up_mf,
            file = paste(sample_name, "_MF_Up_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_up_cc,
            file = paste(sample_name, "_CC_Up_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_down_bp,
            file = paste(sample_name, "_BP_Down_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_down_mf,
            file = paste(sample_name, "_MF_Down_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_down_cc,
            file = paste(sample_name, "_CC_Down_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

## GO data

write.table(GOres_up_bp,
            file = paste(sample_name, "_BP_Up_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_up_mf,
            file = paste(sample_name, "_MF_Up_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_up_cc,
            file = paste(sample_name, "_CC_Up_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_down_bp,
            file = paste(sample_name, "_BP_Down_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_down_mf,
            file = paste(sample_name, "_MF_Down_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_down_cc,
            file = paste(sample_name, "_CC_Down_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")