################################################################################
######################### Differential expression analysis #####################
################################################################################

# This script will perform the differential gene expression analysis using
# DESeq2. It will:
# 1. Start with count matrix input to prepare the input
# (using CDS_counts_processed.txt obtained in previous step with
# 05_01_featureCounts.sh).
# 2. Perform pre-filtering.
# 3. Run the DESeq function and save its output to a variable.
# 4. Perform count data (log) transformation to remove the dependence of the
# variance on the mean.
# 5. Cluster samples based on their gene expression profiles.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
# register(MulticoreParam(8)) # Change this based on your computer core count
# NOTE: MulticoreParam() not supported on Windows


# Setting up color profiles from colorbrewer
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
mycols <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_1 <- "WT"
sample_2 <- "KO"
no_of_reps <- 2

project_name <- paste("RPF", sample_2, "vs", sample_1, sep = "_")

## Create Samples dataframe
condition_column <- c(
  rep(sample_1, no_of_reps),
  rep(sample_2, no_of_reps)
)

run_column <- c(
  paste(sample_1, "1", sep = "_"),
  paste(sample_1, "2", sep = "_"),
  paste(sample_2, "1", sep = "_"),
  paste(sample_2, "2", sep = "_")
)

rep_column <- c(
  "A", "B",
  "A", "B"
)

samples_df <- data.frame(
  condition_column,
  run_column,
  rep_column
)

colnames(samples_df) <- c("condition", "run", "rep")

rownames(samples_df) <- samples_df$run

samples_df$condition <- factor(rep(c(
  rep(sample_1, no_of_reps),
  rep(sample_2, no_of_reps)
)))

featurecount_data <- read.table("CDS_counts_processed.txt", header = TRUE, row.names = 1)

# Reorder to make the order consistent with samples$run
featurecount_data <- featurecount_data[, c(3, 4, 1, 2)]

# Change colnames
colnames(featurecount_data) <- rownames(samples_df)

# Import as DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(
  countData = featurecount_data,
  colData = samples_df,
  design = ~condition
)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Factor levels
# WT columns should be first followed by KO/treatment
dds$condition <- factor(
  c(
    rep(sample_1, no_of_reps),
    rep(sample_2, no_of_reps)
  ),
  levels = c(
    sample_1,
    sample_2
  )
)

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
  margin = c(10, 10), main = "Sample Distance Matrix"
)
dev.off()

# Count matrix heatmap
select <- order(rowMeans(counts(dds, normalized = TRUE)))
df <- as.data.frame(colData(dds)[, c("condition", "rep")])

pdf(paste(project_name, "QC_count_matrix_CDS.pdf", sep = "_"))
pheatmap(assay(rld)[select, ],
  cluster_rows = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  annotation_col = df
)
dev.off()


# PCA plot

pdf(paste(project_name, "QC_PCA_CDS.pdf", sep = "_"))
pcaData <- plotPCA(rld, intgroup = c("condition", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(
  pcaData,
  aes(PC1,
    PC2,
    color = condition,
    shape = rep
  )
) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
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

convertIDs <- function(ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(
    db,
    keys = ids, keytype = from, columns = c(from, to)
  ))
  if (ifMultiple == "putNA") {
    duplicatedIds <- selRes[duplicated(selRes[, 1]), 1]
    selRes <- selRes[!selRes[, 1] %in% duplicatedIds, ]
  }
  return(selRes[match(ids, selRes[, 1]), 2])
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
  ifelse((res_2_vs_1_df$log2FoldChange < -0.5 & res_2_vs_1_df$padj < 0.05), "Downregulated",
    "Unchanged"
  )
)

write.table(res_2_vs_1_df,
  file = paste(project_name, "DESeq2_res.csv", sep = "_"),
  sep = ",",
  row.names = F,
  col.names = T,
  quote = F
)

res_2_vs_1_df$regulation_level <- factor(res_2_vs_1_df$regulation_level, levels = c("Upregulated", "Downregulated", "Unchanged"))

## Remove rows with NA in padj column before plotting
res_2_vs_1_df <- res_2_vs_1_df[!is.na(res_2_vs_1_df$padj), ]

pdf(paste(project_name, "Volcano_plot1.pdf", sep = "_"), width = 4, height = 5)
ggplot(
  res_2_vs_1_df,
  aes(
    x = -log10(padj),
    y = log2FoldChange,
    color = regulation_level
  )
) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
  xlab("-log10(adjusted p-value)") +
  ylab("Log2 fold change") +
  labs(color = "Regulation level") +
  theme_bw()
dev.off()


pdf(paste(project_name, "Volcano_plot2.pdf", sep = "_", paper = "a4"), width = 8, height = 4)
ggplot(
  res_2_vs_1_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = regulation_level
  )
) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
  xlab("Log2 fold change") +
  ylab("-log10(adjusted p-value)") +
  labs(color = "Regulation level") +
  theme_bw()
dev.off()
