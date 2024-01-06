# Set working directory
getwd()

# Clear environment
rm(list=ls())
# Uncomment the following line if you need to install the "BiocManager" package
# Install and load necessary libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
          #install.packages("BiocManager")

# Uncomment the following line if you need to install the "DESeq2" package
#if (!requireNamespace("DESeq2", quietly = TRUE))
          #iocManager::install("DESeq2")
library(DESeq2)

# Uncomment the following line if you need to install the "ComplexHeatmap" package
# BiocManager::install("ComplexHeatmap")
          
# Uncomment the following line if you need to install the "org.Hs.eg.db" package
# BiocManager::install("org.Hs.eg.db")
          
# Load additional libraries or packages if needed
         
          
# Install and load other libraries
install.packages(c("ggplot2", "dplyr", "ashr", "org.Hs.eg.db", "annotables", "ComplexHeatmap"))
library(ggplot2)
library(dplyr)
library(ashr)
library(org.Hs.eg.db)
library(annotables)
library(ComplexHeatmap)



# Read feature counts data
f <- read.table("D:/shared_D_folder/PitPat_analysis/featureCountsP_file/featureCounts_resultsP.txt", header = TRUE)

# Subset data and remove unnecessary columns
dataRow <- subset(f[c(-2:-6), ], select = c(1, 8, 9, 10, 11, 12, 13, 14, 15))
colnames(dataRow) <- c(" ", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")

# Create metadata
condition <- c("growth-hormone", "growth-hormone", "non-functioning", "non-functioning", "non-functioning", "non-functioning", "growth-hormone", "non-functioning")
Sample.names <- c("t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")
pitpat_metadata <- data.frame(Sample.names, condition)
row.names(pitpat_metadata) <- c("t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")

# Create DESeq2 dataset
countData <- dataRow
metaData <- pitpat_metadata
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~condition, tidy = TRUE)
dds <- DESeq(dds)

# Normalize data
dds_norm <- estimateSizeFactors(dds)
counts_norm_dds_1 <- counts(dds_norm, normalized = TRUE)

# DESeq2 analysis

# Filter rows based on row sums (keeping rows with counts >= 10)
keep <- rowSums(counts(dds)) >= 10

# Subset the DESeqDataSet to include only the selected rows
dds_clean <- dds[keep,]

# Run DESeq on the cleaned DESeqDataSet
dds.1 <- DESeq(dds_clean)

# Estimate size factors
dds.1 <- estimateSizeFactors(dds.1)

# Obtain normalized counts using the DESeqDataSet
counts_norm_dds_1 <- counts(dds.1, normalized = TRUE)

# Perform differential gene expression analysis using DESeq2
resNFPA <- results(dds.1, contrast = c("condition", "non-functioning", "growth-hormone"), alpha = 0.05)

# Apply lfcShrink to improve log2 fold change estimates
resNFPA <- lfcShrink(dds.1, contrast = c("condition", "non-functioning", "growth-hormone"), res = resNFPA, type = "ashr")

# Filter DE results based on log2 fold change threshold
resNFPA.lfc0.32 <- results(dds.1, contrast = c("condition", "non-functioning", "growth-hormone"), alpha = 0.05, lfcThreshold = 0.32)

# Annotation with Annotable

# Convert DESeq2 results to a data frame
ensemble.ids1 <- as.data.frame(resNFPA.lfc0.32)

# Add a new column 'ensgene' containing row names
ensemble.ids1$ensgene <- rownames(ensemble.ids1)

# Reset row names to NULL
rownames(ensemble.ids1) <- NULL

# Arrange results by p adjusted value and record all genes

# Annotate DE results with additional information from the 'grch38' data frame
Annotated.res.NFPAlfc0.32 <- ensemble.ids1 %>%
          arrange(padj) %>%                    # Arrange the results by adjusted p-value
          inner_join(grch38, by = "ensgene") %>%  # Inner join with 'grch38' data frame based on 'ensgene'
          dplyr::select(                       # Select specific columns
                    ensgene, entrez, symbol, padj, log2FoldChange, chr, biotype, description, description
          )

# Deduplicate based on 'ensgene' column
Anno.dedup.res.NFPAlfc0.32 <- Annotated.res.NFPAlfc0.32[!duplicated(Annotated.res.NFPAlfc0.32[c('ensgene')]), ]

# Deduplicate based on 'symbol' column within the deduplicated 'ensgene' results
Anno.dedup.res.NFPAlfc0.32 <- Anno.dedup.res.NFPAlfc0.32[!duplicated(Anno.dedup.res.NFPAlfc0.32[c('symbol')]), ]

# Remove white spaces from the 'symbol' column
Anno.dedup.res.NFPAlfc0.32$symbol <- gsub("\\s", "", Anno.dedup.res.NFPAlfc0.32$symbol)


# Set row names of Anno.dedup.res.NFPAlfc0.32 to 'ensgene'
row.names(Anno.dedup.res.NFPAlfc0.32) <- Anno.dedup.res.NFPAlfc0.32$ensgene

# Select the top 30 genes based on adjusted p-value
top30NFPAvGH. <- Anno.dedup.res.NFPAlfc0.32 %>%
          arrange(padj) %>%
          head(30)

# Extract normalized counts for the top 30 genes
top30_norm_counts_NFPA <- counts_norm_dds_1[rownames(top30NFPAvGH.), ]

# Scale the normalized counts
z.top30_norm_counts_NFPA <- t(apply(top30_norm_counts_NFPA, 1, scale))

# Set column names of z.top30_norm_counts_NFPA to sample names
colnames(z.top30_norm_counts_NFPA) <- rownames(metaData)

# Create a heatmap for 30 genes
Heatmap(
          z.top30_norm_counts_NFPA,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          column_labels = colnames(z.top30_norm_counts_NFPA),
          row_labels = top30NFPAvGH.$symbol,
          row_dend_reorder = TRUE,  # Adjusts the order of row labels
          column_dend_reorder = TRUE,  # Adjusts the order of column labels
          column_names_gp = gpar(fontsize = 8, fontface = "bold", col = "blue"),  # Adjusts appearance of column names
          row_names_gp = gpar(fontsize = 7, fontface = "bold", col = "black")  # Adjusts appearance of row names
)
################################################################################

## Create a heatmap for 60 genes
# Scale the normalized counts
z.top60_norm_counts_NFPA <- t(apply(top60_norm_counts_NFPA, 1, scale))

# Set column names of z.top60_norm_counts_NFPA to sample names
colnames(z.top60_norm_counts_NFPA) <- rownames(metaData)

# Create a heatmap with clustering for rows and columns
Heatmap(
          z.top60_norm_counts_NFPA,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          column_labels = colnames(z.top60_norm_counts_NFPA),
          row_labels = top60NFPAvGH.$symbol,
          row_dend_reorder = TRUE,  # Adjusts the order of row labels
          column_dend_reorder = TRUE,  # Adjusts the order of column labels
          column_names_gp = gpar(fontsize = 8, fontface = "bold", col = "blue"),  # Adjusts appearance of column names
          row_names_gp = gpar(fontsize = 6
                              , fontface = "bold", col = "black")  # Adjusts appearance of row names
)
