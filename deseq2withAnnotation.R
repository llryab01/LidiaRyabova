
# DESeq2 analysis and annotation for Lidia Ryabova's MSc Dissertation 
# Title: GH and NFPA PitPat analysis, Pathways in Pituitary Adenoma
# Credits to: Love et al., 2014 (Reference to the paper o& code I am using)
              #Mistry et al., 2021 (Reference to the paper o& code I am using)

# Get the current working directory
getwd()

# Remove all existing objects in the workspace
rm(list=ls())

# Install and load necessary libraries

# DESeq2 for differential expression analysis
#if (!requireNamespace("DESeq2", quietly = TRUE))
         # install.packages("DESeq2")
library(DESeq2)

# htmltools for creating HTML content
#if (!requireNamespace("htmltools", quietly = TRUE))
          #install.packages("htmltools")
library(htmltools)

# ggplot2 for creating visualizations
#if (!requireNamespace("ggplot2", quietly = TRUE))
          #install.packages("ggplot2")
library(ggplot2)

# dplyr for data manipulation
#if (!requireNamespace("dplyr", quietly = TRUE))
          #install.packages("dplyr")
library(dplyr)

# RColorBrewer for color palettes
#if (!requireNamespace("RColorBrewer", quietly = TRUE))
          #install.packages("RColorBrewer")
library(RColorBrewer)

# pheatmap for creating heatmaps
#if (!requireNamespace("pheatmap", quietly = TRUE))
          #install.packages("pheatmap")
library(pheatmap)

# tidyverse for data manipulation and visualization
#if (!requireNamespace("tidyverse", quietly = TRUE))
          #install.packages("tidyverse")
library(tidyverse)

# ashr for multiple testing correction
#if (!requireNamespace("ashr", quietly = TRUE))
          #install.packages("ashr")
library(ashr)

# BiocManager for managing Bioconductor packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
          #install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

# org.Hs.eg.db for mapping between gene IDs and symbols
library(org.Hs.eg.db)


# Read featureCounts results file to 
f <- read.table("D:/shared_D_folder/PitPat_analysis/featureCountsP_file/featureCounts_resultsP.txt", header = TRUE)

# Display column names and preview the data
colnames(f)
head(f)

# Subsetting the feature counts matrix by removing columns 2-6
# Since only one control is available, and DESeq2 requires at least three controls,
# we remove the control sample that does not serve a purpose.
# Instead, two groups are contrasted directly.
dataRow <- subset(f[c(-2:-6), ], select = c(1, 8, 9, 10, 11, 12, 13, 14, 15))

# Display the first few rows of the subsetted data
head(dataRow)

# Rename columns in the 'dataRow' matrix
colnames(dataRow) <- c(" ", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")

# Save the 'dataRow' matrix to a CSV file
write.csv(dataRow, "dataRow.csv")

# Create a metadata data frame
condition <- c("growth-hormone", "growth-hormone", "non-functioning", "non-functioning", "non-functioning", "non-functioning", "growth-hormone", "non-functioning")
Sample.names <- c("t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")

pitpat_metadata <- data.frame(Sample.names, condition)

# Set row names in 'pitpat_metadata'
row.names(pitpat_metadata) <- c("t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")

# Print and display the 'pitpat_metadata' data frame
print(pitpat_metadata)

# Save the 'pitpat_metadata' data frame to a tab-delimited text file
write.table(pitpat_metadata, "metadataExtraCol.txt", col.names = TRUE, row.names = TRUE, sep = "\t")

# Set 'countData' to the 'dataRow' matrix
countData = dataRow

# Assign metadata to 'metaData'
metaData = pitpat_metadata

# Display the first few rows of the 'countData' matrix
head(countData)

# Display the first few rows of the 'metaData' data frame
head(metaData)

# Create a DESeqDataSet object from the count matrix and metadata
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~condition, tidy = TRUE)

# Display information about the 'dds' object
dds

# Run DESeq analysis on the 'dds' object
dds <- DESeq(dds)

# Estimate size factors for normalization
dds_norm <- estimateSizeFactors(dds)

# Obtain normalized counts from the DESeq object
counts_norm <- counts(dds_norm, normalized = TRUE)

# Display the dimensions of the normalized counts matrix
dim(counts_norm)  # 60659 rows, 8 columns

# Filter out transcripts with no/low expression (less than 10 counts)
keep <- rowSums(counts(dds)) >= 10

# Create a new DESeqDataSet by keeping only rows with sufficient counts
dds_clean <- dds[keep,]

# Run DESeq analysis on the cleaned DESeqDataSet
dds.1 <- DESeq(dds_clean)

# View the results of the DESeq analysis
View(dds.1)

# Estimate size factors for normalization on the updated DESeqDataSet
dds.1 <- estimateSizeFactors(dds.1)

# Obtain normalized counts from the updated DESeq object
counts_norm_dds_1 <- counts(dds.1, normalized = TRUE)

# Display the dimensions of the normalized counts matrix
dim(counts_norm_dds_1)  # 28639 rows, 8 columns

# Display the first few rows of the normalized counts matrix
head(counts_norm_dds_1)

# Save normalized counts to a CSV file
write.csv(counts_norm_dds_1, "countsNorm_DDS_1.csv")

# Uncomment the line below if you want to read the CSV file back into R
# normalized_counts_csv <- read.csv("countsNorm_DDS_1.csv")

# Save normalized counts to a tab-delimited text file
write.table(counts_norm_dds_1, "CountsNorm.txt", col.names = TRUE, row.names = TRUE, sep = "\t")


# Perform VST transformation on the DESeqDataSet
vsdata <- vst(dds.1, blind = TRUE)  # 29824 features after transformation

# Plot PCA using the transformed data
plotPCA(vsdata, intgroup = c('Sample.names', 'condition'))

# Save the PCA plot as a PNG file
dev.copy(png, 'PCA_Sample_condition.png')
dev.off()


# Perform variance stabilizing transformation on the DESeqDataSet with blind set to TRUE
vsdata <- vst(dds.1, blind = TRUE)

# Extract the variance-stabilized matrix from the VST object
vsdata_mat <- assay(vsdata)

# Calculate the Pearson correlation matrix of the variance-stabilized matrix
vsdata_mat_cor <- cor(vsdata_mat)

# Create a separate heatmap without annotations (if needed)
plotHeatmap <- pheatmap(vsdata_mat_cor)

# Save the heatmap with annotations as a PNG file
dev.copy(png, "CorPheatmapWithSampleNames.png")
dev.off()

# Create a heatmap using pheatmap with annotations from metaData
pheatmap(vsdata_mat_cor, annotation = dplyr::select(metaData, c("Sample.names", "condition")))
# Save the heatmap with annotations as a PNG file
dev.copy(png, "CorrelationPheatma.SampleNames.Conditions.png")
dev.off()


# Plot gene dispersions using plotDispEsts function
plotDispEsts(dds.1)

# Save the gene dispersions plot as a PNG file
dev.copy(png, 'plotDispEstsdds1.png')
dev.off()
################################################################################

# Perform differential gene expression analysis with DESeq2 using contrasts
# Contrast: Non-Functioning vs. Growth Hormone condition
resNFPA <- results(dds.1, contrast = c("condition", "non-functioning", "growth-hormone"), alpha = 0.05)

# Display a summary of the results
summary(resNFPA)

# Display the first few rows of the results
head(resNFPA)

# Display the dimensions of the results
dim(resNFPA)  # 28639 rows, 6 columns

# Plot a MA plot for the results
plotMA(resNFPA, ylim = c(-8, 8))

# Save the MA plot as a PNG file
dev.copy(png, "PlotMA_resNFPA.data.png")
dev.off()

# Apply shrinkage to log2 fold changes using ashr method
resNFPA <- lfcShrink(dds.1, contrast = c("condition", "non-functioning", "growth-hormone"), res = resNFPA, type = "ashr")

# Display the dimensions of the shrunken results
dim(resNFPA)

# Display the first few rows of the shrunken results
head(resNFPA)

# Plot a MA plot for the shrunken results
plotMA(resNFPA, ylim = c(-10, 10))

# Save the shrunken MA plot as a PNG file
dev.copy(png, "PlotMA_lfcShrinkNFPA.png")
dev.off()
write.csv(resNFPA, "results_NFPAvisGH_pvalue0.05.csv")

##############################################################################
# Adding a fold change threshold:
# With large significant gene lists, extracting meaningful biological relevance can be challenging.
# To increase stringency, a fold change threshold can be added.

# Set the fold change threshold to 0.32, which corresponds to log2(1.25)
resNFPA.lfc0.32 <- results(dds.1, contrast = c("condition", "non-functioning", "growth-hormone"),
                           alpha = 0.05,
                           lfcThreshold = 0.32)

# Display a summary of the results with the fold change threshold
summary(resNFPA.lfc0.32)

# Display the dimensions of the results with the fold change threshold
dim(resNFPA.lfc0.32)  # 28639 rows, 6 columns

# Display the first few rows of the results with the fold change threshold
head(resNFPA.lfc0.32)

# Write the results with the fold change threshold to a CSV file
write.csv(resNFPA.lfc0.32, "results_NFPA.vs.GH.pvalue0.005.lfc0.32.csv")

# Annotation with annotable (see folder Annotable)
## to make : list of up and down regulated genes 9 subset )
## extract Novel genes 


# #############################################################################
#Visualization of the results with a Volcano Plot with resNFPA.lfc0.32

# Generate a logical column indicating significance (padj < 0.05)
resNFPA_all <- data.frame(resNFPA.lfc0.32) %>% mutate(threshold = padj < 0.05)

# Create the volcano plot for all values with padj < 0.05

#Filter out rows with missing values in the color aesthetic
resNFPA_all <- resNFPA_all[complete.cases(resNFPA_all$threshold), ]


# Create the volcano plot for all values with padj < 0.05
  ggplot2::ggplot(resNFPA_all) + 
  ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))
  # Save the Volcano Plot as a PNG file
  
dev.copy(png,"volcanoPlotResNFPAlfc0.32.padj0.05SigGenes.png")
dev.off()
###############################################################################
#Building pheatmap with sig genes and resNFPA.lfc0.32
# Load the dplyr library for data manipulation 
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Arrange the rows of the data frame 'resNFPA.lfc0.32' based on adjusted p-values
resNFPA <- resNFPA.lfc0.32[order(resNFPA.lfc0.32$padj),]

# Display the first few rows of the arranged data frame
head(resNFPA)

# Display a summary of the arranged data frame
summary(resNFPA)

# Convert 'resNFPA' to a data frame
resNFPA_all <- data.frame(resNFPA)

# Subset the data frame to include only rows with padj < 0.05
resNFPA_Sig <- subset(resNFPA_all, padj < 0.05, na.rm = TRUE)

# Display a summary of the subsetted data frame
summary(resNFPA_Sig)

#Subset normalized counts to significant genes
sig_norm_counts_NFPA <- counts_norm_dds_1[rownames(resNFPA_Sig), ]
head(sig_norm_counts_NFPA)
 #Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

pheatmap(sig_norm_counts_NFPA,
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation = dplyr::select(pitpat_metadata, c("Sample.names", "condition")), 
         scale = "row")

dev.copy(png,"pheatmapSigGenesNFPAvGH.png")
dev.off()
################################################################################
#Annotation with ANNOTABLE NFPA vs GH results 

#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
library(annotables)
library(dplyr)

# Read the results file into a data frame
ensemble.ids1 <- as.data.frame(resNFPA.lfc0.32)
head(ensemble.ids1)

# Assign row names to a new column named "ensgene"
ensemble.ids1$ensgene <- rownames(ensemble.ids1)

# Display the first few rows 
head(ensemble.ids1)
# Reset row names to NULL
rownames(ensemble.ids1) <- NULL

# Print the first few rows of ensemble.ids1
head(ensemble.ids1)
colnames(ensemble.ids1)
dim(ensemble.ids1) #28639 7

# Arrange results by p adjusted value and record all genes
ensemble.ids1 %>% 
          arrange(padj) %>% 
          inner_join(grch38, by= "ensgene" ) %>% 
          dplyr:::select(ensgene,entrez,symbol,padj,log2FoldChange, chr, biotype,description, description) -> Annotated.res.NFPAlfc0.32 

# Display the first few rows and dimensions of the annotated data frame
colnames(Annotated.res.NFPAlfc0.32)
head(Annotated.res.NFPAlfc0.32)
dim(Annotated.res.NFPAlfc0.32) #28762     8


# Deduplicate based on 'ensgene'
Anno.dedup.res.NFPAlfc0.32 <- Annotated.res.NFPAlfc0.32[!duplicated(Annotated.res.NFPAlfc0.32[c('ensgene')]), ]
# Display the dimensions of the deduplicated data fra
dim(Anno.dedup.res.NFPAlfc0.32) # 8 28569 
head(Anno.dedup.res.NFPAlfc0.32)

# Save the deduplicated result to a CSV file
write.csv(Anno.dedup.res.NFPAlfc0.32, "Final_res.anno.dedup.NFPA.pvalue0.05.lfc0.32.csv")

# Subset upregulated and downregulated genes
upreg_Anno.dedup.res.NFPAlfc0.32 <-  subset(Anno.dedup.res.NFPAlfc0.32, Anno.dedup.res.NFPAlfc0.32$log2FoldChange > 0 )
downreg_Anno.dedup.res.NFPAlfc0.32 <- subset(Anno.dedup.res.NFPAlfc0.32, Anno.dedup.res.NFPAlfc0.32$log2FoldChange < 0 )

# Display the first few rows of upregulated genes
head(upreg_Anno.dedup.res.NFPAlfc0.32) 

# Save the results to CSV files
write.csv(upreg_Anno.dedup.res.NFPAlfc0.32, "Final_res.upreg_Anno.dedup.NFPAlfc0.32.csv")
write.csv(downreg_Anno.dedup.res.NFPAlfc0.32, "Final_res.downreg_Anno.dedup.NFPAlfc0.3.csv")


################################################################################
# Principal Investigator's Requested Analysis
# The following code performs differential gene expression analysis and visualizations
# as requested by the Principal Investigator for the GH vs. NFPA comparison.
# The generated plots include MA plots, shrunken MA plots, and volcano plots.
# Results are also compared with NFPA vs. GH for discrepancies.

# DESeq2 ANALYSIS WITH CONTRASTS FOR GH VS. NFPA
# Visualized by GH results with contrast, MA plot, shrinked MA plot, volcano, and summary comparisons

# Get the current working directory
getwd()

# Perform differential gene expression analysis for GH vs. NFPA
resGH <- results(dds.1, contrast = c("condition", "growth-hormone", "non-functioning"), alpha = 0.05)

# Display a summary of the GH vs. NFPA results
summary(resGH)

# Display the first few rows of the results
head(resGH)

# Display the dimensions of the results
dim(resGH)  # 28639 rows, 6 columns

# Plot an MA plot for GH vs. NFPA results
plotMA(resGH, ylim = c(-8, 8))
dev.copy(png, "PlotMA_resGHvsNFPA.png")
dev.off()

# Apply shrinkage to log2 fold changes for GH vs. NFPA
resGH <- lfcShrink(dds.1, contrast = c("condition", "growth-hormone", "non-functioning"), res = resGH, type = "ashr")

# Display the dimensions of the shrunken results
dim(resGH)

# Display the first few rows of the shrunken results
head(resGH)

# Plot an MA plot for shrunken GH vs. NFPA results
plotMA(resGH, ylim = c(-10, 10))
dev.copy(png, "PlotMA_GHvsNFPA_lfcShrinkGH.png")
dev.off()

# Save the GH vs. NFPA results to a CSV file
write.csv(resGH, "results_visGHvsNFPA_pvalue0.05.csv")

# Adding a fold change threshold for GH vs. NFPA results
resGH.lfc0.32 <- results(dds.1, contrast = c("condition", "growth-hormone", "non-functioning"), 
                         alpha = 0.05, lfcThreshold = 0.32)

# Display a summary of the results with the fold change threshold
summary(resGH.lfc0.32)

# Display the dimensions of the results with the fold change threshold
dim(resGH.lfc0.32)  # 28639 rows, 6 columns

# Display the first few rows of the results with the fold change threshold
head(resGH.lfc0.32)

# Save the results with the fold change threshold to a CSV file
write.csv(resGHlfc0.32, "results_GHvsNFPApvalue0.005.lfc0.32.csv")

# Comparison of summary results between GH vs. NFPA and NFPA vs. GH
# Display summary results for GH vs. NFPA
summary(resGH.lfc0.32)
# out of 28639 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0.32 (up)    : 1116, 3.9%
# LFC < -0.32 (down) : 1200, 4.2%
# outliers [1]       : 439, 1.5%
# low counts [2]     : 3887, 14%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Display summary results for NFPA vs. GH for comparison
 summary(resNFPA.lfc0.32)
# out of 28639 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0.32 (up)    : 1200, 4.2%
# LFC < -0.32 (down) : 1116, 3.9%
# outliers [1]       : 439, 1.5%
# low counts [2]     : 3887, 14%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


################################################################################

library(dplyr)
resGH <- resGH.lfc0.32[order(resGH.lfc0.32$padj),]
head(resGH)
summary(resGH)
resGH_all <- data.frame(resGH)
resGH_Sig <- subset(resGH_all, padj < 0.05, na.rm=TRUE)
summary(resGH_Sig)

# Create the volcano plot for all values with padj < 0.05 and lfc0.32
resGH_all <- data.frame(resGH) %>% mutate(threshold = padj < 0.05)

# Filter out rows with missing values in padj or log2FoldChange
resGH_all <- resGH_all[complete.cases(resGH_all$padj, resGH_all$log2FoldChange), ]
ggplot2::ggplot(resGH_all) + 
  ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))

dev.copy(png,"volcanoPlotResGHsigGenelfc0.32.Padjless0.05.png")
dev.off()




#########################################################################
library(pheatmap)
library(tidyverse)
library(dbplyr)
library(RColorBrewer)
#Subset normalized counts to significant genes
sig_norm_counts_GH <- counts_norm_dds_1[rownames(resGH_Sig), ]
head(sig_norm_counts_GH)
dim(sig_norm_counts_GH)

# Plot heatmap
pheatmap(sig_norm_counts_GH,
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation = dplyr:::select(pitpat_metadata, c("Sample.names", "condition")), 
         scale = "row")

dev.copy(png,"pheatmapSigGenes_GHvsNFPALFC0.32.PADJ0.05.png")
dev.off()

################################################################################

#RNA expression between GH and NFPA for our top 6 genes.

#reset par
par(mfrow=c(2,3))

plotCounts(dds.1, gene="ENSG00000162670", intgroup="condition")
plotCounts(dds.1, gene="ENSG00000179348", intgroup="condition")
plotCounts(dds.1, gene="ENSG00000188921", intgroup="condition")
plotCounts(dds.1, gene="ENSG00000151062", intgroup="condition")
plotCounts(dds.1, gene="ENSG00000092929", intgroup="condition")
plotCounts(dds.1, gene="ENSG00000117791", intgroup="condition")


###############################################################################################################
getwd()
# Create a data frame with gene symbols and names
gene_data <- data.frame(
          gene_symbol = c("ENSG00000162670", "ENSG00000179348", "ENSG00000188921", 
                          "ENSG00000151062", "ENSG00000092929", "ENSG00000117791"),
          gene_name = c("BRINP3", "GATA2", "HACD4", "CACNA2D4", "UNC13D", "MTARC2")
)

# Set up the layout for the plots
par(mfrow = c(2, 3))

# Loop through each gene and generate plotCounts
for (i in seq_along(gene_data$gene_symbol)) {
          current_gene <- gene_data[i, ]
          plotCounts(dds.1, gene = current_gene$gene_symbol, intgroup = "condition",
                     main = current_gene$gene_name)
          
          # Add legend with gene symbol and name
          legend("topright", legend = paste(current_gene$gene_symbol, current_gene$gene_name),
                 col = "black", pch = 1, cex = 0.8)  # Adjust the font size (0.8 is an example size)
}


################################################################################
#Annotation with ANNOTABLE GH vs NFPA results 

#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
library(annotables)
library(dplyr)

# Read the results file into a data frame
ensemble.ids2 <- as.data.frame(resGH.lfc0.32)
head(ensemble.ids2)

# Assign row names to a new column named "ensgene"
ensemble.ids2$ensgene <- rownames(ensemble.ids2)

# Display the first few rows 
head(ensemble.ids2)
# Reset row names to NULL
rownames(ensemble.ids2) <- NULL

# Print the first few rows of ensemble.ids1
head(ensemble.ids2)
colnames(ensemble.ids2)
dim(ensemble.ids2) #28639 7

# Arrange results by p adjusted value and record all genes
ensemble.ids2 %>% 
          arrange(padj) %>% 
          inner_join(grch38, by= "ensgene" ) %>% 
          dplyr:::select(ensgene,entrez,symbol,padj,log2FoldChange, chr, biotype,description, description) -> Annotated.res.GHlfc0.32 

# Display the first few rows and dimensions of the annotated data frame
colnames(Annotated.res.GHlfc0.32)
head(Annotated.res.GHlfc0.32)
dim(Annotated.res.GHlfc0.32) #28762     8


# Deduplicate based on 'ensgene'
Anno.dedup.res.GHlfc0.32 <- Annotated.res.GHlfc0.32[!duplicated(Annotated.res.GHlfc0.32[c('ensgene')]), ]
# Display the dimensions of the deduplicated data fra
dim(Anno.dedup.res.GHlfc0.32) # 8 28569 

# Save the deduplicated result to a CSV file
write.csv(Anno.dedup.res.GHlfc0.32, "Final_res.anno.dedup.GHvsNFPA.pvalue0.05.lfc0.32.csv")
head(Anno.dedup.res.GHlfc0.32)

# Subset upregulated and downregulated genes
upreg_Anno.dedup.res.GHlfc0.32 <-  subset(Anno.dedup.res.GHlfc0.32, Anno.dedup.res.GHlfc0.32$log2FoldChange > 0 )
downreg_Anno.dedup.res.GHlfc0.32 <- subset(Anno.dedup.res.GHlfc0.32, Anno.dedup.res.GHlfc0.32$log2FoldChange < 0 )

# Display the first few rows of upregulated genes
head(upreg_Anno.dedup.res.GHlfc0.32) 

# Save the results to CSV files
write.csv(upreg_Anno.dedup.res.GHlfc0.32, "Final_res.upreg_Anno.dedup.GH.lfc0.32.csv")
write.csv(downreg_Anno.dedup.res.GHlfc0.32, "Final_res.downreg_Anno.dedup.GH.lfc0.3.csv")

