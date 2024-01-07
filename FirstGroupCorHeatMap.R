
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



rm(list=ls())
# Read featureCounts results file to 
f <- read.table("D:/shared_D_folder/PitPat_analysis/featureCountsP_file/featureCounts_resultsP.txt", header = TRUE)


colnames(f)
head(f)

dataRow = subset(f[c(-2:-6),], select = c(1, 7,8,9,10,11,12,13,14,15))
head(dataRow)
dim(dataRow)

colnames(dataRow) <- c('', 'ctrl','t2','t3','t4','t5','t6','t7','t8','t9')
# remove all rows that contains all '0' in all columns 
colnames(dataRow)
rownames(dataRow)#59659 entries
dim(dataRow) # 60659 by 10
getwd()

head(dataRow)


genotype = c('wt','unknown','unknown', 'unknown','unknown','unknown','unknown','unknown','unknown')

#create condition vector 
condition = c('normal','adenoma','adenoma','adenoma','adenoma','adenoma',
              'adenoma','adenoma','adenoma')
#create dataframe 
pitpat_metadata = data.frame(genotype,condition)
#View(pitpat_metadata)

#nrow(pitpat_metadata)
row.names(pitpat_metadata) = c( 'ctrl','t2','t3','t4','t5','t6','t7','t8','t9')

#write.table(pitpat_metadata, "metadata.txt", col.names = TRUE, row.names = TRUE, sep = "\t")

countData = dataRow
metaData = pitpat_metadata

dds = DESeqDataSetFromMatrix(countData=countData, 
                             colData=metaData, 
                             design=~condition, tidy = TRUE)
dds
dds = DESeq(dds)

res = results(dds, contrast=c("condition", "adenoma","normal")) # perform differential gene expression on row data
summary(res)
View(res)

keep = rowSums(counts(dds)) >= 10 # filter no/low expression transcripts
dds_clean = dds[keep,] # keep those rows that contains read counts = or> than 10
dds.1 = DESeq(dds_clean) # run DESeq function 


res.1 = results(dds.1, contrast=c("condition", "adenoma","normal")) # perform differential gene expression 
#between two classes -fold expression will be A-B


dds_norm <- estimateSizeFactors(dds.1)
# Normalize counts
counts_norm <-counts(dds_norm, normalized=TRUE)

# Variance Stabilizing Transformation (VST)
vsd_dds_norm <- vst(dds_norm, blind=TRUE)
vsd_mat_dds_norm <- assay(vsd_dds_norm)

# Calculate pairwise Pearson correlation coefficients
vsd_cor_dds_norm <- cor(vsd_mat_dds_norm)

# Heatmap without specific annotations
pheatmap(vsd_mat_dds_norm)
# Heatmap with specific annotations
pheatmap(vsd_mat_dds_norm, annotation = select( metadata, condition))

