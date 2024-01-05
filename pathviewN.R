# Install and load necessary packages
install.packages("dplyr")
library(dplyr)
install.packages("pathview")
library(pathview)

# Read and process DESeq2 results NFPA vs GH (log2FC .032)
results_NvG <- read.csv("C:/Users/User/Documents/experiment/all_results_preserved/Repete-theSame_15.10.23/Final_res.anno.dedup.NFPA.pvalue0.05.lfc0.32.csv")

# Display column names of the loaded data
colnames(results_NvG)

# Convert results_NvG to a data frame (may not be necessary)
results_NvG <- data.frame(results_NvG)

# Filter genes with adjusted p-value < 0.05
results_NvG_padj005 <- results_NvG [results_NvG$padj < 0.05,]

# Display dimensions and column names of the filtered data
dim(results_NvG_padj005)
colnames(results_NvG_padj005) #6616 9 

# Create a named vector for log2FoldChange values
f <- results_NvG_padj005$log2FoldChange

# Display the row numbers of the named vector
names(f) <- results_NvG_padj005$entrez
row_number(f)

# Run pathview for the Pre-Initiation Complex (PreIC) formation pathway (KEGG ID: 04110)
PreIC_formation <- pathview(gene.data=f, pathway.id="04110", species="hsa")

# Display the current working directory
getwd()



