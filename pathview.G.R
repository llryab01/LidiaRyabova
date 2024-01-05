# Install and load necessary packages
install.packages("dplyr")
library(dplyr)
install.packages("pathview")
library(pathview)


# Read the GH vs NF results file
results_GvN <- read.csv("C:/Users/User/Documents/experiment/all_results_preserved/Repete-theSame_15.10.23/GHvsNFPArun/sorted_Final_res.anno.dedup.GHvsNFPA.pvalue0.05.lfc0.32.csv")

# Display column names of the loaded data
colnames(results_GvN)

# Convert results_GvN to a data frame
results_GvN <- data.frame(results_GvN)
colnames(results_GvN)

# Display the first few rows of the data
head(results_GvN)

# Filter genes with adjusted p-value < 0.05
results_GvNpadj005 <- results_GvN[results_GvN$padj < 0.05,]

# Display column names of the filtered data
colnames(results_GvNpadj005)

# Display dimensions of the filtered data
dim(results_GvNpadj005)  # 6616 rows, 9 columns

# Create a named vector for log2FoldChange values
g <- results_GvNpadj005$log2FoldChange
names(g) <- results_GvNpadj005$entrez

# Display the row numbers of the named vector
row_number(g)

# Run pathview for Type I Interferon Signaling Pathway (KEGG ID: 04620)
Type1_int0462 <- pathview(gene.data = g, pathway.id = "04620", species = "hsa")

# Run pathview for Type II Interferon to JAK-STAT Signaling Pathway (KEGG ID: 04630)
Type2_int_to_JS <- pathview(gene.data = g, pathway.id = "04630", species = "hsa")

# Run pathview for Interferon RIPK1 3 Signaling Pathway (KEGG ID: 04217)
int_to_RIPK1_3 <- pathview(gene.data = g, pathway.id = "04217", species = "hsa")

