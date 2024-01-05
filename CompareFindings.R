# Install and load necessary packages
#install.packages("dplyr")  # Install 'dplyr' package if not already installed
library(dplyr)  # Load 'dplyr' package for data manipulation
#install.packages("readxl")  # Install 'readxl' package if not already installed
library(readxl)  # Load 'readxl' package for reading Excel files

# Common data with Peculis doi:10.3390/cancers13061395 
# Extracting data from our results and comparing against Peculis Transcriptomic study to identify similarities and differences
results_NvG <- read.csv("C:/Users/User/Documents/experiment/all_results_preserved/Repete-theSame_15.10.23/Final_res.anno.dedup.NFPA.pvalue0.05.lfc0.32.csv")
colnames(results_NvG)
results_NvG <- data.frame(results_NvG)

results_NvG_padj005 <- results_NvG[results_NvG$padj < 0.05,]

Peculis_exe <- readxl::read_excel("Peculis_cancers-13-01395-s001.xlsx")
colnames(Peculis_exe)
head(Peculis_exe)

# Merging our results with Peculis results
merged_ResvPeculispadj0.05 <- inner_join(Peculis_exe, results_NvG_padj005, by = "symbol")
# Creating a file with common genes and similar gene expression.
write.csv(merged_ResvPeculispadj0.05, "merged_ResvPeculispadj0.05.csv")

# Common data with Kim invasive vs Non-invasive PA LFC > 2; DFR 0.01 less top 20 genes
# Comparing data from Kim et al 2019 doi: 10.3803/EnM.2019.34.3.314

Kim_exe <- readxl::read_excel("Kim.Invasive.vs.Non_invasiveLFCmore2FDRless.01.xlsx")
colnames(Kim_exe)
merged_ResvKimpadj0.05 <- inner_join(Kim_exe, results_NvG_padj005, by = "symbol")

colnames(merged_ResvKimpadj0.05)

KimLFC2.FDR0.01v_ourResults_NvG.LFC1.25padj0.05 <- merged_ResvKimpadj0.05[, c("symbol", "ensgene", "entrez", "padj", "log2FoldChange", "KimLog2FC", "KimFDR", "chr", "biotype", "description")]
# Creating a file with common genes and similar gene expression.
write.csv(KimLFC2.FDR0.01v_ourResults_NvG.LFC1.25padj0.05, "KimLFC2.FDR0.01v_ourResults_NvG.LFC1.25padj0.05.csv")

# Common part with Salomon # GH vs non-hormone producing PitNets 
Salomon_exe <- read.csv("Salomon2018.csv")
colnames(Salomon_exe)
# Creating a file with common genes and similar gene expression.
MergeSalomon.vs.OurStudy <- inner_join(Salomon_exe, results_NvG_padj005, by = "symbol")
colnames(MergeSalomon.vs.OurStudy)
# Choosing appropriate columns
Salomon.LFC.more1.GHvNonhormonal_ourStudyGHvNF.LFC.1.25 <- MergeSalomon.vs.OurStudy[, c("symbol", "ensgene", "entrez", "SalomonLFC", "SalomondeltaBeta", "Salomonfeature", "padj", "log2FoldChange", "chr", "biotype", "description")]
write.csv(Salomon.LFC.more1.GHvNonhormonal_ourStudyGHvNF.LFC.1.25, "Salomon.LFC.more1.GHvNonhormonal_ourStudyGHvNF.LFC.1.25.csv")
###############################################################

# Common data with Ronchi et al 2015
Ronchi <- read.csv("Ronchi.csv")
colnames(Ronchi)

Ronchi.vs.OurStudy <- inner_join(Ronchi, results_NvG_padj005, by = "symbol")
# Creating a file with common genes and similar gene expression.
write.csv(Ronchi.vs.OurStudy, "Ronchi.vs.OurStudy.csv")
################################################################

# Falch slow and fast-growing PitNets 
Falch <- readxl::read_excel("Falch.xlsx")
MergeFalch.vs.OurStudy <- inner_join(Falch, results_NvG_padj005, by = "symbol")
colnames(MergeFalch.vs.OurStudy)
# Creating a file with common genes and similar gene expression.
Falch40Top_ourResults_NvG.LFC1.25padj0.05 <- MergeFalch.vs.OurStudy[, c("symbol", "Falchpadjusted", "FalchLC", "EMT", "Cancer", "ensgene", "entrez", "padj", "log2FoldChange", "chr", "biotype", "description")]
write.csv(Falch40Top_ourResults_NvG.LFC1.25padj0.05, "Falch40Top_ourResults_NvG.LFC1.25padj0.05.csv")
