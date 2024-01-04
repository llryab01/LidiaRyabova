#!/bin/bash

# Genome indexing command using STAR
STAR --runThreadN 6 \
     --runMode genomeGenerate \
     --genomeDir /DISK1/research/lryab01/data_pitpat/GenomeIndex \
     --genomeFastaFiles /DISK1/research/lryab01/RefGenomeSeq/GRCh38_latest_genomic.fna \
     --sjdbGTFfile /DISK1/research/lryab01/RefAnnotationSeq/GRCh38_latest_genomic.gff \
     --sjdbOverhang 99

# Explanation:
# --runThreadN: Number of threads or cores to use during genome generation.
# --runMode genomeGenerate: Specifies the mode for genome generation.
# --genomeDir: Directory path to store the generated genome index.
# --genomeFastaFiles: Path to the genomic FASTA file.
# --sjdbGTFfile: Path to the GTF file containing gene annotations.
# --sjdbOverhang: Length of spliced junction overhang. Should be set to read length - 1.

# Additional Notes:
# The above command creates a genome index using STAR for subsequent read alignment.
# Adjust the paths to the genomic FASTA file and GTF file based on your data.
# The sjdbOverhang parameter should be set based on the read length - 1.


# to make file executable use command : "chmod +x fileName.sh"
# to run file use command: "./your_script.sh"
