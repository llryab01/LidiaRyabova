#!/bin/bash

# Set the path to the genome index directory
genome_index="/DISK1/research/lryab01/data_pitpat/GenomeIndex"

# Set the paths to input FASTQ files
read1="t2_1_val_1.fq"
read2="t2_2_val_2.fq"

# Set the output directory and sample name
output_dir="/DISK1/research/lryab01/data_pitpat/2trimedFiles/starNmax1OnePass/start2"
sample_name="Sample_t2"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Run STAR alignment
STAR --runThreadN 4 \
     --genomeDir "$genome_index" \
     --readFilesIn "$read1" "$read2" \
     --outFileNamePrefix "${output_dir}/${sample_name}_" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within KeepPairs \
     --outFilterMultimapNmax 1 \
     --quantMode GeneCounts
