#!/bin/bash

# Trim Galore command for trimming paired-end FASTQ files
trim_galore --paired \
            --fastqc \
            --three_prime_clip_R1 2 \
            --three_prime_clip_R2 2 \
            t2_1.fq t2_2.fq

# Explanation:
# --paired: Specify that the input files are paired-end.
# --fastqc: Output trimmed files in FastQC format.
# --three_prime_clip_R1 2: Trim 2 base pairs from the 3' end of Read 1.
# --three_prime_clip_R2 2: Trim 2 base pairs from the 3' end of Read 2.
# t2_1.fq and t2_2.fq: Replace these with the actual names of your input files.

# Note: FastQC results showed very good quality.

# Additional comments:
# You may want to add any additional information or context specific to your dataset.
