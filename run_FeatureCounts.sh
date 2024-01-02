#!/bin/bash

featureCounts -a "Homo_sapiens.GRCh38.104.gtf" -p -o "featureCounts_resultsP.txt" "ctrl.bam" "t2.bam" "t3.bam" "t4.bam" "t5.bam" "t6.bam" "t7.bam" "t8.bam" "t9.bam"

#Make command executable
#chmod +x run_featureCounts.sh 
#Run command
#./run_FeatureCounts.sh
