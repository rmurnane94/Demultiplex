#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp

/usr/bin/time -v ./demultiplex.py -r1 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz' \
    -r2 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz' \
    -r3 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz' \
    -r4 '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz' \
    -i '/projects/bgmp/shared/2017_sequencing/indexes.txt' \
    -qs 5 \
    -o 'output_files'