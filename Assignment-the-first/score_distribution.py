#!/usr/bin/env python

import bioinfo
import matplotlib.pyplot as plt
import gzip

# r1 = './short_tests/short_R1.fq'
# r2 = './short_tests/short_R2.fq'
# r3 = './short_tests/short_R3.fq'
# r4 = './short_tests/short_R4.fq'

r1 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
r2 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
r3 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
r4 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'


#creating function to calculate mean distribution for each bp in the reads of each file
def mean_distribution(file: str, read_length: int) -> list:
    """Goes through all records of FASTQ file, converts phred quality scores to numbers and adds it to sum list for each base pair.
    Then calculates the mean for each base pair and returns a list with them all"""
   
   #initializing sum list to track quality score values from each bp in reads
    sum_list = []
    for x in range(read_length):
        sum_list.append(0)
    
    #going through the file and pulling out quality scores, converting them to numbers and adding to sum for each bp
    num_lines=0
    # with open (file, "r") as fastq:
    with gzip.open(file, 'rt') as fastq:
        for line in fastq:
            line=line.strip()
            if num_lines%4==3:
                for bp, letter in enumerate(line):
                    sum_list[bp]+=int(bioinfo.convert_phred(letter))
                    
            num_lines+=1

    #converting sums to means
    mean_list = [x/(num_lines/4) for x in sum_list]

    #returns the mean quality score for each bp of the reads
    return mean_list


#calculate the mean_distribution for each read file
#assign the distributions to variables
r1_mean_distribution = mean_distribution(r1,101)
r2_mean_distribution = mean_distribution(r2,8)
r3_mean_distribution = mean_distribution(r3,8)
r4_mean_distribution = mean_distribution(r4,101)

#function for generating plots
def generate_plot(r_file: str, distribution: list):
    """generates distribution plot showing mean quality score at each bp"""
    basepairs = range(len(distribution))

    with open(f'./distribution_plots/{r_file}_quality_distribution.tsv', 'w') as output_file:
    
        for bp in basepairs:
            print(f'{bp}\t{distribution[bp]}', file=output_file)

    plt.bar(x=basepairs, height=distribution, width=1, edgecolor= 'r')
    plt.title(f'{r_file} Quality Score Distribution')
    plt.xlabel('Base Pair')
    plt.ylabel('Mean Quality Score')
    
    plt.savefig(f'./distribution_plots/{r_file}_quality_distribution.png')
    plt.clf()

    

#generates plots
generate_plot('R1', r1_mean_distribution)
generate_plot('R2', r2_mean_distribution)
generate_plot('R3', r3_mean_distribution)
generate_plot('R4', r4_mean_distribution)

