#!/usr/bin/env python

#import dependencies
import bioinfo
import gzip
import argparse

#defining argparse variables
def get_args():
     parser = argparse.ArgumentParser(description="inputs and parameters")
     parser.add_argument("-r1", "--read1", help="read1 file", required=True)
     parser.add_argument("-r2", "--read2", help="read2 file", required=True)
     parser.add_argument("-r3", "--read3", help="read3 file", required=True)
     parser.add_argument("-r4", "--read4", help="read4 file", required=True) 
     parser.add_argument("-i", "--index_file", help="file with indexes", required=True)
     parser.add_argument("-qs", "--quality_score", help="quality score cutoff for indexes", required=True, type=int)
     parser.add_argument("-o", "--output_folder", help="folder for output read files", required=True)
     return parser.parse_args()
	
args = get_args() 



#making function for finding reverse compliments which will be used for finding matching indexes
def reverse_compliment(sequence: str) -> str:
    '''finds reverse compliment of sequence and returns as string'''
    base_matches ={'A': 'T', 'T':'A',
                   'C': 'G', 'G': 'C', 
                   'N': 'N'}
    reversed_sequence = sequence[::-1]
    rev_comp = ''
    for base in reversed_sequence:
         rev_comp += base_matches[base]

    return rev_comp


#make function to find average quality score of a sequence
def average_qc(sequence: str) -> float:
    '''finds mean quality score of bases in sequence and returns as float'''
    sum_qc=0
    seq_length = len(sequence)
    for base in sequence:
         sum_qc += bioinfo.convert_phred(base)

    return sum_qc/seq_length


#make function to write to output files
def write_output(header: str, seq: str, qs: str, file_name_format: str) -> str:
     '''takes in the the record components along with which file to write to and writes to appropriate files'''
     output_file_handles[file_name_format].write(f'{r1_header}\n')
     output_file_handles[file_name_format].write(f'{seq}\n')
     output_file_handles[file_name_format].write(f'+\n')
     output_file_handles[file_name_format].write(f'{qs}\n')


#read files
r1 = args.read1
r2 = args.read2
r3 = args.read3
r4 = args.read4


#making SET of known indexes from the index file
index_set = set()
# index_file = '/projects/bgmp/shared/2017_sequencing/indexes.txt'
index_file = args.index_file
with open(index_file, 'r') as open_file:
     for line in open_file:
          line=line.strip()
          ind_seq = line.split()[4]
          index_set.add(ind_seq)
index_set.remove('index')



#dictionaries for paired index counts. matched and mismatched
match_dict = {}
mismatch_dict = {}


#counters
matched_index_counter = 0
index_hopper_counter = 0
unknown_index_counter = 0

#quality threshold
quality_threshold = args.quality_score

#output folder
output_folder = args.output_folder


#create output file names first and then open them in a dictionary to make them accessible
output_file_names = []
for index in index_set:
    output_file_names.append(f'{index}_R1')
    output_file_names.append(f'{index}_R2')
output_file_names.extend(['nonmatch_R1', 'nonmatch_R2', 'unknown_R1', 'unknown_R2' ])
output_file_handles = {}
for name in output_file_names:
     output_file_handles[name] = open(f'./{output_folder}/{name}.fq','w')




#OPEN files to READ
#for actual files
with gzip.open(r1, 'rt') as open_r1, gzip.open(r2, 'rt') as open_r2,gzip.open(r3, 'rt') as open_r3, gzip.open(r4, 'rt') as open_r4:

#for test files
# with open(r1, 'r') as open_r1, open(r2, 'r') as open_r2, open(r3, 'r') as open_r3, open(r4, 'r') as open_r4:
    i=0
    current_record_r1 = []
    current_record_r2 = []
    current_record_r3 = []
    current_record_r4 = []

#go through each read file and compiles matching records across each file as it goes through. 
    while True:
            line_r1 = open_r1.readline()
            line_r2 = open_r2.readline()
            line_r3 = open_r3.readline()
            line_r4 = open_r4.readline()
            if line_r1 == '':
                break

            line_r1 = line_r1.strip()
            line_r2 = line_r2.strip()
            line_r3 = line_r3.strip()
            line_r4 = line_r4.strip()

            
            current_record_r1.append(line_r1)
            current_record_r2.append(line_r2)
            current_record_r3.append(line_r3)
            current_record_r4.append(line_r4)

            #when you get to a full record
            if i%4 == 3:
                # get both indexes from the current records
                index_r2 = current_record_r2[1]
                index_r3 = current_record_r3[1]

                #reverse compliment of r3 index
                index_r3_rev_comp = reverse_compliment(index_r3) 


                #get quality scores
                quality_r1 = current_record_r1[3]
                quality_r2 = current_record_r2[3]
                quality_r3 = current_record_r3[3]
                quality_r4 = current_record_r4[3]


                #get headers for R1 and R4
                r1_header = current_record_r1[0]
                r4_header = current_record_r4[0]

                #create modified headers
                header_modifier = f' {index_r2}-{index_r3_rev_comp}'
                r1_header += header_modifier
                r4_header += header_modifier    

                #get sequences for r1 and r4
                r1_seq =  current_record_r1[1]
                r4_seq =  current_record_r4[1]  
            

            
            
                #check to see if the indexes are unknown/below threshold
                if index_r2 not in index_set or index_r3_rev_comp not in index_set or average_qc(quality_r2) < quality_threshold or average_qc(quality_r3) < quality_threshold:
                    write_output(r1_header, r1_seq, quality_r1, 'unknown_R1')
                    write_output(r4_header, r4_seq, quality_r4, 'unknown_R2')
                    unknown_index_counter+=1
                     
                 
                #check to see if indexes are matching
                elif index_r2 == reverse_compliment(index_r3):  
                    write_output(r1_header, r1_seq, quality_r1, f'{index_r2}_R1')
                    write_output(r4_header, r4_seq, quality_r4, f'{index_r2}_R2')
                    matched_index_counter+=1

                    if match_dict.get(header_modifier):
                         match_dict[header_modifier] += 1
                    else:
                         match_dict[header_modifier] = 1

                           
                #else they are nonmatching (index hopped)
                else:     
                    write_output(r1_header, r1_seq, quality_r1, 'nonmatch_R1')
                    write_output(r4_header, r4_seq, quality_r4, 'nonmatch_R2')
                    index_hopper_counter+=1

                    if mismatch_dict.get(header_modifier):
                         mismatch_dict[header_modifier] += 1
                    else:
                         mismatch_dict[header_modifier] = 1


                #reset lists to collect next records
                current_record_r1 = []
                current_record_r2 = []
                current_record_r3 = []
                current_record_r4 = []

           #update file line counter
            i+=1


#CLOSE FILES
for x in output_file_handles:
     output_file_handles[x].close()



#write summary report file
total_record_count = matched_index_counter + index_hopper_counter + unknown_index_counter
with open('summary.txt', 'w') as summary_file:
     print(f'Total Records: {total_record_count}', file=summary_file)
     print(f'Number of Matched Indexes: {matched_index_counter} ({round(matched_index_counter/total_record_count*100, 2)}%)', file=summary_file)
     print(f'Number of Index Hopped Pairs: {index_hopper_counter} ({round(index_hopper_counter/total_record_count*100, 2)}%)', file=summary_file)
     print(f'Number of Unkown Index Pairs: {unknown_index_counter} ({round(unknown_index_counter/total_record_count*100, 2)}%)', file=summary_file)
     print(f'Index Quality Score Threshold: {quality_threshold}', file=summary_file)

     print('---------------------------------------------------', file=summary_file)
     print('Counts of Matched Indexes:', file=summary_file)
     
     sorted_matches = dict(sorted(match_dict.items(), key=lambda k: k[1], reverse=True))
     for x in sorted_matches:
          print(f'{x}: {sorted_matches[x]} ({round(sorted_matches[x]/total_record_count*100, 2)}%) ', file=summary_file)
     
     print('---------------------------------------------------', file=summary_file)
     print('Counts of Index Hopped Pairs:', file=summary_file)
     
     sorted_mismatches = dict(sorted(mismatch_dict.items(), key=lambda k: k[1], reverse=True))
     for x in sorted_mismatches:
          print(f'{x}: {sorted_mismatches[x]} ({round(sorted_mismatches[x]/total_record_count*100, 2)}%) ', file=summary_file)
     