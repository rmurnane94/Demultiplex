
**Goals**: Our goal is to look through a lane of sequencing generated from the 2017 BGMP cohort’s library preps and determine the level of index swapping and undetermined index-pairs, before and after quality filtering of index reads. In order to do this, we must first demultiplex the data. In Assignment the First, we will develop a strategy to de-multiplex samples to create 48 FASTQ files that contain acceptable index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped reads-pairs, and two FASTQ files undetermined (non-matching or low quality) index-pairs.

indexes:
```
B1	GTAGCGTA    A5	CGATCGAT    C1	GATCAAGG
B9	AACAGCGA    C9	TAGCCATG    C3	CGGTAATC
B3	CTCTGGAT    C4	TACCGGAT    A11	CTAGCTCA
C7	CACTTCAC    B2	GCTACTCT    A1	ACGATCAG
B7	TATGGCAC    A3	TGTTCCGT    B4	GTCCTAAG
A12	TCGACAAG    C10	TCTTCGAC    A2	ATCATGCG
C2	ATCGTGGT    A10	TCGAGAGT    B8	TCGGATTC
A7	GATCTTGC    B10	AGAGTCCA    A8	AGGATAGC
```

4 fastq files:
```
/projects/bgmp/shared/2017_sequencing/

1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
```

1. Perform some initial data exploration! Record any bash commands you used inside a lab notebook (submit to this repo!).
    1. Determine which files contain the indexes, and which contain the paired end reads containing the biological data of interest. Create a table and label each file with either read1, read2, index1, or index2.

2. Determine the length of the reads in each file.
    1. Determine the phred encoding for these data.
    2. 

running this command for each file to see what the records look like
```zcat 1294_S1_L008_R4_001.fastq.gz | head -4```
 
R1 and R4 files have long sequences, while R2 and R3 have short sequences. R1 and R4 are going to the the read files, while R2 and R3 are the index files. 

running this for each file to determine the sequence lengths
``` zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc -c ```

R1 and R4 have 101 (102 characters -1 for the \n)
R2 and R3 have 8 (9-1)

from the same command as above heading each file, there are a lot of J quality scores for each file. This means the encoding is probably Phred 33 since that would mean the scores were all 10 if it was Phred 64 encoding and that is super low. There are also no lower case letters present which would be the case if it was Phred 64 encoding. 

| File name                    | label  | Read length | Phred encoding |
| ---------------------------- | ------ | ----------- | -------------- |
| 1294_S1_L008_R1_001.fastq.gz | read1  | 101         | 33             |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8           | 33             |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8           | 33             |
| 1294_S1_L008_R4_001.fastq.gz | read2  | 101         | 33             |

Finishing Assignment the first Part 1:

Creating plots for the distributions of the mean quality score per base for each of R1, R2, R3, R4

ran an sbatch script score_distribution_sbatch.sh to run score_distribution.py.

worked through errors, using gzip to read zipped files
runtime results:
```
	User time (seconds): 9737.68
	System time (seconds): 4.35
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:42:45		Exit status: 0

```


Part 2 stuff

Write up a strategy (**NOT A SCRIPT**) for writing an algorithm to de-multiplex files and reporting index-hopping. That is, given four input FASTQ files (2 with biological reads, 2 with index reads) and the 24 known indexes above, demultiplex reads by index-pair, outputting:

- one R1 FASTQ file and one R2 FASTQ file **per** matching index-pair,
- another two FASTQ files for non-matching index-pairs (index-hopping), and
- two additional FASTQ files when one or both index reads are unknown or low quality (do not match the 24 known indexes [this includes indexes with 'N's in them] or do not meet a quality score cutoff)


more initial exploration:
all files have the same number of lines
identifiers from the R1 file match R2 file, match R3 file etc.
all records in each file are in the same order. ie the first record in R1 corresponds to first record in R2,R3,R4



START OF SCRIPT
This algorithm now has to go through and pull out the index pairs that are matched up, the index pairs that are not matched up, and the index pairs where one or both are unknown/low quality (aren't in known indexes, or are below our quality score cutoff)

This means also have to determine a quality score cutoff. TBD

Some barcodes have 'N' in their sequences which are undetermined bases with low quality scores. These barcodes will be categorized as  unknown


might want to make a list/dictionary of all files to read so they are easy to reference: file_list = (R1, R2, R3, R4 ) at the vary least assign variables to each.

Make a dictionary of all indexes: index_dict with 'B1':	'GTAGCGTA' etc. 

Make a function for finding reverse compliments of DNA sequences. The indexes in the R2 file will be forward and match the known indexes. The indexes in R3 will be reverse compliments.
ex: for a match of 'B1' the index  'GTAGCGTA'  will be found in R2, but its reverse compliment 'TACGCTAC' will be found in R3. Having a function will make this a lot easier.
```
def reverse_compliment(sequence: str)
	'takes the DNA sequence and makes reverse compliment'
	return reverse compliment
```

initialize counters to count the number of read-pairs with matched indexes (one for each index pair), the number of read pairs with index-hopping, and the number of read pairs with unknown indexes. 

will want to loop through all output files and empty them in case there is leftovers from last time

With preliminary steps out of the way, could go through the files looking for each index and trying to find matches, but that would necessitate going through the files multiple times. probably best to just go through the files ONCE and then assign records it goes to their new files.

need to open all the files at the same time though
with open(r1, 'r') as file1... do this for all files

since the identifiers between the files match up, can go through and capture each 4 line record from each file at a time

Save each 4 lines of these records in lists r1_record_list, r2_record_list etc.

then go through a bunch of if statements to help categorize which file the records will be output in
if the  index in r2 record list or in r3 record list is not a known index in the dictionary-> unknown
if either index is below the quality threshold -> unknown
if the r2 index doesn't match r3 index (reverse compliment): -> nonmatching
if r2 index matches r3 index (reverse compliment): -> matching

depending on what happens in each case, the next steps will be a little different

in the case they are matching, will need to figure out what the actual index is that is matching. save that as a variable for the current matching index

what output file names are going to look like:
B1_R1 etc for all indexes 
B1_R2 etc for all indexes
nonmatch R1
nonmatch R2
unknown R1
unknown R2

for writing the outputs can use 'a' to append to files without having to rewrite every time
ex if it is unknown: with open (unknown R1, 'a')
then write the record from the saved r1 record list
make sure to add the 'index from r2-index from r3' to the header of the read
repeat for r2


for the index matching files, just use f strings that add in the current matching index into file name ( ex: f'{}_ R1') and write the record r1 and r2 like that in the same manner as above here

also update counters for the properly matched read pairs, index-hopping read pairs, and number of pairs with unknown indexes depending on what the situation is


-----------------------------
submitted the above pseudocode for feedback

"I would recommend creating a dictionary that can track all of the matching, swapping, and unknown indexes. Also, you will need a function that can find the average quality score of the sequence so you can check if it is below the threshold. You did not include an example input and output for the reverse complement function. The logic behind the body of the function is sound, and will properly put the records into their corresponding output files. The only issue is that you should make sure to close the files you are writing to, so you don’t have a bunch of files open at the same time."


"- You don’t need to name the files with the name of the index. Leslie said that we should name the files with the index sequences themselves (i.e. `GTAGCGTA_R1.fq` and `GTAGCGTA_R2.fq`). So it’s not necessary to create the dictionary of indexes.
- In addition to counting the total number of index-hopped reads, you will need to report the counts for each possible permutation of indexes. There are a total of 24*24 possible permutations (24 of them are matched indexes), which is a lot, so it might be useful to have a dictionary to store the counts.
- If you are going to have a quality score cutoff, it would be useful to have a function that will decode the quality scores of the indexes and determine if the cutoff is met or not."


"1. Consider creating a dedicated function to convert Phred quality scores and apply a threshold check. This will make your code cleaner and easier to maintain when evaluating index quality.
2. Be sure to modify the read headers to include barcode information at the time of extraction, regardless of whether the index pair is valid, unknown, or hopped. This ensures all records are appropriately annotated before classification.
3. Don’t forget to close all file handles after writing—especially important when working with multiple output files to prevent data loss or corruption.
4. The reverse complement function should include example inputs and outputs, as explicitly requested in the assignment. This helps with clarity and debugging.
5. Overall, your logic is solid and well-structured. Just tighten up a few implementation details for robustness and completeness!"




Takeaways:

create dictionary and track the matching/swapping/unknown indexes x2

function that can find the average quality score of index sequences to check if they are below threshold x3

for my functions I outline, include input and output examples x2

make sure to CLOSE files that are being written x2

can name files with indexes themselves, do not need to find the name of the index

I mentioned header modifying but might not be explicit enough..do not forget it



Assignment the Third:
Using all my ideas and incorporating feedback for finishing demultiplexing

Per Leslie it turns out it is a better idea to open all files all together once to write to them rather than opening each file independently each time. 

Decided to use a dictionary to hold all the open files so they can easily be accessed. 

Got test version to seemingly produce right outputs for test files

Additional feedback from Leslie: pull indexes from file, not hardcoded

Updating to include dictionary to hold all index pairs counted

Updating code to pull indexes from index file so they don't have to be hardcoded.

Ran on full dataset:
```
	Command being timed: "./demultiplex.py"
	User time (seconds): 4143.48
	System time (seconds): 101.52
	Percent of CPU this job got: 34%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:23:51
	Exit status: 0
```
fixed bug that included output files 'index_R1.fq' and 'index_R2.fq'. this was from pulling the header row from the index file. got rid of it

Talked to Leslie about progress. things to add/fix:
percentages
add the matched pairs to top
use a set, not a list for indexes
use functions for repetitive parts of code

Using a SET, not a list to hold them since it is hashable and easier to search. do not need dictionary since we do not have key-value pairs. updating code to change list to a set.

writing function to write to files to make the code easier to read
the writing part of the code is very repetitive. updated code with the write_output function.

Adding two dictionaries instead of one to keep track of the index pairs. adding one for matches and one for mismatches. this will split my current dictionary of all matches and make it easier to print 

added to summary.txt file to add all summary statistics

ran sbatch script
```
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

```

usr/bin/time -v results
```
User time (seconds): 4385.45
System time (seconds): 54.28
Percent of CPU this job got: 35%
Elapsed (wall clock) time (h:mm:ss or m:ss): 3:29:34

Exit status: 0
```
