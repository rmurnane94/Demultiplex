This algorithm now has to go through and pull out the index pairs that are matched up, the index pairs that are not matched up, and the index pairs where one or both are unknown/low quality (aren't in known indexes, or are below our quality score cutoff). The information for the reads is split between the 4 files. Once the algorithm figures out the situation with the index pairs, the algorithm must add records appropriately to separate output files. There also should be a file with some output stats such as counts of the different types of index pairs etc. 

This means also have to determine a quality score cutoff. TBD

Some barcodes have 'N' in their sequences which are undetermined bases with low quality scores. These barcodes will be categorized as  unknown


might want to make a list/dictionary of all files to read so they are easy to reference: file_list = (R1, R2, R3, R4 ) at the vary least assign variables to each.

Make a dictionary of all indexes: index_dict with 'B1':	'GTAGCGTA' etc. 

Make a function for finding reverse compliments of DNA sequences. The indexes in the R2 file will be forward and match the known indexes. The indexes in R3 will be reverse compliments.
ex: for a match of 'B1' the index  'GTAGCGTA'  will be found in R2, but its reverse compliment 'TACGCTAC' will be found in R3. Having a function will make this a lot easier.
```
def reverse_compliment(sequence: str) -> str:
	'''takes the DNA sequence and makes reverse compliment'''
	return reverse compliment
Example input: AGCT
Example output: AGCT
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
