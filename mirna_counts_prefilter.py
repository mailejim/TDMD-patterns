'''
Input: reads file

Output: filtered reads file with only reads that have the first 18 bases match
    with something in the variants reference dataset
'''

import sys

try:
	reads = open(sys.argv[1],'rU') # file with trimmed,processed reads
except IndexError:
	print("No reads file provided")
	sys.exit()

infilename = reads.name
outfilename = infilename[:-4] + '_prefilter.txt'
outfile = open(outfilename, 'w')

LastLine = False
sequence_list = [] #list of sequences that appear in the matching dictionary

mirna_sequences = open("/lab/solexa_bartel/mailejim/miRNA_mmu_variants_reference.txt",'rU') # file with mirna names / sequences of interest

# references
while LastLine == False:
    line = mirna_sequences.readline()
    if line == '':
        LastLine = True
    else:
        line = line.rstrip('\n')
        line_split = line.split('\t')
        mirna_name = line_split[0]
        sequence = line_split[1][0:18]
        sequence_list.append(sequence) 
LastLine = False

print(sequence_list)

# input file is reads
while LastLine == False:
    line = reads.readline()
    if line == '':
        LastLine = True
    else:
        line = line.rstrip('\n')
        input_sequence = str(line[0:18].rstrip()) # get first 18 characters
        if input_sequence in sequence_list:
            outfile.write(line + '\n')        

outfile.close()
 
