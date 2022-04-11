'''
Input: filtered reads from mirna_counts_mmu_variants_prefilter.py with first 18 bases match
    something in the variants reference dataset

Output: counts file looking at only variants and parent miRNAs

after each round, remove reads from list

find parent matches
find perfect matches
strip to 2 less than the total length of the reference, find matches
    any ambiguities will be assigned to the parent if is one of the options, or split over all possibilities equally
'''


import sys

try:
	reads = open(sys.argv[1],'rU') # file with trimmed,processed reads
except IndexError:
	print("No reads file provided")
	sys.exit()

infilename = reads.name
outfilename = infilename[:-13] + 'counts_var.txt'
outfile = open(outfilename, 'w')

mirna_sequences = open("/lab/solexa_bartel/mailejim/miRNA_mmu_variants_reference.txt",'rU') # file with mirna names / sequences of interest

LastLine = False
sequence_d = dict() #declare dictionary matching miRNA name (key) to sequence (value)
counts_d = dict() #declare dictionary matching miRNA name (key) to count number (value)
reads_list = [] # working list of reads, removed from in each step

# create dictionary of sequences and dictionary of counts, initialize counts at 0
while LastLine == False:
    line = mirna_sequences.readline()
    if line == '':
        LastLine = True
    else:
        line = line.rstrip('\n')
        line_split = line.split('\t')
        mirna_name = line_split[0]
        sequence = line_split[1]
        sequence_d[mirna_name] = sequence
        counts_d[mirna_name] = float(0)
LastLine = False

total = 0
# find perfect matches to parent and variant  miRNAs
while LastLine == False:
    line = reads.readline()
    if line == '':
        LastLine = True
    else:
        
        line = line.rstrip('\n')
        if line in sequence_d.values():
            i = 0
            # count number of times sequence appears
            for n in sequence_d.values():
                if line == n:
                    i += 1
            # add fractionally
            for n in sequence_d:
                if sequence_d[n] == line:
                    counts_d[n] = counts_d[n] + float(1/i)
        else:
            total += 1
            reads_list.append(line)
LastLine = False

#total = 0 # unmatched reads
'''
# count imperfect matches 
for line in reads_list:
    parent = False
    i = 0
    # count number of times sequence appears, and if one of the options is a parent miRNA
    for n in sequence_d:
        sequence = sequence_d[n]
        if sequence[:-2] == line[:len(sequence[:-2])+1]:
            i += 1
            if n[-1] not in 'isd':
                parent = n
    if parent != False:
        counts_d[parent] = counts_d[parent] + 1
    else:
        # add fractionally
        for n in sequence_d:
            if sequence_d[n][:-2] == line[:len(sequence_d[n][:-2])+1]:
                counts_d[n] = counts_d[n] + float(1/i)
    if i == 0:
        total += 1
'''        
for mirna in counts_d:
    outfile.write(mirna + "\t" + str(counts_d[mirna]) + '\n') 
outfile.write(str(total))
outfile.close()







 
