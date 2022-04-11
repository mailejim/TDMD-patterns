#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 10:51:59 2021

@author: mailejim
"""
 
import sys

try:
        mirna_sequences = open(sys.argv[1]) # file with miRNA counts assigned to variants
except IndexError:
        print("No counts file provided")
        sys.exit()


infilename = mirna_sequences.name
outfilename = infilename[:-4] + '_pooled.txt'
outfile = open(outfilename, 'w')
    
# iterate through input file to find the variants
LastLine = False
output = []
variants = {}
while LastLine == False:
    line = mirna_sequences.readline()
    if line == '':
        LastLine = True
    else:
        line = line.rstrip('\n')
        line_split = line.split('\t')
        mirna_name, count = line_split[0], line_split[1]
        
        # add non-variants to the output file
        if mirna_name[-1] != "i" and mirna_name[-1] != "s" and mirna_name[-1] != "d": 
            outfile.write(line + '\n')
        
        # add variants to variants dictionary
        else:
            variants[mirna_name] = count
            

pooled_dict = {} # mutation type as keys, values as number of counts
mutations = []
parent_mirnas = ["miR-7a-5p", "miR-154-3p", "miR-543-3p", "miR-495-3p"]


# find all possible mutations
for n in range(1,25): #each position along the miRNA
    for a in 'ATGC':
        for b in 'ATCG':
            mutations.append(a + str(n) + b + 's')
        mutations.append(a + str(n) + 'i')
        mutations.append(a + str(n) + 'd')

for var in variants: # go through each variant
    for parent in parent_mirnas: # go through each parent mirna
        for mut in mutations: # go through each possible mutation the variant could have
            if parent in var and mut in var:
                
                if mut[-1] == 's':
                    if str(parent) + '-' + mut[1:-2] + mut[-1] in pooled_dict:
                        pooled_dict[str(parent) + '-' + mut[1:-2]+ mut[-1]] = float(pooled_dict[str(parent) + '-' + mut[1:-2]+ mut[-1]]) + float(variants[var])
                    else:
                        pooled_dict[str(parent) + '-' + mut[1:-2]+ mut[-1]] = variants[var]
                else: #insertion or deletion 
                    if str(parent) + '-' + mut[1:] in pooled_dict:
                        pooled_dict[str(parent) + '-' + mut[1:]] = float(pooled_dict[str(parent) + '-' + mut[1:]]) + float(variants[var])
                    else:
                        pooled_dict[str(parent) + '-' + mut[1:]] = variants[var]
             
                
                
for var in pooled_dict:

       outfile.write(str(var) + "\t" + str(pooled_dict[var]) + '\n')

        
outfile.close()
 

