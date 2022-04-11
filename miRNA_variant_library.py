#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 14:30:28 2021

@author: mailejim
"""


from Bio.Seq import Seq
from Bio.Seq import MutableSeq


def get_opposite(base):
    base_seq = Seq(base)
    return str(base_seq.reverse_complement())


miRNA_table = open("miRNA_table.txt", "r").readlines()
miRNA_dict = {}
for i in range(len(miRNA_table)):
    miRNA_table[i] = miRNA_table[i].split("\t")
    miRNA_table[i][1] = Seq(miRNA_table[i][1].strip("\n"))
    miRNA_dict[miRNA_table[i][0]] = miRNA_table[i][1]


constant1 = Seq("CTGCATGTATCAATTCTTGTGTATCTCGCCG")
constant2 = Seq("GTGTTATTCTTG")
constant3 = Seq("TGGCGCGGTACACTTATTCGTCCAATTTTTC")

miRNA_list = ["miR-7a-5p", "miR-154-3p", "miR-543-3p", "miR-495-3p"]


miRNA_guide_list  = []  
sequences = []
for ii in miRNA_list:
    miRNA_guide = miRNA_dict[ii]
    print(miRNA_guide)
    miRNA_guide_list.append(miRNA_guide)
    
    # 1 nt substitions
    for i in range(1,len(miRNA_guide)):
        for n in ["A", "T", "G", "C"]:    
            miRNA_guide_list.append(miRNA_guide[0:i] + Seq(get_opposite(n)) + miRNA_guide[i+1:])

           
    # 2 nt substitutions
    for i in range(len(miRNA_guide)-1): 
        for n in ["AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC", "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC"]:
            miRNA_guide_list.append(miRNA_guide[0:i] + Seq(get_opposite(n)) + miRNA_guide[i+2:])
    
    if len(miRNA_guide) == 24:
        # deletions
        for i in range(len(miRNA_guide)):
            if i != len(miRNA_guide)-1:
                miRNA_guide_list.append(miRNA_guide[0:i] + miRNA_guide[i+1:])
            else:
                miRNA_guide_list.append(miRNA_guide[0:i])
        
        for i in range(len(miRNA_guide)-1):
            if i != len(miRNA_guide)-2:
                miRNA_guide_list.append(miRNA_guide[0:i] + miRNA_guide[i+2:])
            else:
                miRNA_guide_list.append(miRNA_guide[0:i])
    
    if len(miRNA_guide) == 22:           
        # insertions
        for i in range(2,len(miRNA_guide)+1):
            for base in ["A", "T", "C", "G"]:
                miRNA_guide_list.append(miRNA_guide[0:i] + Seq(base) + miRNA_guide[i:])
        
        for i in range(2,len(miRNA_guide)-1):
            for base in ["AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC", "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC"]:
                miRNA_guide_list.append(miRNA_guide[0:i] + Seq(base) + miRNA_guide[i:])
 
print(len(miRNA_guide_list))

def bulge(sequence, position):
    #print(sequence, position)
    if str(Seq(sequence[position]).reverse_complement()) != "G":
        return sequence[0:position] + sequence[position+1:]
    else:

        if str(Seq(sequence[position-1]).reverse_complement()) != "G":
            
            return sequence[0:position-1] + sequence[position:]
        
        elif str(Seq(sequence[position+1]).reverse_complement()) != "G":
            return sequence[0:position+1] + sequence[position+2:]
        
        elif str(Seq(sequence[position-2]).reverse_complement()) != "G":
            return sequence[0:position-2] + sequence[position-1:]
    return 'ERROR'

for miRNA_guide in miRNA_guide_list:
    miRNA_star = MutableSeq(str(miRNA_guide.reverse_complement()))
    
    if len(miRNA_star) == 22:
        miRNA_star  = miRNA_star[0:6] + MutableSeq(str(Seq(miRNA_star[6]).reverse_complement())) + miRNA_star[7:] # bulge the 6th position
        miRNA_star  = miRNA_star[0:12] + MutableSeq(str(Seq(miRNA_star[12]).reverse_complement())) + miRNA_star[13:] # bulge the 12th position
        miRNA_star  = miRNA_star[0:17] + MutableSeq(str(Seq(miRNA_star[17]).reverse_complement())) + miRNA_star[18:] # bulge the 17th position

    if len(miRNA_star) == 23:
        miRNA_star  = bulge(miRNA_star, 6) #miRNA_star[0:6] + miRNA_star[7:] # delete the 6th position
        miRNA_star  = miRNA_star[0:13] + MutableSeq(str(Seq(miRNA_star[13]).reverse_complement())) + miRNA_star[14:] # bulge the 13th position
        miRNA_star  = miRNA_star[0:18] + MutableSeq(str(Seq(miRNA_star[18]).reverse_complement())) + miRNA_star[19:] # bulge the 18th position

    if len(miRNA_star) == 24:
        miRNA_star  = bulge(miRNA_star, 6) #miRNA_star[0:6] + miRNA_star[7:] # delete the 6th position
        miRNA_star  = bulge(miRNA_star, 14) # miRNA_star[0:14] + miRNA_star[15:] # delete the 14th position
        miRNA_star  = miRNA_star[0:19] + MutableSeq(str(Seq(miRNA_star[19]).reverse_complement())) + miRNA_star[20:] # bulge the 19th position

    sequences.append(constant1 + miRNA_guide + constant2 + miRNA_star + constant3)
 

n = 0
lines_seen = set()
outfile = open("miRNA_seq_variants_no_duplicates.txt", "w")
for line in sequences:
    if line not in lines_seen: # not a duplicate
        outfile.write(str(line)+"\n")
        lines_seen.add(line)
        n +=1
outfile.close()
print(n)


