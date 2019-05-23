#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:16:03 2019

@author: anna

This script reads a table with SNPs with chromosome, position, allele frequency (AF) 
and genotype coded as 0, 1 or ./. for missing data. 
According to AF it converts the genotype to F (frequent - major allele), R (rare - minor allele).


USAGE:

python convertMajor2Minor.py
requires:   genotypes1.tab
output:     convert2MajorMinor.tab
"""


#filename = "sample_genotypes.tab"
filename = "genotypes1.tab"

def recode(biglist,switch):
    newlist = []
    if switch < 0.5:
        for ele in biglist:
            if ele == './.':
                #p = 'Missing'
                p = 'M'
            elif ele == "0":
                #p = "Major"
                p = 'F'
            elif ele == "1":
                #p = "Minor"
                p = 'R'
            newlist.append(p)
    elif switch >= 0.5:
        for ele in biglist:
            if ele == './.':
                #p = 'Missing'
                p = 'M'
            elif ele == '1':
                #p = 'Major'
                p = 'F'
            elif ele == '0':
                #p = 'Minor'
                p = 'R'
            newlist.append(p)        
    return newlist


def readMe(fname):
    fh = open(fname,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    n = [(i[0],i[1],i[2],i[3:]) for i in k]
    return n


if __name__ == '__main__':
    r = readMe(filename)
    #print(r)
    N = []
    eh = open('convert2MajorMinor.tab','w')
    i = 1
    for ele in r:
        af = float(ele[2])
        recoded = recode(ele[3],af)
        recoded_2w = '\t'.join(recoded)
        #N.append((ele[0],ele[1],recoded))
        eh.write(ele[0]+'\t'+ele[1]+'\t'+str(i)+'\t'+recoded_2w+'\n')
        i+=1
    eh.flush()
    eh.close()
    #print(N)
