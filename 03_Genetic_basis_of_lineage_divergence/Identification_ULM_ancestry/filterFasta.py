#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 17:02:37 2019

@author: anna


This script reads fasta alignments and a table with list of strains to select 
from each alignment, and filters out only strains which are in the list.

USAGE:

python filterFasta.py
requirements:   plotTree_good_novoulmi.csv
                jgi.p|Ophnu1|<gene_name>_OphioH327chr_1_flt.fasta
output:         Gene_<gene_name>.fasta

"""

import glob
from Bio import SeqIO


fasta = glob.glob("*.fasta")

def readMe(fname):
    fh = open(fname, 'r')
    linie = fh.readlines()
    k = [ele.split()[0] for ele in linie]
    d = [ele.split(',') for ele in k]
    return d

if __name__ == '__main__':
    
    r = readMe("plotTree_good_novoulmi.csv")
    D =  {a[0].split("_")[3]:a[1:] for a in r}
    
    for plik in fasta:
        geneID = plik.split('|')[-1].split("_")[0]
        print(geneID)
        if geneID in D.keys():
            gList = D[geneID]
            R = []
            for record in SeqIO.parse(plik, "fasta"):
                name = record.id
                if name in gList:
                    R.append(record)
    
            SeqIO.write(R,"Gene_"+geneID+".fasta","fasta")