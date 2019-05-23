#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 15:56:16 2018

@author: anna

This script reads a table with per site estimates of Pi, performs reshuffling 
of all positions in the genome and recalulates Pi for each window;
after x times, 95th and 99th percentile of bootstrapped Pi distribution is obtained


USAGE:
    
python bootstrapPi.py

requirements:   sample_table_tP_persite.txt
output:         bootstrapPi.out

"""


import pandas as pd
import numpy as np


#chroms = {'OphioH327chr_1':6937932,'OphioH327chr_2':6817711,'OphioH327chr_3':3669772,'OphioH327chr_4':3419703,\
#          'OphioH327chr_5':2848703,'OphioH327chr_6':2801594,'OphioH327chr_7':2758224,'OphioH327chr_8':2531247}

chroms = {'OphioH327chr_1':400000}

winSize=5000
stepSize=5000

replicates = 100


#dd = pd.read_parquet('../tables_angsd/ame1_table_tP_persite.txt',engine='fastparquet')

print("Getting number of snps in each window\n")

tab = pd.read_csv("sample_table_tP_persite.txt",sep="\t",header=None,names=["Chr","WinCenter","Site","Pi"])


Nsnp = []
i=0
for chr in sorted(chroms.keys()):
    print(chr)
    x=1
    y=winSize
    step = stepSize
    CH = tab[(tab['Chr'] == chr)]
    while y < chroms[chr]:
        #print(x)
        nf = CH[(CH['WinCenter'] >= x) & (CH['WinCenter'] <= y)]
        f = chr, x, y
        fbed = chr, x-1, y, i,len(nf)
        Nsnp.append(fbed)
        x+=step
        y+=step
        i+=1
 

OUT = []
for repl in list(range(replicates)):
    print("\nReplicate "+str(repl)+"\n")
    #ame1 = pd.read_table("../tables_angsd/ame1_table_tP_persite.txt",sep="\t",header=None,names=["Chr","WinCenter","Site","Pi"])
    print("Sampling")
    df = tab.sample(frac=1).reset_index(drop=True)
    #D = tab.merge(df, how='outer', left_index=True, right_index=True)
    print("Merging")
    tab["Pi_rand"] = pd.Series(df["Pi"],index=tab.index)
    
    R = []
    P = []

    i=0
    for chr in sorted(chroms.keys()):
        print(chr)
        x=1
        y=winSize
        step = stepSize
        CH = tab[(tab['Chr'] == chr)]
        while y < chroms[chr]:
            #print(x)
            nf = CH[(CH['WinCenter'] >= x) & (CH['WinCenter'] <= y)]
            stat = nf["Pi_rand"].mean()
            f = chr, x, y
            fbed = chr, x-1, y, stat
            R.append(fbed)
            P.append(stat)
            x+=step
            y+=step
            i+=1
            
    OUT.append(P)
    
    
newO = zip(*OUT)

M = []
for ele in list(newO):
    perc = np.percentile(ele,[95,99])
    M.append(perc)
    
wh = open("bootstrapPi.out","w")
for ele in list(range(len(Nsnp))):
    wh.write(Nsnp[ele][0]+"\t"+str(Nsnp[ele][1])+"\t"+str(Nsnp[ele][2])+"\t"+str(M[ele][0])+"\t"+str(M[ele][1])+"\n")
    
wh.flush()
wh.close()
    
    
    
