#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:45:00 2018

@author: anna

This script reads trees in newick format, identifies Ophiostoma lineages, 
and identifies O. novo-ulmi strains which cluster together with O. ulmi,
to be later filtered out from fasta alignments

USAGE:

python filter_ULMLike.py
requirements:       Species_distinction_PCA.txt
                    RAxML_bipartitions_50coll.Gene_<gene_name>
output:             plotTree_good_novoulmi.csv
                    plotTree_bad_novoulmi.csv
                    plotTree_pop_counts.txt
"""


#from ete3 import Tree, TreeStyle, TextFace, NodeStyle, AttrFace, faces
from ete3 import Tree
import pandas as pd
import glob


### Getting list of all novo-ulmi, and ulmi

samples_fn = "Species_distinction_PCA.txt"
df = pd.read_csv(samples_fn,sep='\t',header=0)
ame1_df = df.loc[df["GroupSamtools"]=="A1"]
ame2_df = df.loc[df["GroupSamtools"]=="A2"]
nov_df = df.loc[df["GroupSamtools"]=="N"]
ulm_df = df.loc[df["GroupSamtools"]=="U"]
him_df = df.loc[df["GroupSamtools"]=="H"]
novoulmi = pd.concat([ame1_df, ame2_df, nov_df])
novoulmi_list = list(novoulmi['SamName'])
nov_list = list(nov_df['SamName'])
ame1_list = list(ame1_df['SamName'])
ame2_list = list(ame2_df['SamName'])
u_list = list(ulm_df['SamName'])
h_list = list(him_df['SamName'])
ulm_list = u_list + h_list
all_samples = novoulmi_list + ulm_list



### Getting closest relatives of one leaf (It means getting all descendants of a parental node of a leaf)

def getRelatedLeafs(node):
    nodeup = node.up
    n = nodeup.get_descendants()
    N = []
    for nod in n:
        if nod.is_leaf():
            N.append(nod.name)
    return N, nodeup


### Getting trees

trees = glob.glob("RAxML_bipartitions_50coll.Gene_*")

### Loop for trees and find ULM-like novo-ulmi strains 

T = {}
Good_samples = []
Pop_counts = []
for tree2 in trees:
    
    win = tree2
    print(win)
       
    t = Tree(tree2,format=2)
    t.get_leaf_names()
    
    ### Loop for each novo-ulmi
        
    ulm_like_strains = []
    for nu in novoulmi_list:
        #print(nu)   
        # eaxample for one sample
        #nu = "H332_N"
        
        ### in a loop collect strain names as long as they are >=8
        node = t&nu
        relatives, uppernode = getRelatedLeafs(node)
        leaf_nodes = relatives
        
        while len(leaf_nodes) <8:
            #print(len(leaf_nodes))
            #print(uppernode.name)
            relatives, uppernode = getRelatedLeafs(uppernode)
            leaf_nodes = relatives
        
        ### count how many novo-ulmi and ulmi there are
        if nu in leaf_nodes:
            leaf_nodes.remove(nu)
        
        ulmi = set(leaf_nodes) & set(ulm_list)
        nulm = set(leaf_nodes) & set(novoulmi_list)
        frac_ulm = float(len(ulmi))/(len(leaf_nodes))
        
        ### condition if there is more than 50% of ulmi among relatives - mark this strain as ulmi-like
        if frac_ulm > 0.5:
            ulm_like_strains.append(nu)
    
        T[win] = ulm_like_strains
    
    ### collect only good samples
    good_samples = list(set(all_samples) - set(ulm_like_strains))
    Good_samples.append([win] + good_samples)
    
    ### collect info about lineage counts after filtering [gene Name1 Name2 Nnov Nulm Nhim]
    lin_ame1 = list(set(ame1_list) - set(ulm_like_strains))
    lin_ame2 = list(set(ame2_list) - set(ulm_like_strains))
    lin_nov = list(set(nov_list) - set(ulm_like_strains))

    pop_counts = [str(ele) for ele in [win, len(lin_ame1), len(lin_ame2), len(lin_nov), len(u_list), len(h_list)]]
    Pop_counts.append(pop_counts)


### writing bad strains - these novo-ulmi which are ulm-like

wh = open("plotTree_bad_novoulmi.csv","w")
for ele in T.keys():
    if T[ele]:
        f = ','.join([str(i) for i in T[ele]])
        wh.write(ele+','+f+'\n')
wh.flush()
wh.close()

### writing good strains - these novo-ulmi which are not ulm-like + the rest

wh2 = open("plotTree_good_novoulmi.csv","w")
for ele in Good_samples:
    f = ','.join(ele)
    wh2.write(f+'\n')
wh2.flush()
wh2.close()

### writing pop count info

wh3 = open("plotTree_pop_counts.txt","w")
wh3.write("Tree\tAME1\tAME2\tNOV\tULM\tHIM\n")
for ele in Pop_counts:
    f = '\t'.join(ele)
    wh3.write(f+'\n')
wh3.flush()
wh3.close()




