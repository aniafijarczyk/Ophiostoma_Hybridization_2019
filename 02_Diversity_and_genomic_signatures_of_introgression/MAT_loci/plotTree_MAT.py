#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:45:00 2018

@author: anna
"""

from ete3 import Tree, TreeStyle, NodeStyle
import pandas as pd

tree = "RAxML_bipartitions_50coll.Gene_1981842_OphioH327chr_2"
#tree = "RAxML_bipartitions_50coll.Gene_1981846_OphioH327chr_2"

t = Tree(tree,format=2)
t.get_leaf_names()


samples_fn = "Species_distinction_PCA.txt"
df = pd.read_csv(samples_fn,sep='\t',header=0)
ame1_df = df.loc[df["GroupSamtools"]=="A1"]
ame2_df = df.loc[df["GroupSamtools"]=="A2"]
nov_df = df.loc[df["GroupSamtools"]=="N"]
ulm_df = df.loc[df["GroupSamtools"]=="U"]
him_df = df.loc[df["GroupSamtools"]=="H"]


# General tree outline
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_support = True
ts.tree_width = 550


# General style for the nodes
for n in t.traverse():
   nstyle = NodeStyle()
   nstyle["shape"] = "circle"
   nstyle["fgcolor"] = "black"
   nstyle["size"] = 0
   n.set_style(nstyle)

# Node style for each lineage
A1style = NodeStyle()
A1style["shape"] = "circle"
A1style["fgcolor"] = "#CD2626"
A1style["size"] = 14
A1 = [t.get_leaves_by_name(ele)[0] for ele in ame1_df.SamName]
for n in A1:
    n.set_style(A1style)

A2style = NodeStyle()
A2style["shape"] = "circle"
A2style["fgcolor"] = "#EEB422"
A2style["size"] = 14
A2 = [t.get_leaves_by_name(ele)[0] for ele in ame2_df.SamName]
for n in A2:
    n.set_style(A2style)

Nstyle = NodeStyle()
Nstyle["shape"] = "circle"
Nstyle["fgcolor"] = "#698b22"
Nstyle["size"] = 14
N = [t.get_leaves_by_name(ele)[0] for ele in nov_df.SamName]
for n in N:
    n.set_style(Nstyle)

Ustyle = NodeStyle()
Ustyle["shape"] = "circle"
Ustyle["fgcolor"] = "#1874cd"
Ustyle["size"] = 14
U = [t.get_leaves_by_name(ele)[0] for ele in ulm_df.SamName]
for n in U:
    n.set_style(Ustyle)

Hstyle = NodeStyle()
Hstyle["shape"] = "circle"
Hstyle["fgcolor"] = "#5d478b"
Hstyle["size"] = 14
H = [t.get_leaves_by_name(ele)[0] for ele in him_df.SamName]
for n in H:
    n.set_style(Hstyle)


treeID = tree.split('_')[3]
t.render("Tree_"+treeID+".pdf",tree_style=ts)
