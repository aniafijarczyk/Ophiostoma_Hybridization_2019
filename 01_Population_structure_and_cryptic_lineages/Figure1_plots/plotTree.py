#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:45:00 2018

@author: anna
"""


from ete3 import Tree, TreeStyle, NodeStyle, faces
import pandas as pd
import glob


trees = glob.glob("Newick_Export_ordered_root_75coll.nwk")


for tree2 in trees:

    win = tree2.split("_")[1]
    print(win)
    
    t = Tree(tree2,format=2)
    t.get_leaf_names()
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)
    
    samples_fn = "Species_distinction_PCA.txt"
    df = pd.read_csv(samples_fn, sep = '\t', dtype= {'NameBam':'str','Region':'str'},keep_default_na=False)
    ame1_df = df.loc[df["GroupSamtools"]=="A1"]
    ame2_df = df.loc[df["GroupSamtools"]=="A2"]
    nov_df = df.loc[df["GroupSamtools"]=="N"]
    ulm_df = df.loc[df["GroupSamtools"]=="U"]
    him_df = df.loc[df["GroupSamtools"]=="H"]

    regs = list(df.Region)
    regcols = []
    for ele in regs:
        if ele == "EA":
            regcols.append("black")
        elif ele == "NA":
            regcols.append("grey")
        else:
            regcols.append("lightgrey")
    df['RegCols'] = regcols

    # General tree outline    
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.tree_width = 7500
    ts.draw_guiding_lines = True
    

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
    A1 = [t.get_leaves_by_name(ele)[0] for ele in ame1_df.NameBam]
    for n in A1:
        n.set_style(A1style)
    
    A2style = NodeStyle()
    A2style["shape"] = "circle"
    A2style["fgcolor"] = "#EEB422"
    A2style["size"] = 14
    A2 = [t.get_leaves_by_name(ele)[0] for ele in ame2_df.NameBam]
    for n in A2:
        n.set_style(A2style)
    
    Nstyle = NodeStyle()
    Nstyle["shape"] = "circle"
    Nstyle["fgcolor"] = "#698b22"
    Nstyle["size"] = 14
    N = [t.get_leaves_by_name(ele)[0] for ele in nov_df.NameBam]
    for n in N:
        n.set_style(Nstyle)
    
    Ustyle = NodeStyle()
    Ustyle["shape"] = "circle"
    Ustyle["fgcolor"] = "#1874cd"
    Ustyle["size"] = 14
    U = [t.get_leaves_by_name(ele)[0] for ele in ulm_df.NameBam]
    for n in U:
        n.set_style(Ustyle)


    def mylayout(node):
        if node.is_leaf():
            row = df.loc[df["NameBam"]==node.name]
            ind = row.index[0]
            regioncol = row.at[ind,'RegCols']
            shapeFace = faces.RectFace(width=12, height=12, fgcolor='white', bgcolor=regioncol)
            faces.add_face_to_node(shapeFace, node, column=0, position="aligned")
                
    ts.draw_guiding_lines=True
    ts.layout_fn = mylayout
    t.render("Tree_"+win+".pdf",tree_style=ts)
    
    leaf_names = t.get_leaf_names()
    wh = open('plotTree_strain_order.txt','w')
    f = '\n'.join(leaf_names)
    wh.write(f+'\n')
    wh.flush()
    wh.close()