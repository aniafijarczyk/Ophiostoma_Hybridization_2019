#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:45:00 2018

@author: anna
"""


from ete3 import Tree, TreeStyle, NodeStyle, faces, PhyloNode
import pandas as pd
import glob


#trees = glob.glob("Newick_Export_ordered_root_75coll.nwk")
#trees = glob.glob("concatenateGenes_tree_doublepres")
trees = glob.glob("NewickExport_fasttree.nwk")

for tree2 in trees:

    win = tree2.split("_")[1]
    print(win)
    
   
    t = Tree(tree2,format=2)
    t.get_leaf_names()
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    #R = t.get_common_ancestor("HP30_H","HP32_H")
    t.set_outgroup(R)
    
    # changing bootstraps & clustering
    for i in t.traverse():
        #supp = i.support
        i.support = i.support*100
        if i.support < 75.0:
            i.dist = 0
            i.support = 0
    #t.write(format=2,outfile="test_tree")
    

    samples_fn = "Species_distinction_PCA.txt"
    df = pd.read_csv(samples_fn, sep = '\t', dtype= {'SamName':'str','Region':'str'},keep_default_na=False)
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
        elif ele == "NZ":
            regcols.append("lightgrey")
        elif ele == "unk":
            regcols.append("white")
    df['RegCols'] = regcols

    # General tree outline    
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.tree_width = 300
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
    A1style["fgcolor"] = "#d55e00"
    A1style["size"] = 14
    A1 = [t.get_leaves_by_name(ele)[0] for ele in ame1_df.SamName]
    for n in A1:
        n.set_style(A1style)
    
    A2style = NodeStyle()
    A2style["shape"] = "circle"
    A2style["fgcolor"] = "#f0e442"
    A2style["size"] = 14
    A2 = [t.get_leaves_by_name(ele)[0] for ele in ame2_df.SamName]
    for n in A2:
        n.set_style(A2style)
    
    Nstyle = NodeStyle()
    Nstyle["shape"] = "circle"
    Nstyle["fgcolor"] = "#009e73"
    Nstyle["size"] = 14
    N = [t.get_leaves_by_name(ele)[0] for ele in nov_df.SamName]
    for n in N:
        n.set_style(Nstyle)
    
    Ustyle = NodeStyle()
    Ustyle["shape"] = "circle"
    Ustyle["fgcolor"] = "#0072b2"
    Ustyle["size"] = 14
    U = [t.get_leaves_by_name(ele)[0] for ele in ulm_df.SamName]
    for n in U:
        n.set_style(Ustyle)

    Hstyle = NodeStyle()
    Hstyle["shape"] = "circle"
    Hstyle["fgcolor"] = "#cc79a7"
    Hstyle["size"] = 14
    H = [t.get_leaves_by_name(ele)[0] for ele in him_df.SamName]
    for n in H:
        n.set_style(Hstyle)
        
        
    def mylayout(node):
        if node.is_leaf():
            row = df.loc[df["SamName"]==node.name]
            ind = row.index[0]
            regioncol = row.at[ind,'RegCols']
            if regioncol == 'black':
                shapeFace = faces.RectFace(width=24, height=12, fgcolor='white', bgcolor=regioncol)
                faces.add_face_to_node(shapeFace, node, column=0, position="aligned")
            elif regioncol == 'grey':
                #shapeFace = faces.CircleFace(radius=6,color='lightgrey')
                #faces.add_face_to_node(shapeFace, node, column=0, position="aligned")
                shapeFace = faces.RectFace(width=24, height=12, fgcolor='white', bgcolor="lightgrey")
                faces.add_face_to_node(shapeFace, node, column=0, position="aligned")                
            elif regioncol == 'lightgrey':
                shapeFace = faces.TextFace(text = "NZ",fsize=10, ftype = "Arial")
                faces.add_face_to_node(shapeFace, node, column=0, position="aligned")
                #shapeFace = faces.RectFace(width=24, height=12, fgcolor='black', bgcolor="white")
                #faces.add_face_to_node(shapeFace, node, column=0, position="aligned")
            elif regioncol == 'white':
                shapeFace = faces.RectFace(width=24, height=12, fgcolor='white', bgcolor="white")
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
