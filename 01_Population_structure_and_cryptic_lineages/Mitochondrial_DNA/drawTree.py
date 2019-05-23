from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import pandas as pd


# Reading population info
samples_fn = "Species_distinction_PCA.txt"
df = pd.read_csv(samples_fn,sep='\t')
ame1_df = df.loc[df["GroupSamtools"]=="A1"]
ame2_df = df.loc[df["GroupSamtools"]=="A2"]
nov_df = df.loc[df["GroupSamtools"]=="N"]
ulm_df = df.loc[df["GroupSamtools"]=="U"]
him_df = df.loc[df["GroupSamtools"]=="H"]

# Reading tree
tree_collapsed_50 = "RAxML_bipartitions_50coll.mtDNA"
t = Tree(tree_collapsed_50,format=2)

# General tree outline
#anc = t.get_common_ancestor("HP30_H","DDS82_A2")
ts = TreeStyle()
ts.show_leaf_name = True
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



#t.show(tree_style=ts)

# Saving tree
#t.render("mtDNATree.png",tree_style=ts,units="mm",w=300)
t.render("mtDNATree.pdf",tree_style=ts)
