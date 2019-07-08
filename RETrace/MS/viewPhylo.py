#!/usr/bin/env python3
import os
os.environ['QT_QPA_PLATFORM']='offscreen' #Fix for remote ssh tree render <https://github.com/etetoolkit/ete/issues/387>
from RETrace.MS.utilities import import_sampleDict
from ete3 import Tree, TreeStyle, TextFace, NodeStyle #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>

def viewPhylo(sample_info, tree_file, prefix, bootstrap):
    '''
    This script will utilize ete3 to draw the phylogenetic tree calculated from buildPhylo
    '''

    #Import the samples that we used for buildPhylo
    sampleDict = import_sampleDict(sample_info)
    MStree = Tree(tree_file)
    #We want to assign colors to the leaves of the tree (this is based on <https://stackoverflow.com/questions/39380907/how-to-color-leaves-on-ete3-tree-python-3>)
    for node in MStree.traverse():
        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 2
        nstyle["hz_line_width"] = 2
        # node.img_style['size'] = 0 #Hide circles at nodes
        if node.is_leaf():
            nstyle = NodeStyle()
            color = sampleDict[node.name]["clone_color"]
            nstyle["fgcolor"] = color
            nstyle["size"] = 15
            node.set_style(nstyle)
        else:
            nstyle["size"] = 0
        node.img_style = nstyle
            # node.img_style["bgcolor"] = color
            # name_face = TextFace(node.name, fgcolor=color, fsize=10)
            # node.add_face(name_face, column=0, position='branch-right')

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_leaf_name = True
    # ts.mode = "c"
    # ts.arc_start = -180 # 0 degrees = 3 o'clock
    # ts.arc_span = 180
    if bootstrap is True:
        ts.show_branch_support = True
    MStree.render(prefix + ".viewPhylo.pdf", tree_style=ts)
