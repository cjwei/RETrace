#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict
from ete3 import Tree, TreeStyle #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>

def viewPhylo.py(sample_info, tree_file, prefix, bootstrap):
    '''
    This script will utilize ete3 to draw the phylogenetic tree calculated from buildPhylo
    '''

    #Import the samples that we used for buildPhylo
    sampleDict = import_sampleDict(sample_info)
    MStree = Tree(tree_file)
    ts = TreeStyle()
    ts.show_leaf_name = True
    if bootstrap is True:
        ts.show_branch_support = True
    MStree.render(prefix + ".buildPhylo.svg", tree_style=ts)
