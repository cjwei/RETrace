#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict
from tqdm import tqdm
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

def calcTriplets(NJTree, exvivo_pd, sampleDict, prefix):
    analyzed_triplets = []
    errorDict = {}
    #Iterate through all nodes in tree and calculate all possible triplets (where each triple comes from at least two different clones)
    for sample1 in tqdm(sorted(sampleDict.keys())):
        tripletDict = {}
        clone1 = sampleDict[sample1]["clone"]
        for sample2 in sorted(sampleDict.keys()):
            clone2 = sampleDict[sample2]["clone"]
            for sample3 in sorted(sampleDict.keys()):
                clone3 = sampleDict[sample3]["clone"]
                triplet = tuple(sorted([sample1, sample2, sample3]))
                if len(set([clone1, clone2, clone3])) >= 2 and triplet not in analyzed_triplets and len(set(triplet)) == 3:
                # if len(set([clone1, clone2, clone3])) == 3 and triplet not in analyzed_triplets: #We want only triplets where each of the three leaves stem from different clones
                    analyzed_triplets.append(triplet)
                    if triplet not in tripletDict.keys():
                        tripletDict[triplet] = []
                    tree_dist = []
                    ref_dist = []
                    for sample_pair in ([0,1], [0,2], [1,2]):
                        MRCA = NJTree.get_common_ancestor([triplet[sample_pair[0]], triplet[sample_pair[1]]])
                        MRCA_dist = NJTree.get_distance(MRCA)
                        tree_dist.append(MRCA_dist)
                        ref_dist.append(exvivo_pd.loc[sampleDict[triplet[sample_pair[0]]]["clone"]][sampleDict[triplet[sample_pair[1]]]["clone"]])
                    if max(ref_dist) not in errorDict.keys():
                        errorDict[max(ref_dist)] = {}
                        errorDict[max(ref_dist)]["Correct"] = 0
                        errorDict[max(ref_dist)]["Total"] = 0
                    if tree_dist.index(max(tree_dist)) == ref_dist.index(max(ref_dist)): #If using MRCA depth as comparison metric
                        errorDict[max(ref_dist)]["Correct"] += 1 #If using MRCA depth as comparison metric
                    errorDict[max(ref_dist)]["Total"] += 1 #If using MRCA depth as comparison metric
    #Plot proportion of correct triplets per distance metric
    dist_list = sorted(errorDict.keys())
    corr_list = []
    [total_corr, total_triplets] = [0,0]
    for dist in dist_list:
        corr_list.append(errorDict[dist]["Correct"]/errorDict[dist]["Total"])
        total_corr += errorDict[dist]["Correct"]
        total_triplets += errorDict[dist]["Total"]
    dist_list.append("All Triplets")
    corr_list.append(total_corr/total_triplets)
    corr_df = pd.DataFrame({"Distance": dist_list, "Correct Triplets": corr_list})
    sns_barplot = sns.barplot(x="Distance", y="Correct Triplets", data=corr_df)
    fig = sns_barplot.get_figure()
    fig.savefig(prefix + ".evalPhylo.eps")
    return errorDict

def evalPhylo(sample_info, prefix, exVivo_dist, tree_file):
    '''
    This script is made specifically for our known ex vivo HCT116 tree.  It will be used to determine the accuracy of any phylogenetic tree we calculate.  To do this, we need to input the following files:
        1) sample_info = tab-delimited file containing file location of sample bam, sample name [same as leaves in tree file], sex, and sample type (i.e. clone in ex vivo tree [2-1-G10_3-1-A2, 2-1-H7_3-6-C6, 2-1-G10_3-1-B1, 2-2-B1_3-2-A6])
        2) exVivo_dist = csv file containing MRCA distance from root, as approximated in units of cell divisions
        2) newick_tree = file containing Newick tree output
        3) prefix = output prefix for error calculation statistics comparing calculated Newick tree to given ex vivo tree
    We will then use the above input to calculate an errorDict which contains the following structure:
        errorDict
            cell_div = reference cell division difference between nodes (ex: [2-1-G10_3-1-A2, 2-1-G10_3-1-B1, 2-2-B1_3-2-A6] = abs(max()))
                "Correct" = number of triplets correct
                "Total" = total number of triplets analyzed
    '''
    sampleDict = import_sampleDict(sample_info)

    #Import exVivo tree distances into pandas dataframe
    exVivo_pd = pd.read_csv(exVivo_dist, delimiter=',', index_col=0)

    #Import Newick tree and calculate all triplets of leaves in tree [sample] (with at least two clones per triple).  Determine whether distances between leaves is correct based on exVivo_pd
    with open(tree_file, 'r') as f:
        newick_tree = f.read().replace("\n",'')
    NJTree = Tree(newick_tree)
    errorDict = calcTriplets(NJTree, exVivo_pd, sampleDict, prefix)

if __name__ == "__main__":
    main()
