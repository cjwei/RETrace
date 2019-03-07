#!/usr/bin/env python3
import argparse
from tqdm import tqdm
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
import pandas as pd

'''
This script is made specifically for our known ex vivo HCT116 tree.  It will be used to determine the accuracy of any phylogenetic tree we calculate.  To do this, we need to input the following files:
    1) sample_info = tab-delimited file containing file location of sample bam, sample name [same as leaves in tree file], and sample type (i.e. clone in ex vivo tree [2-1-G10_3-1-A2, 2-1-H7_3-6-C6, 2-1-G10_3-1-B1, 2-2-B1_3-2-A6])
    2) exVivo_dist = csv file containing matrix of distances of each pair of clones to the root (in the case of our HCT116 ex vivo tree, root is 0-1-E6)
    2) newick_tree = file containing Newick tree outputted from runHipSTR.py
    3) prefix = output prefix for error calculation statistics comparing calculated Newick tree to given ex vivo tree
We will then use the above input to calculate an errorDict which contains the following structure:
    errorDict
        cell_div = reference cell division difference between nodes (ex: [2-1-G10_3-1-A2, 2-1-G10_3-1-B1, 2-2-B1_3-2-A6] = abs(max()))
            "Correct" = number of triplets correct
            "Total" = total number of triplets analyzed

'''

def calcTriplets(NJTree, exvivo_pd, sampleDict):
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
                    if tree_dist.index(max(tree_dist)) == ref_dist.index(max(ref_dist)):
                        errorDict[max(ref_dist)]["Correct"] += 1
                    errorDict[max(ref_dist)]["Total"] += 1
        for dist in sorted(errorDict.keys()):
            print(str(dist) + "\t" + str(errorDict[dist]["Correct"]) + "\t" + str(errorDict[dist]["Total"]) + "\n")
    return errorDict

def main():
    parser = argparse.ArgumentParser(description="Evaluate phylogenetic tree versus known ex vivo tree by comparing all triplets")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delmited file containing sample information (save as used for runHipSTR.py)")
    parser.add_argument('--exVivo', action="store", dest="exVivo_dist", help="csv file containing pairwise cell division distances across ex vivo tree")
    parser.add_argument('--tree', action="store", dest="tree_file", help="Newick tree file")
    args = parser.parse_args()

    #We want to first parse the sample_info file
    sampleDict = {}
    with open(args.sample_info) as f:
        for line in f:
            if len(line.split()) == 3: #If clone is specified
                (bam, sample, clone) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
                sampleDict[sample]["clone"] = clone
            elif len(line.split()) == 2: #If clone is not specified
                (bam, sampleName) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
            else:
                print("Incorrect formmating for line (bam, sampleName, [optional] clone):\n" + "\t".join(line.split))
                return

    #Import exVivo tree distances into pandas dataframe
    exVivo_pd = pd.read_csv(args.exVivo_dist, delimiter=',', index_col=0)

    #Import Newick tree and calculate all triplets of leaves in tree [sample] (with at least two clones per triple).  Determine whether distances between leaves is correct based on exVivo_pd
    with open(args.tree_file, 'r') as f:
        newick_tree = f.read().replace("\n",'')
    NJTree = Tree(newick_tree)
    errorDict = calcTriplets(NJTree, exVivo_pd, sampleDict)

if __name__ == "__main__":
    main()
