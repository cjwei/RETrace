#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict
from tqdm import tqdm
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import multiprocessing
import random
manager = multiprocessing.Manager()

'''
This script is similar to the evalPhylo.py script in the main RETrace program.  However, the goal of this script is to compare number of targets per sample versus sample accuracy (as calculated by triplet accuracy).
Ideally, what we should see is there should be a positive correlation between number of targets captured per sample versus sample accuracy.
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

def multi_calcTriplets(triplet_group, NJTree, exvivo_pd, sampleDict, tripletDict):
    for triplet in tqdm(triplet_group):
        tree_dist = []
        ref_dist = []
        for sample_pair in ([0,1], [0,2], [1,2]):
            MRCA = NJTree.get_common_ancestor([triplet[sample_pair[0]], triplet[sample_pair[1]]])
            MRCA_dist = NJTree.get_distance(MRCA)
            tree_dist.append(MRCA_dist)
            ref_dist.append(exvivo_pd.loc[sampleDict[triplet[sample_pair[0]]]["clone"]][sampleDict[triplet[sample_pair[1]]]["clone"]])
        #We want to save the count information in encoded string treeDict[triplet]: ref_dist,[1 (correct) or 0 (incorrect)]
        if tree_dist.index(max(tree_dist)) == ref_dist.index(max(ref_dist)): #If using MRCA depth as comparison metric
            tripletDict[triplet] = str(max(ref_dist)) + ",1"
        else:
            tripletDict[triplet] = str(max(ref_dist)) + ",0"
    return

def plotTriplets(tripletDict, prefix):
    errorDict = {}
    for triplet in tqdm(sorted(tripletDict.keys())):
        (ref_dist, correct_bool) = tripletDict[triplet].split(',')
        if ref_dist not in errorDict.keys():
            errorDict[ref_dist] = {}
            errorDict[ref_dist]["Correct"] = 0
            errorDict[ref_dist]["Total"] = 0
        errorDict[ref_dist]["Correct"] += int(correct_bool)
        errorDict[ref_dist]["Total"] += 1
    #Plot propotion of correct triplets per distance metric
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

    #Write correct triplet rate into text file
    f_output = open(prefix + ".evalPhylo.txt", 'w')
    f_output.write("\n".join([str(dist_list[i]) + "\t" + str(corr_list[i]) for i in range(len(dist_list))]) + "\n")
    f_output.close()

    return errorDict

def evalPhylo(sample_info, prefix, exVivo_dist, tree_file, nproc):
    sampleDict = import_sampleDict(sample_info)

    #Import exVivo tree distances into pandas dataframe
    exVivo_pd = pd.read_csv(exVivo_dist, delimiter=',', index_col=0)

    #Import Newick tree and calculate all triplets of leaves in tree [sample] (with at least two clones per triple).  Determine whether distances between leaves is correct based on exVivo_pd
    with open(tree_file, 'r') as f:
        newick_tree = f.read().replace("\n",'')
    NJTree = Tree(newick_tree)

    #Create list of all triplets and run through to determine whether each triplet is correct
    tripletDict = manager.dict() #This allows for parallel processing of triplet errors for multiple triplets at once
    jobs = []
    triplet_set = set()
    print("Naming all triplets in tree")
    for sample1 in tqdm(sorted(sampleDict.keys())):
        clone1 = sampleDict[sample1]["clone"]
        for sample2 in sorted(sampleDict.keys()):
            clone2 = sampleDict[sample2]["clone"]
            for sample3 in sorted(sampleDict.keys()):
                clone3 = sampleDict[sample3]["clone"]
                triplet = tuple(sorted([sample1, sample2, sample3]))
                if len(set([clone1, clone2, clone3])) >= 2 and len(set(triplet)) == 3:
                # if len(set([clone1, clone2, clone3])) == 3 and triplet not in triplet_list: #We want only triplets where each of the three leaves stem from different clones
                    triplet_set.add(triplet)
    triplet_list = list(triplet_set)

    print("Calculating correct triplet rate")
    random.shuffle(triplet_list) #We want to randomize in order to even processing time
    for triplet_group in [triplet_list[i::nproc] for i in range(nproc)]:
        p = multiprocessing.Process(target = multi_calcTriplets, args = (triplet_group, NJTree, exVivo_pd, sampleDict, tripletDict))
        jobs.append(p)
        p.start()
    #Join tripletDict
    for p in jobs:
        p.join()

    #Plot number of correct triplets
    print("Plotting/printing correct triplet rate")
    plotTriplets(tripletDict, prefix)

if __name__ == "__main__":
    main()
