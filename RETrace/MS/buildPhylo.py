#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict, import_targetDict
import pickle
import more_itertools
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
from tqdm import tqdm
import numpy as np

def calcBulk(alleleDict, targetDict):
    '''
    Cluster all alleles belonging to given target across all sample.  Use these as possible "allele_groups"
    '''
    for target_id in sorted(alleleDict.keys()):
        if target_id in targetDict.keys(): #We only want to analyze target_id specified in targetDict
            sub_len = len(targetDict[target_id]["sub_seq"])
            all_alleles = set()
            for sample in sorted(alleleDict[target_id]["sample"].keys()):
                try:
                    all_alleles.update([int(allele/sub_len) for allele in alleleDict[target_id]["sample"][sample]["allelotype"]]) #Convert alleles to number of subunits (because we want to group all consecutive integer number subunit alleles)
                except:
                    continue
            alleleDict[target_id]["allele_groups"] = {}
            for indx, group in enumerate(more_itertools.consecutive_groups(sorted(set(all_alleles)))): #Need to use "set(filtered_alleles)" in order to prevent repeated int in filtered_alleles
                alleleDict[target_id]["allele_groups"][indx] = [allele * sub_len for allele in list(group)] #Save allele groups in terms of raw number of bases difference from ref
    return alleleDict

def calcDist(alleleDict, distDict, sample_pair, sample1, sample2, shared_targets, dist_metric):
    '''
    Calculate distance between two samples given their shared_targets and dist_metric
    '''
    total_dist = 0
    num_alleles = 0
    for target_id in shared_targets:
        target_dist = 0
        for allele_indx, allele_group in alleleDict[target_id]["allele_groups"].items():
            allelotype1 = list(set(alleleDict[target_id]["sample"][sample1]["allelotype"]).intersection(set(allele_group)))
            allelotype2 = list(set(alleleDict[target_id]["sample"][sample2]["allelotype"]).intersection(set(allele_group)))
            if len(allelotype1) == 2 and len(allelotype2) == 2:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]) + abs(allelotype1[1] - allelotype2[1]), abs(allelotype1[1] - allelotype2[0]) + abs(allelotype1[0]- allelotype2[1]))
                elif dist_metric == "EqorNot":
                    target_dist += int(len(set(allelotype1).symmetric_difference(set(allelotype2)))/2)
                num_alleles += 2
            elif len(allelotype1) == 1 and len(allelotype2) == 2:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]), abs(allelotype1[0] - allelotype2[1]))
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0] and allelotype1[0] != allelotype2[1]:
                        target_dist += 1
                num_alleles += 1
            elif len(allelotype1) == 2 and len(allelotype2) == 1:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]), abs(allelotype1[1] - allelotype2[0]))
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0] and allelotype1[1] != allelotype2[0]:
                        target_dist += 1
                num_alleles += 1
            elif len(allelotype1) == 1 and len(allelotype2) == 1:
                if dist_metric == "Abs":
                    target_dist += abs(allelotype1[0] - allelotype2[0])
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0]:
                        target_dist += 1
                num_alleles += 1
        total_dist += target_dist
    distDict["sampleComp"][sample_pair] = {}
    distDict["sampleComp"][sample_pair]["dist"] = float(total_dist/num_alleles)
    distDict["sampleComp"][sample_pair]["num_targets"] = len(shared_targets)

    return distDict

def makeDistMatrix(sample_list, sharedDict, alleleDict, dist_metric):
    '''
    Wrapper for calculating distance matrix (saved within distDict), which has the following structure:
        distDict
            "sampleComp"
                sample1
                    sample2
                        "dist"
                        "num_targets"
    '''
    distDict = {}
    distDict["sampleComp"] = {}
    #Calculate pairwise distance for all samples depending on shared targets (follow order specified in target_list [esp for bootstrapping, which may have duplicates due to sampling w/ replacement])
    for sample1 in sorted(sample_list):
        for sample2 in sorted(sample_list):
            sample_pair = tuple(sorted([sample1,sample2]))
            if sample_pair not in distDict["sampleComp"].keys():
                distDict = calcDist(alleleDict, distDict, sample_pair, sample1, sample2, sharedDict[sample_pair], dist_metric)
    return distDict

def drawTree(distDict, sample_list, outgroup, prefix, bootstrap):
    '''
    Run neighbor-joining phylogenetic tree building algorithm on pairwise cell distance (saved in distDict)
    '''
    distMatrix = []
    targetMatrix = []
    for sample1 in sorted(sample_list):
        sample1_dist = []
        sample1_targets = []
        for sample2 in sorted(sample_list):
            sample_pair = tuple(sorted([sample1, sample2]))
            sample1_dist.append(distDict["sampleComp"][sample_pair]["dist"])
            sample1_targets.append(distDict["sampleComp"][sample_pair]["num_targets"])
        distMatrix.append(sample1_dist)
        targetMatrix.append(sample1_targets)
    if bootstrap is False: #Only output statistics for distance and number targets shared if for original tree (don't output for bootstrap resampling)
        statsOutput = open(prefix + ".buildPhylo.stats.txt", 'w')
        for dist_indx,dist_list in enumerate(distMatrix): #Print matrix containing distances
            statsOutput.write(sorted(sample_list)[dist_indx] + "," + ",".join(str(round(i,3)) for i in dist_list) + "\n")
        for target_indx,target_list in enumerate(targetMatrix): #Print matrix containing number targets shared between each pair
            statsOutput.write(sorted(sample_list)[target_indx] + "," + ",".join(str(j) for j in target_list) + "\n")
        statsOutput.close()
    distObj = DistanceMatrix(distMatrix,sorted(sample_list))
    skbio_tree = nj(distObj, result_constructor=str)
    ete_tree = Tree(skbio_tree) #We use skbio to first make a tree from distance matrix then convert to ete tree
    if outgroup is "NA":
        return ete_tree
    else:
        if outgroup == "Midpoint":
            tree_midpoint = ete_tree.get_midpoint_outgroup()
            ete_tree.set_outgroup(tree_midpoint)
        else:
            ete_tree.set_outgroup(outgroup)
    return ete_tree

def bootstrapTree(nodeDict, treeTemp, sample_list):
    '''
    Evaluate whether original node (keys of nodeDict) are found within treeTemp (tree from bootstrapped samples).
    Be sure to only consider sample_name within sample_list, which may be different than original because of bootstrap resampling
    '''
    #We want to create a list containing all nodes found in treeTemp
    tempNodes = set()
    for node in treeTemp.search_nodes():
        leaf_list = []
        for leaf in node:
            leaf_list.append(leaf.name)
        tempNodes.add(tuple(sorted(set(leaf_list))))
    #We next want to loop through all nodes in the original tree to determine if it is found within bootstrapped tree (disregard sample_names thrown out)
    for node in nodeDict.keys():
        node_intersect = tuple(sorted(set(node).intersection(sample_list)))
        if node_intersect == node:
            nodeDict[node]["Num_sampled"] += 1
            if node_intersect in tempNodes:
                nodeDict[node]["Num_verified"] += 1
    return nodeDict

def buildPhylo(sample_info, prefix, target_info, alleleDict_file, dist_metric, outgroup, bootstrap):
    '''
    Draw phylogenetic tree based on calculated allelotype of single cells.  This is done by:
        1) Filtering out only likely alleles by creating a "pseudo"-bulk in which we cluster all single cell alleles together (calcBulk)
        2) Make distance matrix usin given dist_metric (Abs or EqorNot)
        3) Draw tree using neighbor-joining method
        4) Repeating 2-3 but bootstrapping the samples in order to determine the support of each node of the original tree
    '''

    sampleDict = import_sampleDict(sample_info)
    sample_list = sorted(sampleDict.keys()) #We want to create a list containing all sample_name (for original tree prior to bootstrapping)

    targetDict = import_targetDict(target_info)

    #Calculate "bulk" allelotype counts by merging all single cell calls
    alleleDict = pickle.load(open(alleleDict_file, 'rb'))
    alleleDict = calcBulk(alleleDict, targetDict)

    #Pre-calculate shared targets between each sample
    sharedDict = {} #Contains all shared targets between each pairwise sample
    for sample1 in sorted(sample_list):
        for sample2 in sorted(sample_list):
            sample_pair = tuple(sorted([sample1,sample2]))
            if sample_pair not in sharedDict.keys():
                shared_targets = []
                for target_id in sorted(alleleDict.keys()):
                    if sample1 in alleleDict[target_id]["sample"].keys() and sample2 in alleleDict[target_id]["sample"].keys() and target_id in targetDict.keys(): #Only want to analyze target_id in targetDict
                        if "allelotype" in alleleDict[target_id]["sample"][sample1].keys() and "allelotype" in alleleDict[target_id]["sample"][sample2].keys():
                            shared_targets.append(target_id)
                sharedDict[sample_pair] = shared_targets

    #Calculate original tree using all samples found within sampleDict
    distDict_original = makeDistMatrix(sample_list, sharedDict, alleleDict, dist_metric) #Calculate pairwise distance between each sample
    tree_original = drawTree(distDict_original, sample_list, outgroup, prefix, False) #Draw neighbor joining tree.  We want to declare Fase for bootstrap because we want to output stats file for original tree
    f_tree_original = open(prefix + '.buildPhylo.newick-original.txt', 'w')
    f_tree_original.write(tree_original.write(format = 0))
    f_tree_original.write("\n")
    f_tree_original.close()
    #Bootstrap resample to create new distance matrices/trees and add support values to internal nodes of original tree
    if bootstrap is True:
        #Determine dictionary of tree nodes from tree_original
        nodeDict = {}
        for node in tree_original.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name) #We need to compare leaves in a node cluster without regard to order
            nodeDict[tuple(sorted(leaf_list))] = {}
            nodeDict[tuple(sorted(leaf_list))]["Num_verified"] = 0 #Contains values for number of times random bootstrap tree (tree_temp) contains given node
            nodeDict[tuple(sorted(leaf_list))]["Num_sampled"] = 0 #Contains number of times the node occured during bootstrap resampling
            # nodeDict[tuple(sorted(leaf_list))]["NodeID"] = node.write(format = 9)
        for i in tqdm(range(10000)): #Bootstrap resample 10,000 times
            #Random downsample from pool of available targets (target_list) to use for distance calculation
            sample_list = set(np.random.choice(list(sampleDict.keys()), len(sampleDict.keys()), replace=True))
            if outgroup not in ["Midpoint", "NA"]: #We want to make sure our outgroup remains in the tree even during bootstraping
                sample_list.add(outgroup)
            distDict_temp = makeDistMatrix(sample_list, sharedDict, alleleDict, dist_metric)
            tree_temp = drawTree(distDict_temp, sample_list, outgroup, prefix, bootstrap)
            nodeDict = bootstrapTree(nodeDict, tree_temp, sample_list) #Determine whether each node in original tree is found in tree_temp
        #Add support information to original tree
        for node in tree_original.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name)
            if nodeDict[tuple(sorted(leaf_list))]["Num_sampled"] > 0:
                node_support = float(nodeDict[tuple(sorted(leaf_list))]["Num_verified"] / nodeDict[tuple(sorted(leaf_list))]["Num_sampled"])
            else:
                node_support = 2.0 #Assign nodes that were not present in any bootstrap simulation a value of 2.0
            node.add_features(support = node_support)
        #Output tree with optional support values
        f_tree_bootstrap = open(prefix + '.buildPhylo.newick-bootstrap.txt', 'w')
        f_tree_bootstrap.write(tree_original.write(format = 0))
        f_tree_bootstrap.write("\n")
        f_tree_bootstrap.close()
