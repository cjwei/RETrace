#!/usr/bin/env python3
import argparse
from importer import import_sampleDict
import pickle
import more_itertools
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
from tqdm import tqdm
import numpy as np

def calcBulk(alleleDict, targetDict):
    for target_id in sorted(alleleDict.keys()):
        sub_len = len(targetDict[target_id]["sub_seq"])
        all_alleles = set()
        for sample in sorted(alleleDict[target_id]["sample"].keys()):
            all_alleles.update([int(allele/sub_len) for allele in alleleDict[target_id]["sample"][sample]]) #Convert alleles to number of subunits (because we want to group all consecutive integer number subunit alleles)
        alleleDict[target_id]["allele_groups"] = {}
        for indx, group in enumerate(more_itertools.consecutive_groups(sorted(set(all_alleles)))): #Need to use "set(filtered_alleles)" in order to prevent repeated int in filtered_alleles
            alleleDict[target_id]["allele_groups"][indx] = [allele * sub_len for allele in list(group)] #Save allele groups in terms of raw number of bases difference from ref
    return alleleDict

def calcDist(alleleDict, distDict, sample1, sample2, shared_targets, dist_metric):
    total_dist = 0
    num_alleles = 0
    for target_id in shared_targets:
        target_dist = 0
        for allele_indx, allele_group in alleleDict[target_id]["allele_groups"].items():
            allelotype1 = list(set(alleleDict[target_id]["sample"][sample1]).intersection(set(allele_group)))
            allelotype2 = list(set(alleleDict[target_id]["sample"][sample2]).intersection(set(allele_group)))
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
        if target_id not in distDict["targetComp"].keys():
            distDict["targetComp"][target_id] = {}
        distDict["targetComp"][target_id][tuple(sorted([sample1,sample2]))] = target_dist #Save contribution to distance from each target_id
        total_dist += target_dist
    distDict["sampleComp"][sample1][sample2]["dist"] = float(total_dist/num_alleles)
    distDict["sampleComp"][sample1][sample2]["num_targets"] = len(shared_targets)
    return distDict

def makeDistMatrix(sampleDict, alleleDict, target_list, dist_metric, bootstrap):
    distDict = {}
    distDict["sampleComp"] = {}
    distDict["targetComp"] = {}
    #Calculate pairwise distance for all samples depending on shared targets (follow order specified in target_list [esp for bootstrapping, which may have duplicates due to sampling w/ replacement])
    for sample1 in tqdm(sorted(sampleDict.keys())):
        distDict["sampleComp"][sample1] = {}
        for sample2 in sorted(sampleDict.keys()):
            distDict["sampleComp"][sample1][sample2] = {}
            shared_targets = set()
            for target_id in sorted(alleleDict.keys()):
                if sample1 in alleleDict[target_id]["sample"].keys() and sample2 in alleleDict[target_id]["sample"].keys() and target_id in target_list:
                    shared_targets.add(target_id)
            #Use calcDist function to calculate genotype_dist between samples and contribution of each target to distance
            distDict = calcDist(alleleDict, distDict, sample1, sample2, shared_targets, dist_metric)
    return distDict

def drawTree(distDict, sampleDict, outgroup, prefix, bootstrap):
    distMatrix = []
    targetMatrix = []
    for sample1 in sorted(sampleDict.keys()):
        sample1_dist = []
        sample1_targets = []
        for sample2 in sorted(sampleDict.keys()):
            sample1_dist.append(distDict["sampleComp"][sample1][sample2]["dist"])
            sample1_targets.append(distDict["sampleComp"][sample1][sample2]["num_targets"])
        distMatrix.append(sample1_dist)
        targetMatrix.append(sample1_targets)
    if bootstrap is False: #Only output statistics for distance and number targets shared if for original tree (don't output for bootstrap resampling)
        statsOutput = open(prefix + ".stats.out", 'w')
        for dist_indx,dist_list in enumerate(distMatrix): #Print matrix containing distances
            statsOutput.write(sorted(sampleDict.keys())[dist_indx] + "," + ",".join(str(round(i,3)) for i in dist_list) + "\n")
        for target_indx,target_list in enumerate(targetMatrix): #Print matrix containing number targets shared between each pair
            statsOutput.write(sorted(sampleDict.keys())[target_indx] + "," + ",".join(str(j) for j in target_list) + "\n")
        statsOutput.close()
    distObj = DistanceMatrix(distMatrix,sorted(sampleDict.keys()))
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

def bootstrapTree(nodeDict, treeTemp):
    for temp_node in treeTemp.search_nodes():
        leaf_list = []
        for leaf in temp_node:
            leaf_list.append(leaf.name)
        if tuple(sorted(leaf_list)) in nodeDict.keys():
            nodeDict[tuple(sorted(leaf_list))]["Bootstrap"] += 1
    return nodeDict

def buildPhylo():
    parser = argparse.ArgumentParser(description="Draw phylogenetic tree based on calculated allelotype of single cells")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specify output prefix (for targetDict and alleleDict, along with any stats or plot files)")
    parser.add_argument('--targets', action="store", dest="target_file", default="All", help="[Optional] Specify list of targets used for distance calculation")
    parser.add_argument('--dist', action="store", dest="dist_metric", help="Specify distance metric for pairwise comparisons [Abs, EqorNot]", default="EqorNot")
    parser.add_argument('--outgroup', action="store", dest="outgroup", default="NA", help="[Optional] Specify outgroup for rooted NJ tree (if use midpoint, specify 'Midpoint')")
    parser.add_argument('-bootstrap', action="store_true", help="Flag for indicating whether we want to bootstrap the tree to determine node support (random sampling across all targets)")
    args = parser.parse_args()

    #Parse sample_info file
    sampleDict = import_sampleDict(args.sample_info)
    if not sampleDict:
        return

    #Import targetDict from pickle file
    targetDict = pickle.load(open(args.prefix + ".targetDict.pkl", "rb"))

    #Import alleleDict from pickle file
    alleleDict = pickle.load(open(args.prefix + ".alleleDict.pkl", "rb"))

    #Calculate "bulk" allelotype counts by merging all single cell calls
    alleleDict = calcBulk(alleleDict, targetDict)

    #Calculate original tre eusing all targets found within alleleDict (or subset as specified in args.target_file)
    if args.target_file is "All":
        target_list = sorted(alleleDict.keys())
    else:
        with open(target_file) as f:
            target_list = f.read().splitlines()
    #Calculate pairwise distance between each sample
    distDict_original = makeDistMatrix(sampleDict, alleleDict, target_list, args.dist_metric, False)
    #Draw neighbor joining tree
    tree_original = drawTree(distDict_original, sampleDict, args.outgroup, args.prefix, False) #We want to declare Fase for args.bootstrap because we want to output stats file for original tree
    print("----------Original Tree----------")
    print(tree_original.write(format = 0))

    #Bootstrap resample to create new distance matrices/trees and add support values to internal nodes of original tree
    if args.bootstrap is True:
        #Determine dictionary of tree nodes from tree_original
        nodeDict = {}
        for node in tree_original.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name) #We need to compare leaves in a node cluster without regard to order
            nodeDict[tuple(sorted(leaf_list))] = {}
            nodeDict[tuple(sorted(leaf_list))]["Bootstrap"] = 0 #Contains values for number of times random bootstrap tree (tree_temp) contains given node
            nodeDict[tuple(sorted(leaf_list))]["NodeID"] = node.write(format = 9)
        for i in tqdm(range(10)): #Bootstrap resample 10 times
            #Random downsample from pool of available targets (target_list) to use for distance calculation
            bootstrap_targets = list(np.random.choice(target_list, len(target_list), replace=True))
            distDict_temp = makeDistMatrix(sampleDict, alleleDict, bootstrap_targets, args.dist_metric, args.bootstrap)
            tree_temp = drawTree(distDict_temp, sampleDict, args.outgroup, args.prefix, args.bootstrap)
            nodeDict = bootstrapTree(nodeDict, tree_temp) #Determine whether each node in original tree is found in tree_temp
        #Add support information to original tree
        for node in tree_original.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name)
            node.add_feature(support = round(nodeDict[tuple(sorted(leaf_list))]["Bootstrap"]/10, 2))

    #Output tree (if bootstrapped, output with support values)
    print("----------Tree with Bootstrap Values----------")
    print(tree_original.write(format = 0))

if __name__ == "__main__":
    buildPhylo()
