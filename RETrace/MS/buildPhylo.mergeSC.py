#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict, import_targetDict
import pickle
import more_itertools
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree, TreeStyle #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
from tqdm import tqdm
import numpy as np
import random
import scipy
import math

def mergeSC(raw_alleleDict, targetDict, sampleDict, filtered_samples, n_merge):
    '''
    Randomly merge allelotype of n_merge number of singe cells that belong to the same clone
    We also want to cluster all alleles belonging to given target across all sample.  Use these as possible "allele_clusters"
    '''
    #We want to first determine the sample pairings that we want to use for the merge
    cloneDict = {} #Contains all clones/samples
    for sample in filtered_samples:
        clone = sampleDict[sample]["clone"]
        if sampleDict[sample]["clone"] in cloneDict:
            cloneDict[clone].append(sample)
        else:
            cloneDict[clone] = [sample]
    merge_list = []
    for clone in sorted(cloneDict.keys()):
        shuff_samples = cloneDict[clone]
        random.shuffle(shuff_samples)
        for indx in range(0, len(shuff_samples), n_merge):
            group = tuple(sorted(shuff_samples[indx:indx + n_merge]))
            merge_list.append(group)
    #Merge allelotypes from raw_alleleDict
    alleleDict = {}
    alleleDict["mergeSC"] = ['&'.join(group) for group in merge_list] #We want to save the list of merged SC
    for target_id in sorted(raw_alleleDict.keys()):
        if target_id in targetDict.keys(): #We only want to analye target_id specified in targetDict
            alleleDict[target_id] = {}
            alleleDict[target_id]["sample"] = {}
            alleleDict[target_id]["allele_groups"] = {}
            sub_len = len(targetDict[target_id]["sub_seq"]) #For calculating allele_clusters
            all_alleles = set() #For calculating allele_clusters
            for sample in sorted(raw_alleleDict[target_id]["sample"].keys()):
                for group in merge_list:
                    if sample in group:
                        group_name = '&'.join(group)
                        if group_name not in alleleDict[target_id]["sample"].keys():
                            alleleDict[target_id]["sample"][group_name] = {}
                            alleleDict[target_id]["sample"][group_name]["msCount"] = []
                            alleleDict[target_id]["sample"][group_name]["allelotype"] = []
                        alleleDict[target_id]["sample"][group_name]["msCount"].extend(raw_alleleDict[target_id]["sample"][sample]["msCount"])
                        alleleDict[target_id]["sample"][group_name]["allelotype"].extend(raw_alleleDict[target_id]["sample"][sample]["allelotype"])
                        try:
                            all_alleles.update([int(allele/sub_len) for allele in raw_alleleDict[target_id]["sample"][sample]["allelotype"]]) #Convert alleles to number of subunits (because we want to group all consecutive integer number subunit alleles)
                        except:
                            continue
            for indx, allele_cluster in enumerate(more_itertools.consecutive_groups(sorted(set(all_alleles)))): #Need to use "set(filtered_alleles)" in order to prevent repeated int in filtered_alleles
                alleleDict[target_id]["allele_groups"][indx] = [allele * sub_len for allele in list(allele_cluster)] #Save allele_cluster in terms of raw number of bases difference from ref
            # for group in alleleDict[target_id]["sample"].keys():
            #     print(target_id + "\t" + ','.join(group) + "\t" + ','.join(str(x) for x in alleleDict[target_id]["sample"][group]["allelotype"]))
    return alleleDict

def calcDist(alleleDict, distDict, sample_pair, sample1, sample2, shared_targets, dist_metric):
    '''
    Calculate distance between two samples given their shared_targets and dist_metric
    '''
    if sample1 != sample2:
        total_dist = 0
        num_comp = 0
        for target_id in shared_targets:
            for allele_indx, allele_group in alleleDict[target_id]["allele_groups"].items():
                dist_list = [] #We want to keep track of all possible distances between equal chunks of allelotype1 and allelotype2 and choose the minimum
                allelotype1 = list(set(alleleDict[target_id]["sample"][sample1]["allelotype"]).intersection(set(allele_group)))
                allelotype2 = list(set(alleleDict[target_id]["sample"][sample2]["allelotype"]).intersection(set(allele_group)))
                n = min(len(allelotype1), len(allelotype2)) #We want to look the number of alleles used for matching
                if n == 0:
                    continue #skip if at least one of the single cells don't have alleles found within allele_group
                for i in range(0, len(allelotype1), n):
                    temp_allelotype1 = allelotype1[i:i + n]
                    for j in range(0, len(allelotype2), n):
                        temp_allelotype2 = allelotype2[j:j + n]
                        if dist_metric == "EqorNot":
                            dist =  int(len(set(temp_allelotype1).symmetric_difference(set(temp_allelotype2))) / 2)
                        elif dist_metric == "Abs":
                            dist = sum(abs(x - y) for x, y in zip(sorted(temp_allelotype1), sorted(temp_allelotype2))) #This is based on <https://stackoverflow.com/questions/41229052/smallest-sum-of-difference-between-elements-in-two-lists>
                        dist_list.append(dist)
                # print(target_id + "\t" + sample1 + "\t" + sample2 + "\t" + str(allele_indx) + "\t" + ','.join(str(w) for w in allele_group) + "\t" + ','.join(str(x) for x in allelotype1) + "\t" + ','.join(str(y) for y in allelotype2) + "\t" + str(min(dist_list)) + "\t" + str(n))
                total_dist += min(dist_list)
                num_comp += n
    else:
        total_dist = 0
        num_comp = len(shared_targets)
    distDict["sampleComp"][sample_pair] = {}
    distDict["sampleComp"][sample_pair]["dist"] = float(total_dist/num_comp)
    distDict["sampleComp"][sample_pair]["num_targets"] = len(shared_targets)
    # distDict["sampleComp"][sample_pair]["num_alleles"] = int(num_comp)
    distDict["sampleComp"][sample_pair]["num_diff"] = int(total_dist)
    return distDict

def makeDistMatrix(sharedDict, alleleDict, mergeSC_list, dist_metric):
    '''
    Wrapper for calculating distance matrix (saved within distDict), which has the following structure:
        distDict
            "mergeSC"
            "sampleComp"
                sample1
                    sample2
                        "dist"
                        "num_targets"
    '''
    distDict = {}
    distDict["mergeSC"] = mergeSC_list #We want to save the names of each mergeSC into distDict
    distDict["sampleComp"] = {}
    #Calculate pairwise distance for all samples depending on shared targets (follow order specified in target_list [esp for bootstrapping, which may have duplicates due to sampling w/ replacement])
    for sample1 in tqdm(sorted(mergeSC_list)):
        for sample2 in sorted(mergeSC_list):
            sample_pair = tuple(sorted([sample1,sample2]))
            if sample_pair not in distDict["sampleComp"].keys():
                distDict = calcDist(alleleDict, distDict, sample_pair, sample1, sample2, sharedDict[sample_pair], dist_metric)
    return distDict

def drawTree(distDict, alleleDict, mergeSC_list, outgroup, prefix, bootstrap):
    '''
    Run neighbor-joining phylogenetic tree building algorithm on pairwise cell distance (saved in distDict)
    '''
    distMatrix = []
    targetMatrix = []
    pairwise_numTargets = []
    sample_numTargets = []
    for sample1 in sorted(mergeSC_list):
        sample1_dist = []
        sample1_targets = []
        for sample2 in sorted(mergeSC_list):
            sample_pair = tuple(sorted([sample1, sample2]))
            sample1_dist.append(distDict["sampleComp"][sample_pair]["dist"])
            sample1_targets.append(distDict["sampleComp"][sample_pair]["num_targets"])
        distMatrix.append(sample1_dist)
        targetMatrix.append(sample1_targets)
        if sample1 != sample2:
            pairwise_numTargets.append(distDict["sampleComp"][sample_pair]["num_targets"])
        else:
            sample_numTargets.append(distDict["sampleComp"][sample_pair]["num_targets"])
    if bootstrap is False: #Only output statistics for distance and number targets shared if for original tree (don't output for bootstrap resampling)
        statsOutput = open(prefix + ".buildPhylo.stats.txt", 'w')
        statsOutput.write("Avg targets shared per pair of cells:\t" + str(float(sum(pairwise_numTargets) / len(pairwise_numTargets))) + "\n")
        statsOutput.write("Avg targets captured per single cell:\t" + str(float(sum(sample_numTargets) / len(sample_numTargets))) + "\n")
        for dist_indx,dist_list in enumerate(distMatrix): #Print matrix containing distances
            statsOutput.write(sorted(mergeSC_list)[dist_indx] + "," + ",".join(str(round(i,3)) for i in dist_list) + "\n")
        for target_indx,target_list in enumerate(targetMatrix): #Print matrix containing number targets shared between each pair
            statsOutput.write(sorted(mergeSC_list)[target_indx] + "," + ",".join(str(j) for j in target_list) + "\n")
        statsOutput.close()
        pickle.dump(distDict, open(prefix + ".buildPhylo.distDict.pkl", "wb")) #We want to print out the distance information for each single cell pair that was used to buildPhylo (this will be useful for downstream statistics)
    distObj = DistanceMatrix(distMatrix,sorted(mergeSC_list))
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

def bootstrapTree(nodeDict, treeTemp, bootstrap_samples):
    '''
    Evaluate whether original node (keys of nodeDict) are found within treeTemp (tree from bootstrapped samples).
    Be sure to only consider sample_name within bootstrap_samples, which may be different than original because of bootstrap resampling
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
        node_intersect = tuple(sorted(set(node).intersection(bootstrap_samples)))
        if node_intersect == node:
            nodeDict[node]["Num_sampled"] += 1
            if node_intersect in tempNodes:
                nodeDict[node]["Num_verified"] += 1
    return nodeDict

def buildPhylo(sample_info, sample_list, prefix, target_info, alleleDict_file, dist_metric, outgroup, bootstrap, n_merge):
    '''
    Draw phylogenetic tree based on calculated allelotype of single cells.  This is done by:
        1) Filtering out only likely alleles by creating a "pseudo"-bulk in which we cluster all single cell alleles together (calcBulk)
        2) Make distance matrix usin given dist_metric (Abs or EqorNot)
        3) Draw tree using neighbor-joining method
        4) Repeating 2-3 but bootstrapping the samples in order to determine the support of each node of the original tree
    '''

    #Import the samples that we want to keep for buildPhylo
    sampleDict = import_sampleDict(sample_info)
    filtered_samples = open(sample_list, 'r').read().splitlines()

    targetDict = import_targetDict(target_info)

    #Calculate "bulk" allelotype counts by merging all single cell calls
    raw_alleleDict = pickle.load(open(alleleDict_file, 'rb'))

    #We want to merge n number of single cells together
    print("\tMerging SC data")
    alleleDict = mergeSC(raw_alleleDict, targetDict, sampleDict, filtered_samples, n_merge)
    mergeSC_list = alleleDict["mergeSC"]

    #Pre-calculate shared targets between each sample
    print("\tPre-computing shared target_id between each pairwise sample")
    sharedDict = {} #Contains all shared targets between each pairwise sample
    for sample1 in tqdm(mergeSC_list):
        for sample2 in sorted(mergeSC_list):
            sample_pair = tuple(sorted([sample1,sample2]))
            if sample_pair not in sharedDict.keys():
                shared_targets = []
                for target_id in sorted(alleleDict.keys()):
                    if target_id == "mergeSC":
                        continue
                    if sample1 in alleleDict[target_id]["sample"].keys() and sample2 in alleleDict[target_id]["sample"].keys() and target_id in targetDict.keys(): #Only want to analyze target_id in targetDict
                        if "allelotype" in alleleDict[target_id]["sample"][sample1].keys() and "allelotype" in alleleDict[target_id]["sample"][sample2].keys():
                            shared_targets.append(target_id)
                sharedDict[sample_pair] = shared_targets
                # print("\t".join(sample_pair) + "\t" + str(len(shared_targets)))

    #Calculate original tree using all samples found within sampleDict
    print("\tCalculating distance matrix for each pairwise comparison and drawing original Newick tree (no bootstrap)")
    distDict_original = makeDistMatrix(sharedDict, alleleDict, mergeSC_list, dist_metric) #Calculate pairwise distance between each sample
    MStree = drawTree(distDict_original, alleleDict, mergeSC_list, outgroup, prefix, False) #Draw neighbor joining tree.  We want to declare Fase for bootstrap because we want to output stats file for original tree
    f_MStree = open(prefix + '.buildPhylo.newick-original.txt', 'w')
    f_MStree.write(MStree.write(format = 0))
    f_MStree.write("\n")
    f_MStree.close()

    #Bootstrap resample to create new distance matrices/trees and add support values to internal nodes of original tree
    if bootstrap is True:
        #Determine dictionary of tree nodes from MStree
        nodeDict = {}
        for node in MStree.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name) #We need to compare leaves in a node cluster without regard to order
            nodeDict[tuple(sorted(leaf_list))] = {}
            nodeDict[tuple(sorted(leaf_list))]["Num_verified"] = 0 #Contains values for number of times random bootstrap tree (tree_temp) contains given node
            nodeDict[tuple(sorted(leaf_list))]["Num_sampled"] = 0 #Contains number of times the node occured during bootstrap resampling
            # nodeDict[tuple(sorted(leaf_list))]["NodeID"] = node.write(format = 9)
        for i in tqdm(range(10000)): #Bootstrap resample 10,000 times
            #Random downsample from pool of available targets (target_list) to use for distance calculation
            bootstrap_samples = set(np.random.choice(filtered_samples, len(filtered_samples), replace=True))
            if outgroup not in ["Midpoint", "NA"]: #We want to make sure our outgroup remains in the tree even during bootstraping
                bootstrap_samples.add(outgroup)
            tree_temp = drawTree(distDict_original, alleleDict, bootstrap_samples, outgroup, prefix, bootstrap)
            nodeDict = bootstrapTree(nodeDict, tree_temp, bootstrap_samples) #Determine whether each node in original tree is found in tree_temp
        #Add support information to original tree
        for node in MStree.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name)
            if nodeDict[tuple(sorted(leaf_list))]["Num_sampled"] > 0:
                node_support = round(float(nodeDict[tuple(sorted(leaf_list))]["Num_verified"] / nodeDict[tuple(sorted(leaf_list))]["Num_sampled"]),2)
            else:
                node_support = 0.00 #Assign nodes that were not present in any bootstrap simulation a value of 2.0
            node.add_features(support = node_support)
        #Output tree with optional support values
        f_tree_bootstrap = open(prefix + '.buildPhylo.newick-bootstrap.txt', 'w')
        f_tree_bootstrap.write(MStree.write(format = 0))
        f_tree_bootstrap.write("\n")
        f_tree_bootstrap.close()
