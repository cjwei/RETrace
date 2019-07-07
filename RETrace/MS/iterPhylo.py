#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict, import_targetDict
import pickle
import more_itertools
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
from tqdm import tqdm
import numpy as np
import random
import scipy
import math
import operator

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

def calcDist(mergeDict, alleleDict, sample1, sample2, dist_metric):
    '''
    Calculate distance between two samples given desired dist_metric
    '''
    total_dist = 0
    num_comp = 0
    for target_id in sorted(alleleDict.keys()):
        sample1_list = set(mergeDict["leafDict"][sample1]).intersection(set(alleleDict[target_id]["sample"]))
        sample2_list = set(mergeDict["leafDict"][sample2]).intersection(set(alleleDict[target_id]["sample"]))
        if len(sample1_list) > 0 and len(sample2_list) > 0: #We want to analyze target_id that was covered by both sample1 and sample2
            for allele_indx, allele_group in alleleDict[target_id]["allele_groups"].items():
                dist_list = [] #We want to keep track of all possible distances between equal chunks of allelotype1 and allelotype2 and choose the minimum
                allelotype1 = []
                for sample in sample1_list:
                    allelotype1.extend(list(set(alleleDict[target_id]["sample"][sample]["allelotype"]).intersection(set(allele_group))))
                allelotype2 = []
                for sample in sample2_list:
                    allelotype2.extend(list(set(alleleDict[target_id]["sample"][sample]["allelotype"]).intersection(set(allele_group))))
                for allele1 in sorted(allelotype1):
                    for allele2 in sorted(allelotype2):
                        if dist_metric == "EqorNot":
                            if allele1 != allele2:
                                total_dist += 1
                        elif dist_metric == "Abs":
                            total_dist += abs(allele1 - allele2)
                        num_comp += 1
                # if '&' in sample1 or '&' in sample2:
                #     print(sample1 + ":" + ','.join(str(x) for x in sorted(allelotype1)) + "\t" + sample2 + ":" + ','.join(str(y) for y in sorted(allelotype2)))
                # n = min(len(allelotype1), len(allelotype2)) #We want to look the number of alleles used for matching
                # if n == 0:
                #     continue #skip if at least one of the single cells don't have alleles found within allele_group
                # for i in range(0, len(allelotype1), n):
                #     temp_allelotype1 = allelotype1[i:i + n]
                #     for j in range(0, len(allelotype2), n):
                #         temp_allelotype2 = allelotype2[j:j + n]
                #         if dist_metric == "EqorNot":
                #             dist =  int(len(set(temp_allelotype1).symmetric_difference(set(temp_allelotype2))) / 2)
                #         elif dist_metric == "Abs":
                #             dist = sum(abs(x - y) for x, y in zip(sorted(temp_allelotype1), sorted(temp_allelotype2))) #This is based on <https://stackoverflow.com/questions/41229052/smallest-sum-of-difference-between-elements-in-two-lists>
                #         dist_list.append(dist)
                # # print(target_id + "\t" + sample1 + "\t" + sample2 + "\t" + str(allele_indx) + "\t" + ','.join(str(w) for w in allele_group) + "\t" + ','.join(str(x) for x in allelotype1) + "\t" + ','.join(str(y) for y in allelotype2) + "\t" + str(min(dist_list)) + "\t" + str(n))
                # total_dist += min(dist_list)
                # num_comp += n
    #We want to save this into distDict
    mergeDict["distDict"]['&'.join(sorted([sample1, sample2]))] = float(total_dist / num_comp)
    return mergeDict

def mergeClosest(mergeDict, alleleDict, dist_metric):
    min_dist = float("inf")
    for sample1 in mergeDict["availNodes"]:
        for sample2 in mergeDict["availNodes"]:
            if sample1 != sample2:
                potential_node = '&'.join(sorted([sample1, sample2]))
                if potential_node not in mergeDict["distDict"].keys():
                    mergeDict = calcDist(mergeDict, alleleDict, sample1, sample2, dist_metric)
                if mergeDict["distDict"][potential_node] <= min_dist:
                    # print(potential_node + "\t" + str(mergeDict["distDict"][potential_node]))
                    min_dist = mergeDict["distDict"][potential_node]
                    mergeNode = potential_node
                    mergeChildren = sorted([sample1, sample2])
                    mergeLeaves = mergeDict["leafDict"][sample1] | mergeDict["leafDict"][sample2]
    #Update availNodes by removing children and replacing with the final merged node
    for child in mergeChildren:
        mergeDict["availNodes"].remove(child)
    mergeDict["availNodes"].append(mergeNode)
    #Update leafDict with the leaves for the merged node
    mergeDict["leafDict"][mergeNode] = mergeLeaves
    #Update nodeDict with the children of hte merged node
    mergeDict["nodeDict"][mergeNode] = mergeChildren
    # print("----------Chosen node:----------\n" + mergeNode + "\n" + ','.join(mergeChildren))
    return mergeDict

def iterPhylo(sample_info, sample_list, prefix, target_info, alleleDict_file, dist_metric):
    '''
    Iteratively build phylogenetic tree based on calculated allelotype of single cells.  This is different than buildPhylo because we are iteratively merging single cells using the following steps:
        1) Filtering out only likely alleles by creating a "pseudo"-bulk in which we cluster all single cell alleles together (calcBulk)
        2) Make distance matrix usin given dist_metric (Abs or EqorNot)
        3) Determine the pair of single cells (or merged cells) that are closest to one another and merge into new merged
        4) Repeat steps 2 & 3 until reach convergence (number of merged cells is 1)
    '''

    #Import the samples that we want to keep for buildPhylo
    sampleDict = import_sampleDict(sample_info)
    filtered_samples = open(sample_list, 'r').read().splitlines()

    targetDict = import_targetDict(target_info)

    #Calculate "bulk" allelotype counts by merging all single cell calls
    alleleDict = pickle.load(open(alleleDict_file, 'rb'))
    alleleDict = calcBulk(alleleDict, targetDict)

    #We want to iterateively merge single cell together
    mergeDict = {}
    mergeDict["nodeDict"] = {} #This contains the final tree relationships between parent and two children
    mergeDict["leafDict"] = {} #This contains all leaves of the parents found in nodeDict
    mergeDict["distDict"] = {} #Initialize dictionary containing all distances between possible nodes
    mergeDict["availNodes"] = [] #This contains all possible available nodes for merging
    for sample in sorted(filtered_samples): #We want to initialize the mergeDict and availNodes list to contain all leaves of the tree (the individual samples in filtered_samples)
        mergeDict["availNodes"].append(sample)
        mergeDict["leafDict"][sample] = set([sample])
    while len(mergeDict["availNodes"]) > 1: #We want to continue th eanalysis until all availNodes have been collapsed into a single node
        mergeDict = mergeClosest(mergeDict, alleleDict, dist_metric)
