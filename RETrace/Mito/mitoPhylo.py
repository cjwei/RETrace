#!/usr/bin/env python3
import argparse
import re
from tqdm import tqdm
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap
import math
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
import pickle

def import_sampleDict(info_file):
    sampleDict = {}
    with open(info_file) as f:
        for line in f:
            if len(line.split()) == 4: #If clone is specified
                (pileup, sample, sex, clone) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["pileup"] = pileup
                sampleDict[sample]["sex"] = sex
                sampleDict[sample]["clone"] = clone
            elif len(line.split()) == 3: #If clone is not specified
                (pileup, sample, sex) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["pileup"] = pileup
                sampleDict[sample]["sex"] = sex
            else:
                print("Incorrect formatting for line (pileup, sample, sex, [optional] clone):\n" + "\t".join(line.split()))
    return sampleDict

def import_pileup(sampleDict):
    '''
    This script will take as input sampleDict and collect base information for each sample
    '''
    baseDict = {}
    for sample in tqdm(sorted(sampleDict.keys())):
        with open(sampleDict[sample]["pileup"]) as f_pileup:
            for line in f_pileup:
                (chrom, loc, ref_base, num_reads, bases, quals) = line.split()
                base_loc = chrom + ":" + loc
                if base_loc not in baseDict.keys():
                    baseDict[base_loc] = {}
                baseDict[base_loc][sample] = [0, 0, 0, 0] #We want to keep a tuple containing the number of reads for each possible base (A,C,G,T)
                base_id = ("A", "C", "G", "T")
                #We want to remove additional characters (this is similar to <https://github.com/riverlee/pileup2base/blob/master/pileup2baseindel_no_strand.pl>)
                #1) remove ^. pattern
                bases = re.sub(r"\^.", "", bases)
                #2) remove $ pattern
                bases = re.sub(r"\$.", "", bases)
                #3) remove [-+][0-9]+[ACGTNacgtn]+ pattern (indel)
                bases = re.sub(r"[-+](\d+)(\w+)", "", bases)
                for read in list(bases):
                    if read in ('.', ','):
                        baseDict[base_loc][sample][base_id.index(ref_base.upper())] += 1
                    elif read.upper() == 'A':
                        baseDict[base_loc][sample][0] += 1
                    elif read.upper() == 'C':
                        baseDict[base_loc][sample][1] += 1
                    elif read.upper() == 'G':
                        baseDict[base_loc][sample][2] += 1
                    elif read.upper() == 'T':
                        baseDict[base_loc][sample][3] += 1
    return baseDict

def filterBases(sampleDict, baseDict, min_cov, heteroplasmy):
    #We want to first filter out only bases that have a minimum of 2.5% heteroplasmy in at least one sample
    filtered_baseDict = {}
    for base_loc in sorted(baseDict.keys()):
        for indx in range(4):
            num_pass = 0
            for sample in sorted(baseDict[base_loc].keys()):
                if (sum(baseDict[base_loc][sample]) >= min_cov) and (baseDict[base_loc][sample][indx] / sum(baseDict[base_loc][sample])) >= heteroplasmy:
                    base_info = tuple([base_loc, indx])
                    if base_info in filtered_baseDict.keys():
                        filtered_baseDict[base_info].append(sample)
                    else:
                        filtered_baseDict[base_info] = [sample]
    return filtered_baseDict

def plot_baseUMAP(sampleDict, baseDict, filtered_baseDict, base_prefix):
    '''
    We want to plot the basecalls for each single cell using UMAP.
    '''
    baseMatrix = []
    for base_info in tqdm(sorted(filtered_baseDict.keys())):
        (base_loc, indx) = base_info
        base_freq = []
        for sample in sorted(sampleDict.keys()):
            if sample in filtered_baseDict[base_info]:
                base_freq.append(baseDict[base_loc][sample][indx] / sum(baseDict[base_loc][sample]))
            else:
                base_freq.append(None)
        if None not in base_freq:
            baseMatrix.append(base_freq)
    baseArray = np.array(baseMatrix)
    print(baseArray.shape)
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(baseArray.transpose())

    #Create UMAP plot for each clone.  In order to do so, we want to group UMAP embedding array with each respective clone
    fig, ax = plt.subplots()
    clone_set = set()
    for sample in sorted(sampleDict.keys()):
        clone_set.add(sampleDict[sample]["clone"])
    color_indx = 0
    for clone in sorted(clone_set):
        sub_matrix = []
        for indx in range(len(sampleDict.keys())):
            sample = list(sorted(sampleDict.keys()))[indx]
            if sampleDict[sample]["clone"] == clone:
                sub_matrix.append([embedding[indx, 0], embedding[indx, 1]])
        sub_embedding = np.array(sub_matrix)
        ax.scatter(sub_embedding[:, 0], sub_embedding[:, 1], c=sns.color_palette()[color_indx], label=clone)
        color_indx += 1
    plt.legend(loc = 'best')
    plt.gca().set_aspect('equal', 'datalim')
    plt.savefig(base_prefix + ".UMAP.png")
    return

def plot_distUMAP(sampleDict, baseDict, filtered_baseDict, dist_prefix):
    '''
    We want to plot a UMAP using distance matrix
    '''
    #We first want to calculate a distDict for each pairwise sample comparison
    distDict = {}
    distDict["mergeSC"] = sorted(sampleDict.keys()) #We want to specify the mergeSC to make this compatible to legacy evalPhylo
    distDict["sampleComp"] = {}
    num_bases_list = []
    for sample1 in tqdm(sorted(sampleDict.keys())):
        for sample2 in sorted(sampleDict.keys()):
            sample_pair = tuple(sorted([sample1, sample2]))
            if sample_pair not in distDict["sampleComp"].keys():
                distDict["sampleComp"][sample_pair] = {}
                dist = 0
                num_bases = 0
                for base_info in sorted(filtered_baseDict.keys()):
                    if sample1 in filtered_baseDict[base_info] and sample2 in filtered_baseDict[base_info]:
                        (base_loc, indx) = base_info
                        num_reads1 = sum(baseDict[base_loc][sample1])
                        num_reads2 = sum(baseDict[base_loc][sample2])
                        dist += math.sqrt(abs(baseDict[base_loc][sample1][indx] / num_reads1 - abs(baseDict[base_loc][sample2][indx] / num_reads2)))
                        num_bases += 1
                distDict["sampleComp"][sample_pair]["dist"] = float(dist / num_bases)
                distDict["sampleComp"][sample_pair]["num_bases"] = num_bases
                num_bases_list.append(num_bases)
    print("Average bases shared between each pair of SC:\t" + str(sum(num_bases_list) / len(num_bases_list)))
    #We want to create a distMatrix summarizing the pairwise distances
    distMatrix = []
    for sample1 in sorted(sampleDict.keys()):
        sample1_dist = []
        for sample2 in sorted(sampleDict.keys()):
            sample_pair = tuple(sorted([sample1, sample2]))
            sample1_dist.append(distDict["sampleComp"][sample_pair]["dist"])
        distMatrix.append(sample1_dist)
    distObj = DistanceMatrix(distMatrix, sorted(sampleDict.keys()))
    skbio_tree = nj(distObj, result_constructor=str)
    ete_tree = Tree(skbio_tree)
    tree_midpoint = ete_tree.get_midpoint_outgroup()
    ete_tree.set_outgroup(tree_midpoint)
    f_treeOut = open(dist_prefix + ".newick.txt", 'w')
    f_treeOut.write(ete_tree.write(format = 0))
    #We want to save the distDict for use in evalPhylo
    pickle.dump(distDict, open(dist_prefix + ".mitoPhylo.pkl", "wb"))

    #We also want to draw UMAP using the pairwise distance matrix
    distArray = np.array(distMatrix)
    print(distArray.shape)
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(distArray.transpose())

    #Create UMAP plot for each clone.  In order to do so, we want to group UMAP embedding array with each respective clone
    fig, ax = plt.subplots()
    clone_set = set()
    for sample in sorted(sampleDict.keys()):
        clone_set.add(sampleDict[sample]["clone"])
    color_indx = 0
    for clone in sorted(clone_set):
        sub_matrix = []
        for indx in range(len(sampleDict.keys())):
            sample = list(sorted(sampleDict.keys()))[indx]
            if sampleDict[sample]["clone"] == clone:
                sub_matrix.append([embedding[indx, 0], embedding[indx, 1]])
        sub_embedding = np.array(sub_matrix)
        ax.scatter(sub_embedding[:, 0], sub_embedding[:, 1], c=sns.color_palette()[color_indx], label=clone)
        color_indx += 1
    plt.legend(loc = 'best')
    plt.gca().set_aspect('equal', 'datalim')
    plt.savefig(dist_prefix + ".UMAP.png")

    return

def main():
    parser = argparse.ArgumentParser(description="Parse pileup file to obtain ratio of bases for each mitochondrial DNA position")
    parser.add_argument("--sample_info", action="store", dest="sample_info", help="Sample info file containing SC pileupe locations and and clone information")
    parser.add_argument("--base_prefix", action="store", dest="base_prefix", help="Specify UMAP basecall plot file")
    parser.add_argument("--dist_prefix", action="store", dest="dist_prefix", help="Specify UMAP distance metric plot file")
    parser.add_argument("--heteroplasmy", action="store", type=float, default=0.025, dest="heteroplasmy", help="Specify the minimum percentage of heteroplasmy for at least one cell (defaul [from Ludwig et al] = 0.025)")
    parser.add_argument("--min_cov", action="store", type=int, default=100, dest="min_cov", help="Specify miminum number of read coverage for given base in single cell")
    args = parser.parse_args()

    sampleDict = import_sampleDict(args.sample_info)

    baseDict = import_pileup(sampleDict)

    #We want to filter the baseDict according to heteroplasmy and min_cov
    filtered_baseDict = filterBases(sampleDict, baseDict, args.min_cov, args.heteroplasmy)

    plot_baseUMAP(sampleDict, baseDict, filtered_baseDict, args.base_prefix)

    plot_distUMAP(sampleDict, baseDict, filtered_baseDict, args.dist_prefix)

if __name__ == "__main__":
    main()
