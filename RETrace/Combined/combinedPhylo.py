#!/usr/bin/env python3
from RETrace.Combined.make_MSDist import make_MSDist
from RETrace.Combined.make_MethylDist import make_MethylDist
from RETrace.MS.utilities import import_targetDict
import gzip
import pickle
import os
import pandas as pd
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
from skbio import DistanceMatrix
from skbio.tree import nj

def importReg(Ensemble_gff): #High memory requirement (~40gb for hg19 Ensembl Reg Build)
    regDict = {}
    #We want to import the Ensembl Regulatory Build windows from the gff file
    with gzip.GzipFile(Ensemble_gff, 'rb') as f_gff:
        for line in f_gff:
            line = line.decode('utf-8')
            (chr, source, reg_type, start, end, score, strand, frame, attribute) = line.split("\t") #Based on gff format on Ensembl <https://uswest.ensembl.org/info/website/upload/gff.html>
            reg_name = reg_type + ':' + chr + ':' + start + '-' + end
            if chr not in regDict.keys():
                if chr.startswith('GL'): #We need to make a special condition for chromosomes that start with "GL" in order to match against sampleDict
                    chr = chr.replace('GL', 'gl').replace('.1', '')
                regDict[chr] = {}
            for pos in range(int(start), int(end) + 1): #We want to iterate through all base positions within reg window range
                regDict[chr][int(pos)] = reg_name
    return regDict

def drawTree(MS_distDict, Methyl_distDict, filtered_samples, ratio, outgroup):
    '''
    Merge MS and Methyl distance matrices
    '''
    merged_distMatrix = []
    for sample1 in sorted(filtered_samples):
        sample1_dist = []
        for sample2 in sorted(filtered_samples):
            merged_dist = (MS_distDict[sample1][sample2] * ratio) + (Methyl_distDict[sample1][sample2] * (1 - ratio)) / 100 #We want to scale methyl PD dist properly because PD is calculated from a 0-100 scale while MS dist is 0-1 scale
            sample1_dist.append(merged_dist)
        merged_distMatrix.append(sample1_dist)
    '''
    Run neighbor-joining phylogenetic tree building algorithm on pairwise cell distance (saved in distDict)
    '''
    distObj = DistanceMatrix(merged_distMatrix, sorted(filtered_samples))
    print(distObj.data)
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

def combinedPhylo(sample_list, target_info, alleleDict_file, methDict_file, Ensemble_gff, prefix, dist_metric, ratio, outgroup):
    '''
    This script will combine both MS and Methyl calls for each single cell in order to build a phylogenetic tree.  The way it does this is as follows:
        1) Utilize the alleleDict provided for each single cell within the sample_list calculate microsatellite distance between each cell
        2) Calculate pairwise dissimilarity between each pair of single cell (whether across all shared CpGs or those found in regulatory build windows)
        3) Using both MS and Methyl distance matrices, merge each pairwise distance for combined
    '''

    #Import the samples that we want to keep for combinedPhylo.  This is from sample_list file that contains a single sample_name per line
    filtered_samples = open(sample_list, 'r').read().splitlines()

    if os.path.isfile(prefix + ".MSDist.csv") is False:
        #Import necessary data for make_MSDist
        print("Importing MS Data")
        targetDict = import_targetDict(target_info)
        alleleDict = pickle.load(open(alleleDict_file, 'rb'))

        #Calculate pairwise MS distance
        print("Calculating MSDist")
        make_MSDist(filtered_samples, targetDict, alleleDict, dist_metric, prefix)
    else:
        print(prefix + ".MSDist.csv already computed")

    if os.path.isfile(prefix + ".MethylDist.csv") is False:
        #Import necessary data for make_MethylDist
        print("Importing Methyl Data")
        methDict = pickle.load(open(methDict_file, 'rb'))
        regDict = importReg(Ensemble_gff)

        #Calculate pairwise Methyl distance
        print("Calculating MethylDist")
        make_MethylDist(filtered_samples, methDict, regDict, prefix)
    else:
        print(prefix + ".MethylDist.csv already computed")

    #Merge MS and Methyl dist, then draw phylogenetic tree
    MS_distDict = pd.read_csv(prefix + ".MSDist.csv", index_col=0).to_dict()
    Methyl_distDict = pd.read_csv(prefix + ".MethylDist.csv", index_col=0).to_dict()
    tree_original = drawTree(MS_distDict, Methyl_distDict, filtered_samples, ratio, outgroup)
    f_tree_original = open(prefix + ".newick-original.txt", 'w')
    f_tree_original.write(tree_original.write(format = 0))
