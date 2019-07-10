#!/usr/bin/env python3
import argparse
import pickle
from statistics import mean
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

'''
python script.py --alleleDict ... --samples ... --prefix ...
We want to take as input alleleDict outputted by either Custom_allelotype or HipSTR_allelotype (or can be merged of both).  Also input a list of samples that we want to compare (must be >2 samples, 2 groups):
    group1  sample1
    group1  sample2
    ...
    group2  sample3
    group2  sample4
    ...
The structure of alleleDicts is as follows:
    alleleDict
        target_id (from targetDict)
            "sample"
                sample (from sampleDict, which is already defined when labeling readGroups prior to HipSTR)
                    "msCount"
                        list of msCounts
                    "allelotype"
                        list of alleles (2 alleles)
It will then plot the number of reads per target_id comparing the two groups.  This will be across average samples within a group (can have a '-ignoreEmpty' to not factor in samples with zero coverage)
The purpose for this script is to compare experiments from two different days or single cells belonging to two clones
'''

def comp_seqDepth(alleleDict, groupDict, ignoreEmpty, min_reads):
    depthDict = {}
    for target_id in sorted(alleleDict.keys()):
        depthDict[target_id] = {}
        for group_name in sorted(groupDict["sample_list"].keys()):
            depth_list = []
            for sample_name in sorted(groupDict["sample_list"][group_name]):
                if sample_name in alleleDict[target_id]["sample"].keys() and len(alleleDict[target_id]["sample"][sample_name]["msCount"]) >= min_reads:
                    depth_list.append(len(alleleDict[target_id]["sample"][sample_name]["msCount"]))
                else:
                    if ignoreEmpty is False:
                        depth_list.append(0)
            depthDict[target_id][group_name] = float(sum(depth_list) / max(len(depth_list), 1))
    #Save summary of seqDepth in groupDict
    groupDict["seqDepth"] = {}
    for group_name in sorted(groupDict["sample_list"].keys()):
        groupDict["seqDepth"][group_name] = []
    for target_id in sorted(depthDict.keys()):
        for group_name in sorted(groupDict["sample_list"].keys()):
            if ignoreEmpty is True:
                if len(depthDict[target_id].keys()) == len(groupDict["sample_list"].keys()): #If ignoreEmpty is true, we only want to save target_id seqDepth if exists in both group_names
                    groupDict["seqDepth"][group_name].append(depthDict[target_id][group_name])
            else:
                if group_name in depthDict[target_id].keys():
                    groupDict["seqDepth"][group_name].append(depthDict[target_id][group_name])
                else:
                    groupDict["seqDepth"][group_name].append(0)
    return groupDict

def plot_seqDepth(groupDict, prefix):
    group_list = sorted(list(groupDict["seqDepth"].keys()))
    seqDepth_df = pd.DataFrame(data = groupDict["seqDepth"])
    fig, ax = plt.subplots()
    # sns.set(font_scale=0.5)
    # seqDepth_plot = sns.jointplot(x = group_list[0], y = group_list[1], data = seqDepth_df, xlim=(0,100), ylim=(0,100))
    # seqDepth_plot.savefig(prefix + "." + ".".join(str(group_name) for group_name in group_list) + ".seqDepth.eps", format="eps", dpi=1000)

    fig, ax = plt.subplots()
    ax = sns.scatterplot(x = group_list[0], y = group_list[1], data = seqDepth_df)
    points = np.linspace(-10000, 10000, 100) #We want to also plot an identiy line for comparison
    plt.gca().plot(points, points, color='k', marker=None, linestyle='--', linewidth=1.0)
    ax.set_xlim([0,10000])
    ax.set_ylim([0,10000])
    plt.savefig(prefix + "." + ".".join(str(group_name) for group_name in group_list) + ".seqDepth.eps", format="eps", dpi=1000)
    return

def main():
    parser = argparse.ArgumentParser(description="Compare sequencing depth per target across sample groups")
    parser.add_argument('--alleleDict', action="store", dest="f_alleleDict", help="Merged alleleDict containing all samples we want to analyze")
    parser.add_argument('--samples', action="store", dest="f_samples", help="List of groups/sample names to be analyzed")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Prefix output plot files")
    parser.add_argument('-ignoreEmpty', action="store_true", default=False, help="Flag indicating we want to ")
    parser.add_argument('--min_reads', action="store", dest="min_reads", type=int, default=0, help="Minimum number of reads for target to be considered")
    args = parser.parse_args()

    alleleDict = pickle.load(open(args.f_alleleDict, 'rb'))

    #Import sample/group info
    groupDict = {}
    groupDict["sample_list"] = {}
    with open(args.f_samples, 'r') as f:
        for line in f:
            (group_name, sample_name) = line.split()
            if group_name not in groupDict["sample_list"].keys():
                groupDict["sample_list"][group_name] = []
            groupDict["sample_list"][group_name].append(sample_name)

    #We next want to create the two lists containing read depth for shared targets
    groupDict = comp_seqDepth(alleleDict, groupDict, args.ignoreEmpty, args.min_reads)

    #Plot sequencing depth for each target_id across the two groups
    plot_seqDepth(groupDict, args.prefix)

if __name__ == "__main__":
    main()
