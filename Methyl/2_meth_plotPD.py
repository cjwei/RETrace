#!/usr/bin/env python3
import argparse
import pickle
import pandas
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import numpy as np

'''
Usage: python script.py --PDdict ... --prefix
This script will go through and plot the PDdict file comparing single cell samples to cell type.  It will output heatmap plots summarizing the PD calculated
'''

def plotPD(PDdict_pkl, labelDict, prefix):
    with open(PDdict_pkl, 'rb') as PDdict_file:
        PDdict = pickle.load(PDdict_file)

    for sample_name in sorted(PDdict["PD"].keys()):
        #We want to check whether we specify the samples to include in histogram
        if any(labelDict):
            if sample_name not in labelDict.keys():
                PDdict["PD"].pop(sample_name, None)
                continue
            else:
                label = labelDict[sample_name]
                PDdict["PD"][label] = PDdict["PD"].pop(sample_name)
                sample_name = label
        #We want to ensure that the dataframe has the same number of columns (remove sample_name that don't have PD calculated for all reference cell types)
        if len(PDdict["PD"][sample_name]) < len(PDdict["index"]):
            print(sample_name + "\t" + ','.join(str(i) for i in PDdict["PD"][sample_name]))
            PDdict["PD"].pop(sample_name, None)
    PD_df = pandas.DataFrame(PDdict["PD"], index=sorted(PDdict["index"]))

    sns.set(font_scale=0.5)
    PD_clustermap = sns.clustermap(PD_df, xticklabels=True)
    PD_clustermap.savefig(prefix + ".PD.eps", format="eps", dpi=1000)

    PD_clustermap_z_sample = sns.clustermap(PD_df, z_score=1, xticklabels=True) #Draw clustermap based on z-scores of each column (i.e. sample input) with lower z-score indicating higher similarity (lower PD)
    PD_clustermap_z_sample.savefig(prefix + ".z_sample.eps", format="eps", dpi=1000)

    PD_clustermap_z_cellType = sns.clustermap(PD_df, z_score=0, xticklabels=True) #Draw clustermap based on z-scores of each row (i.e. cell type)
    PD_clustermap_z_cellType.savefig(prefix + ".z_cellType.eps", format="eps", dpi=1000)
    return

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--PDdict', action="store", dest="PDdict_pkl", help="Pre-computed PDdict containing pairwise dissimilarity between samples and cell types")
    parser.add_argument('--samples', action="store", dest="sample_file", default=None, help="[Optional] Specify the samples that we want to include in heatmap along with new labels (i.e. all single cells, etc)")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specifies prefix for output files")
    args = parser.parse_args()

    #We want to import the new sample labels for plotting
    labelDict = {}
    if args.sample_file is not None:
        with open(args.sample_file) as f:
            for line in f:
                [sample_name, label] = line.split()
                labelDict[sample_name] = label

    plotPD(args.PDdict_pkl, labelDict, args.prefix)

#%%
if __name__ == "__main__":
    main()
