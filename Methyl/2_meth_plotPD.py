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

def plotPD(PDdict_pkl, prefix):
    with open(PDdict_pkl, 'rb') as PDdict_file:
        PDdict = pickle.load(PDdict_file)

    PD_df = pandas.DataFrame(PDdict["PD"], index=sorted(PDdict["index"]))

    sns.set(font_scale=2)
    PD_clustermap = sns.clustermap(PD_df)
    PD_clustermap.savefig(prefix + ".PD.png")

    PD_clustermap_z = sns.clustermap(PDdf, z_score=1) #Draw clustermap based on z-scores of each column (i.e. sample input) with lower z-score indicating higher similarity (lower PD)
    PD_clustermap_z.savefig(prefix + ".z_norm.png")
    return

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--PDdict', action="store", dest="PDdict_pkl", help="Pre-computed PDdict containing pairwise dissimilarity between samples and cell types")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specifies prefix for output files")
    args = parser.parse_args()

    plotPD(args.PDdict_pkl, args.prefix)

#%%
if __name__ == "__main__":
    main()
