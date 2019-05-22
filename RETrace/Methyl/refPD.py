#!/usr/bin/env python3
import os
import pickle
from tqdm import tqdm
import pickle
import pandas
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

def calcPD(sampleDict, typeDict, filtered_samples, DMR_bool, min_shared, min_rate, prefix):
    #We want to iterate through all cell types and calculate pairwise dissimilarity when compared to all samples
    PDdict_temp = {} #Keep all pairwise dissimilarity comparisons in temp dictionary
    for cellType_name in sorted(typeDict["cellType"].keys()):
        print("Intersecting:\t" + cellType_name)
        f_cellType = open(typeDict["cellType"][cellType_name]["file_loc"])
        for i, line in enumerate(tqdm(f_cellType)):
            if line.startswith('#'):
                continue
            base_info = line.split()
            (chr, pos, chain, ctype, meth, total_reads) = base_info[0:6]
            if ctype.startswith('CG'): #We only care about CGN methylated cytosine calls
                base_loc = chr.replace("chr","") + ":" + pos + ";" + chain
                cellType_methRate = float(int(meth)/int(total_reads))
            else:
                continue
            if DMR_bool is True and 'DMR' not in line:
                continue #Skip base_loc if not contained within DMR and DMR_bool is set to True
            if not min(min_rate, 1 - min_rate) < cellType_methRate < max(min_rate, 1 - min_rate): #Only consider base locations in which reference methRate is not between (min_rate, 1-min_rate)
                if base_loc in sampleDict["base"].keys(): #We next want to search matching base_loc in sampleDict
                    for sample_name in sampleDict["base"][base_loc]:
                        if sample_name in filtered_samples: #Only analyze samples of interest
                            #We want to calculate pairwise dissimilarity between sample and cellType and save into PDdict_temp
                            sample_methRate = float(sampleDict["base"][base_loc][sample_name][0] / sampleDict["base"][base_loc][sample_name][1])
                            dis = abs(sample_methRate - cellType_methRate) * 100
                            if sample_name not in PDdict_temp.keys():
                                PDdict_temp[sample_name] = {}
                            if cellType_name not in PDdict_temp[sample_name].keys():
                                PDdict_temp[sample_name][cellType_name] = []
                            PDdict_temp[sample_name][cellType_name].append(dis)
        f_cellType.close()
    #Save pairwise dissimilarity into PDdict for future processing/plotting
    cellType_list = []
    PDdict = {}
    PDdict["PD"] = {}
    PD_output = open(prefix + ".PD.txt", 'w')
    PD_output.write("Sample\tCell Type\tPairwise Dissimilarity\tNum CpG Shared\n")
    for sample_name in sorted(PDdict_temp.keys()):
        PD_list = []
        for cellType_name in sorted(PDdict_temp[sample_name].keys()):
            if len(PDdict_temp[sample_name][cellType_name]) >= min_shared:
                pairwise_dist = float(sum(PDdict_temp[sample_name][cellType_name])/len(PDdict_temp[sample_name][cellType_name]))
                PD_output.write(sample_name + "\t" + cellType_name + "\t" + str(round(pairwise_dist, 4)) + "\t" + str(len(PDdict_temp[sample_name][cellType_name])) + "\n")
                PD_list.append(pairwise_dist)
                cellType_list.append(cellType_name)
            else:
                PD_list.append(0)
        PDdict["PD"][sample_name] = PD_list
    PD_output.close()
    PDdict["index"] = sorted(list(set(cellType_list))) #This contains cell type names used for pandas dataframe index
    #Export PDdict to file using pickle
    with open(prefix + ".PD.pkl", 'wb') as PDdict_file:
        pickle.dump(PDdict, PDdict_file, protocol=pickle.HIGHEST_PROTOCOL)
    return PDdict

def plotPD(PDdict, prefix):
    #We want to ensure that the dataframe has the same number of columns (remove sample_name that don't have PD calculated for all reference cell types)
    for sample_name in sorted(PDdict["PD"].keys()):
        if len(PDdict["PD"][sample_name]) < len(PDdict["index"]):
            # print(sample_name + "\t" + ','.join(str(i) for i in PDdict["PD"][sample_name]))
            PDdict["PD"].pop(sample_name, None)
    PD_df = pandas.DataFrame(PDdict["PD"], index=sorted(PDdict["index"]))

    sns.set(font_scale=0.5)
    print("Drawing clustermap for raw PD")
    PD_clustermap = sns.clustermap(PD_df, xticklabels=True)
    PD_clustermap.savefig(prefix + ".PD.eps", format="eps", dpi=1000)

    print("Drawing clustermap with Z-score across samples")
    PD_clustermap_z_sample = sns.clustermap(PD_df, z_score=1, xticklabels=True) #Draw clustermap based on z-scores of each column (i.e. sample input) with lower z-score indicating higher similarity (lower PD)
    PD_clustermap_z_sample.savefig(prefix + ".PD.z_sample.eps", format="eps", dpi=1000)

    print("Drawing clustermap for Z-score across cell types")
    PD_clustermap_z_cellType = sns.clustermap(PD_df, z_score=0, xticklabels=True) #Draw clustermap based on z-scores of each row (i.e. cell type)
    PD_clustermap_z_cellType.savefig(prefix + ".PD.z_cellType.eps", format="eps", dpi=1000)
    return

def refPD(sample_list, ref_info, sample_methDict, prefix, DMR, min_shared, min_rate):
    '''
    This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039> and Paiwise Dissimilarity calculations from Hui 2018 Stem Cell Reports paper <https://doi.org/10.1016/j.stemcr.2018.07.003>.
    It will take as input methylation info for single cell samples and reference cellt ypes in order to calculate the pairwise dissimilarity between each
    '''
    if os.path.exists(prefix + ".PD.pkl") is False:
        #Import the samples that we want to keep for PD analysis
        filtered_samples = open(sample_list, 'r').read().splitlines()

        #Import pre-computed sample methDict
        print("Import sampleDict")
        with open(sample_methDict, 'rb') as sample_methDict_file:
            sampleDict = pickle.load(sample_methDict_file)

        #Import cell type file information into cellDict
        print("Import reference cell type info")
        typeDict = {}
        typeDict["cellType"] = {}
        with open(ref_info) as f_info:
            for cellType in f_info:
                (file_loc, cellType_name) = cellType.split()
                typeDict["cellType"][cellType_name] = {}
                typeDict["cellType"][cellType_name]["file_loc"] = file_loc

        #Calculate pairwise dissimilarity matrix
        print("Calculating PD between single cells and cell type")
        PDdict = calcPD(sampleDict, typeDict, filtered_samples, DMR, int(min_shared), float(min_rate), prefix)
    else:
        print("Importing pre-computed PDdict")
        with open(prefix + ".PD.pkl", 'rb') as PDdict_f:
            PDdict = pickle.load(PDdict_f)

    #Plot PD
    print("Plotting PD between single cells and cell type")
    plotPD(PDdict, prefix)
