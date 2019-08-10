#!/usr/bin/env python3
import os
import pickle
from tqdm import tqdm
import pickle
import pandas
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from statistics import mean
import numpy as np
import re

def calcPD(sampleDict, typeDict, filtered_samples, DMR_bool, reg_bool, min_shared, min_rate, min_reads, prefix):
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
            region_set = set()
            if DMR_bool is True:
                if 'DMR' not in line:
                    continue #Skip base_loc if not contained within DMR and DMR_bool is set to True
                else:
                    region_set.update([region for region in base_info[7:] if 'DMR' in region])
            if reg_bool is True:
                if 'reg' not in line:
                    continue #Skip base_loc if not contained within regulatory build window and reg_bool is set to True
                else:
                    region_set.update([region for region in base_info[7:] if 'reg' in region])
            if not min(min_rate, 1 - min_rate) < cellType_methRate < max(min_rate, 1 - min_rate): #Only consider base locations in which reference methRate is not between (min_rate, 1-min_rate)
                if base_loc in sampleDict["base"].keys(): #We next want to search matching base_loc in sampleDict
                    for sample_name in sampleDict["base"][base_loc]:
                        if sample_name in filtered_samples: #Only analyze samples of interest
                            if sampleDict["base"][base_loc][sample_name][1] >= min_reads:
                                #We want to calculate pairwise dissimilarity between sample and cellType and save into PDdict_temp
                                sample_methRate = float(sampleDict["base"][base_loc][sample_name][0] / sampleDict["base"][base_loc][sample_name][1])
                                # if not min(min_rate, 1 - min_rate) < sample_methRate < max(min_rate, 1 - min_rate): #We also want to only consider base locations in which sample is not between (min_rate, 1-min_rate)
                                dis = abs(sample_methRate - cellType_methRate) * 100
                                if sample_name not in PDdict_temp.keys():
                                    PDdict_temp[sample_name] = {}
                                if cellType_name not in PDdict_temp[sample_name].keys():
                                    PDdict_temp[sample_name][cellType_name] = {}
                                    PDdict_temp[sample_name][cellType_name]["PD"] = []
                                    PDdict_temp[sample_name][cellType_name]["region"] = set()
                                PDdict_temp[sample_name][cellType_name]["PD"].append(dis)
                                PDdict_temp[sample_name][cellType_name]["region"].update(region_set)
        f_cellType.close()
    #Save pairwise dissimilarity into PDdict for future processing/plotting
    cellType_list = []
    PDdict = {}
    PDdict["PD"] = {}
    PDdict["numCpG"] = {}
    PDdict["numRegion"] = {}
    PD_output = open(prefix + ".PD.txt", 'w')
    PD_output.write("Sample\tCell Type\tPairwise Dissimilarity\tNum CpG Shared\tNum Region Shared\n")
    for sample_name in sorted(PDdict_temp.keys()):
        PD_list = []
        num_CpG = [] #This contains the number of CpG sites compared
        num_region = [] #This contains the number of regions compared
        for cellType_name in sorted(PDdict_temp[sample_name].keys()):
            if len(PDdict_temp[sample_name][cellType_name]["PD"]) >= min_shared:
                pairwise_dist = float(sum(PDdict_temp[sample_name][cellType_name]["PD"])/len(PDdict_temp[sample_name][cellType_name]["PD"]))
                PD_output.write(sample_name + "\t" + cellType_name + "\t" + str(round(pairwise_dist, 4)) + "\t" + str(len(PDdict_temp[sample_name][cellType_name]["PD"])) + "\t" + str(len(PDdict_temp[sample_name][cellType_name]["region"])) + "\n")
                PD_list.append(pairwise_dist)
                num_CpG.append(len(PDdict_temp[sample_name][cellType_name]["PD"]))
                num_region.append(len(PDdict_temp[sample_name][cellType_name]["region"]))
                cellType_list.append(cellType_name)

            else:
                PD_list.append("NA")
                num_CpG.append(0)
                num_region.append(0)
        PDdict["PD"][sample_name] = PD_list
        PDdict["numCpG"][sample_name] = num_CpG
        PDdict["numRegion"][sample_name] = num_region
    PD_output.close()
    PDdict["index"] = sorted(list(set(cellType_list))) #This contains cell type names used for pandas dataframe index
    #Export PDdict to file using pickle
    with open(prefix + ".PD.pkl", 'wb') as PDdict_file:
        pickle.dump(PDdict, PDdict_file, protocol=pickle.HIGHEST_PROTOCOL)
    return PDdict

def calcPD_merge(sampleDict, typeDict, filtered_samples, DMR_bool, reg_bool, min_shared, min_rate, min_reads, min_CpG, prefix):
    '''
    We want to create a module that merges and calls average methylation calls across region of interest (DMR or reg or both)
    This average methylation is then used to compare single cells against reference with shared DMR or reg windows
    '''
    print("Importing cellType methylation rates per region (DMR or regulatory build)")
    cell_mergeDict = {} #Contains regions as keys and either base locations or cellType methylation rate as values for all cellType
    for cellType_name in tqdm(sorted(typeDict["cellType"].keys())):
        cell_mergeDict[cellType_name] = {}
        f_cellType = open(typeDict["cellType"][cellType_name]["file_loc"])
        for line in f_cellType:
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
            if reg_bool is True and 'reg' not in line:
                continue #Skip base_loc if not contained within regulatory build window and reg_bool is set to True
            region_list = base_info[7:]
            if int(total_reads) >= min_reads:
                for region in region_list:
                    if region not in cell_mergeDict[cellType_name].keys():
                        cell_mergeDict[cellType_name][region] = {}
                        cell_mergeDict[cellType_name][region]["base_loc"] = []
                        cell_mergeDict[cellType_name][region]["methRate"] = []
                    cell_mergeDict[cellType_name][region]["base_loc"].append(base_loc)
                    cell_mergeDict[cellType_name][region]["methRate"].append(cellType_methRate)
        f_cellType.close()

    print("Importing sample methylation rates per region (DMR or regulatory build) and intersecting with cellType")
    PDdict_temp = {} #Keep all pairwise dissimilarity comparisons in temp dictionary
    for cellType_name in tqdm(sorted(cell_mergeDict.keys())):
        for region in sorted(cell_mergeDict[cellType_name].keys()):
            sample_mergeDict = {} #Contains temporary methRate information for each sample_name for given region
            for base_loc in cell_mergeDict[cellType_name][region]["base_loc"]:
                if base_loc in sampleDict["base"].keys():
                    for sample_name in sampleDict["base"][base_loc]:
                        if sampleDict["base"][base_loc][sample_name][1] >= min_reads and sample_name in filtered_samples:
                            sample_methRate = float(int(sampleDict["base"][base_loc][sample_name][0]) / int(sampleDict["base"][base_loc][sample_name][1]))
                            #We want to save the methRate at given base for sample_name
                            if sample_name not in sample_mergeDict.keys():
                                sample_mergeDict[sample_name] = []
                            # sample_mergeDict[sample_name][region]["base_loc"].append(base_loc)
                            sample_mergeDict[sample_name].append(sample_methRate)
            for sample_name in sorted(sample_mergeDict.keys()):
                if sample_name not in PDdict_temp.keys():
                    PDdict_temp[sample_name] = {}
                if cellType_name not in PDdict_temp[sample_name].keys():
                    PDdict_temp[sample_name][cellType_name] = {}
                if len(sample_mergeDict[sample_name]) >= min_CpG and len(cell_mergeDict[cellType_name][region]["methRate"]) >= min_CpG:
                    PDdict_temp[sample_name][cellType_name][region] = {}
                    PDdict_temp[sample_name][cellType_name][region]["sample"] = np.asarray(sample_mergeDict[sample_name])
                    PDdict_temp[sample_name][cellType_name][region]["ref"] = np.asarray(cell_mergeDict[cellType_name][region]["methRate"])
    cell_mergeDict.clear() #Clear dictionary to save memory
    sample_mergeDict.clear()

    #Save pairwise dissimilarity into PDdict for future processing/plotting
    cellType_list = []
    PDdict = {}
    PDdict["PD"] = {}
    PDdict["numRegion"] = {}
    PDdict["numCpG"] = {}
    PD_output = open(prefix + ".PD.merge.txt", 'w')
    PD_output.write("Sample\tCell Type\tPairwise Dissimilarity\tNum Region Shared\tNum CpG Shared in Regions\n")
    print("Calculating PD for all samples")
    os.makedirs(prefix + ".stats")
    for sample_name in tqdm(sorted(PDdict_temp.keys())):
        PD_list = []
        numRegion_list = []
        numCpG_list = []
        stats_output = open(prefix + ".stats/" + prefix + "." + sample_name + ".PD.merge.stats.txt", 'w') #We want to save the regions that are
        for cellType_name in sorted(PDdict_temp[sample_name].keys()):
            numCpG = 0
            pairwise_dist = 0
            for region in sorted(PDdict_temp[sample_name][cellType_name].keys()):
                numCpG += len(PDdict_temp[sample_name][cellType_name][region]["sample"])
                dist = (np.mean(PDdict_temp[sample_name][cellType_name][region]["sample"]) - np.mean(PDdict_temp[sample_name][cellType_name][region]["ref"])) * 100
                if cellType_name == "HCT116":
                    if "DMR" in region:
                        (chrom, chromStart, chromEnd) = re.split(':|-', region)[1:]
                    elif "reg" in region:
                        (chrom, chromStart, chromEnd) = re.split(':|-', region)[2:]
                    stats_output.write("chr" + chrom.replace('chr','')+ "\t" + chromStart + "\t" + chromEnd + "\t" + str(round(dist, 2)) + "\t" + \
                        str(round(np.mean(PDdict_temp[sample_name][cellType_name][region]["sample"]), 2)) + "\t" + str(round(np.mean(PDdict_temp[sample_name][cellType_name][region]["ref"]), 2)) + "\t" + \
                        ','.join(str(round(i, 2)) for i in PDdict_temp[sample_name][cellType_name][region]["sample"]) + "\t" + \
                        ','.join(str(round(j, 2)) for j in PDdict_temp[sample_name][cellType_name][region]["ref"]) + "\n")
                pairwise_dist += np.absolute(dist)
            numRegion = len(PDdict_temp[sample_name][cellType_name].keys())
            if numRegion > 0:
                PD_list.append(pairwise_dist / numRegion)
            numRegion_list.append(numRegion)
            numCpG_list.append(numCpG)
            cellType_list.append(cellType_name)
            PD_output.write(sample_name + "\t" + cellType_name + "\t" + str(pairwise_dist) + "\t" + str(numRegion) + "\t" + str(numCpG) + "\n")
        PDdict["PD"][sample_name] = PD_list
        PDdict["numRegion"][sample_name] = numRegion_list
        PDdict["numCpG"][sample_name] = numCpG_list
        stats_output.close()
    PD_output.close()
    PDdict["index"] = sorted(list(set(cellType_list))) #This contains cell type names used for pandas dataframe index
    #Export PDdict to file using pickle
    with open(prefix + ".PD.merge.pkl", 'wb') as PDdict_file:
        pickle.dump(PDdict, PDdict_file, protocol=pickle.HIGHEST_PROTOCOL)
    return PDdict

def plotPD(PDdict, prefix, merge):
    #We want to ensure that the dataframe has the same number of columns (remove sample_name that don't have PD calculated for all reference cell types)
    print("Making PDdict dataframe")
    for sample_name in tqdm(sorted(PDdict["PD"].keys())):
        if len(PDdict["PD"][sample_name]) < len(PDdict["index"]):
            # print(sample_name + "\t" + ','.join(str(i) for i in PDdict["PD"][sample_name]))
            PDdict["PD"].pop(sample_name, None)
            PDdict["numCpG"].pop(sample_name, None)
            if merge is True:
                PDdict["numRegion"].pop(sample_name, None)
    PD_df = pandas.DataFrame(PDdict["PD"], index=sorted(PDdict["index"]))
    if merge is False:
        numCpG_df = pandas.DataFrame(PDdict["numCpG"], index=sorted(PDdict["index"]))
    else:
        numRegion_df = pandas.DataFrame(PDdict["numRegion"], index=sorted(PDdict["index"]))

    #We want to plot heatmaps of the PD comparison between cells and cell type
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

    #We want to plot heatmaps of the numCpG comparison between cells and cell type
    if merge is False:
        print("Drawing clustermap for numCpG")
        numCpG_clustermap = sns.clustermap(numCpG_df, xticklabels=True)
        numCpG_clustermap.savefig(prefix + ".numCpG.eps", format="eps", dpi=1000)

        print("Drawing clustermap for numCpG with Z-score across samples")
        numCpG_clustermap_z_sample = sns.clustermap(numCpG_df, z_score=1, xticklabels=True)
        numCpG_clustermap_z_sample.savefig(prefix + ".numCpG.z_sample.eps", format="eps", dpi=1000)

        print("Drawing clustermap for numCpG with Z-score across cellType")
        numCpG_clustermap_z_cellType = sns.clustermap(numCpG_df, z_score=0, xticklabels=True)
        numCpG_clustermap_z_cellType.savefig(prefix + ".numCpG.z_cellType.eps", format="eps", dpi=1000)
    else:
        print("Drawing clustermap for numRegion")
        numRegion_clustermap = sns.clustermap(numRegion_df, xticklabels=True)
        numRegion_clustermap.savefig(prefix + ".numRegion.eps", format="eps", dpi=1000)

        print("Drawing clustermap for numRegion with Z-score across samples")
        numRegion_clustermap_z_sample = sns.clustermap(numRegion_df, z_score=1, xticklabels=True)
        numRegion_clustermap_z_sample.savefig(prefix + ".numRegion.z_sample.eps", format="eps", dpi=1000)

        print("Drawing clustermap for numRegion with Z-score across cellType")
        numRegion_clustermap_z_cellType = sns.clustermap(numRegion_df, z_score=0, xticklabels=True)
        numRegion_clustermap_z_cellType.savefig(prefix + ".numRegion.z_cellType.eps", format="eps", dpi=1000)
    return

def refPD(sample_list, ref_info, sample_methDict, prefix, DMR, reg, merge, min_shared, min_rate, min_reads, min_CpG):
    '''
    This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039> and Paiwise Dissimilarity calculations from Hui 2018 Stem Cell Reports paper <https://doi.org/10.1016/j.stemcr.2018.07.003>.
    It will take as input methylation info for single cell samples and reference cellt ypes in order to calculate the pairwise dissimilarity between each
    '''
    if merge is False:
        f_pickle = prefix + ".PD.pkl"
    else:
        f_pickle = prefix + ".PD.merge.pkl"
    if os.path.exists(f_pickle) is False:
        #Import the samples that we want to keep for PD analysis
        filtered_samples = open(sample_list, 'r').read().splitlines()

        #Import pre-computed sample methDict
        print("Import sampleDict")
        with open(sample_methDict, 'rb') as sample_methDict_file:
            sampleDict = pickle.load(sample_methDict_file)

        #Import cell type file information into typeDict
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
        if merge is True:
            if True not in (DMR, reg):
                print("We cannot calculate merged pairwise dissimilarity without DMR or reg windows specified")
                return
            else:
                PDdict = calcPD_merge(sampleDict, typeDict, filtered_samples, DMR, reg, int(min_shared), float(min_rate), int(min_reads), int(min_CpG), prefix)
        else:
            PDdict = calcPD(sampleDict, typeDict, filtered_samples, DMR, reg, int(min_shared), float(min_rate), int(min_reads), prefix)
    else:
        print("Importing pre-computed PDdict")
        with open(f_pickle, 'rb') as PDdict_f:
            PDdict = pickle.load(PDdict_f)

    #Plot PD
    print("Plotting PD between single cells and cell type")
    plotPD(PDdict, prefix, merge)
