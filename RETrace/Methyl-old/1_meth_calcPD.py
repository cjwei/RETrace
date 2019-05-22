#!/usr/bin/env python3
import argparse
import pickle
from tqdm import tqdm

'''
Usage: python script.py --sample ... --cellType ... --prefix ... -stats --ref_CGI ...
This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039> and Paiwise Dissimilarity calculations from Hui 2018 Stem Cell Reports paper <https://doi.org/10.1016/j.stemcr.2018.07.003>
To run this script, we must include the following input:
    1) sample_methDict = pre-computed methDict containing single cell sample information.  The structure of methDict is as follows:
    methDict
        "sample"
            sample_name
                "file_loc" = file_loc
        "base"
            base_loc (chr:pos;chain)
                sample_name = (num_meth, total_reads)
    2) cellType_info = tab-delimited file containing file locations and cell type information for reference tsv files
    3) ref_CGI = reference bed file containing CGI locations across the desired reference genome
    4) prefix = prefix for output files
The script will then go through and output the following statistics/files:
    1) prefix.PD.txt = summarizes pairwise dissimilarity scores between samples and given reference cell-types
    2) prefix.covStats.txt = contains basic statistics of read coverage (i.e. num CpG covered, avg read coverage).  It will also include CpG island coverage statistics
    3) prefix.PD.pkl = exported dictionary containing pairwise-dissimilarity between samples and cell types
'''

def calcPD(sampleDict, typeDict, DMR_bool, seqDepth, min_shared, min_rate, prefix):
    #We want to iterate through all cell types and calculate pairwise dissimilarity when compared to all samples
    PDdict_temp = {} #Keep all pairwise dissimilarity comparisons in temp dictionary
    covDict = {} #We also want to keep track of all CpG base_loc and DMR's covered per sample
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
            if 'DMR' in line:
                DMR_loc = base_info[7].replace("DMR:", "")
            else:
                DMR_loc = None
                if DMR_bool is True:
                    continue #Skip base_loc if not contained within DMR and DMR_bool is set to True
            if not min(min_rate, 1 - min_rate) < cellType_methRate < max(min_rate, 1 - min_rate): #Only consider base locations in which reference methRate is not between (min_rate, 1-min_rate)
                if base_loc in sampleDict["base"].keys(): #We next want to search matching base_loc in sampleDict
                    for sample_name in sampleDict["base"][base_loc]:
                        #We want to save coverage data into covDict
                        if sample_name not in covDict.keys():
                            covDict[sample_name] = {}
                            covDict[sample_name]["all_CpG"] = set()
                            covDict[sample_name]["CpG"] = {}
                            covDict[sample_name]["DMR"] = {}
                        if cellType_name not in covDict[sample_name]["CpG"].keys():
                            covDict[sample_name]["CpG"][cellType_name] = set()
                        covDict[sample_name]["all_CpG"].add(base_loc) #Save all CpG's across all cellTypes captued per sample
                        covDict[sample_name]["CpG"][cellType_name].add(base_loc) #Save CpG's for each cellType individually
                        if DMR_loc is not None:
                            if cellType_name not in covDict[sample_name]["DMR"].keys():
                                covDict[sample_name]["DMR"][cellType_name] = set()
                            covDict[sample_name]["DMR"][cellType_name].add(DMR_loc) #Save DMR locations for each cellType individually
                        #Next we want to calculate pairwise dissimilarity between sample and cellType and save into PDdict_temp
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
    PD_output.write("Sample\tCell Type\tPairwise Dissimilarity\tNum Shared\n")
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
        if len(PD_list) > 0:
            PDdict["PD"][sample_name] = PD_list
    PD_output.close()

    #Export PDdict to file using pickle
    PDdict["index"] = sorted(list(set(cellType_list))) #This contains cell type names used for pandas dataframe index

    with open(prefix + ".PD.pkl", 'wb') as PDdict_file:
        pickle.dump(PDdict, PDdict_file, protocol=pickle.HIGHEST_PROTOCOL)

    return covDict

def printStats(covDict, DMR_bool, prefix):
    f_stats = open(prefix + ".PD.stats.txt", 'w')
    #Gather information for all sample_name and cellType_name
    sample_list = covDict.keys()
    cellType_list = []
    for sample_name in sample_list:
        cellType_list.extend(list(covDict[sample_name]["DMR"].keys()))
    cellType_list = list(set(cellType_list))
    #Make sure that there are no empty values in the covDict (fill in with an empty set if there is)
    for sample_name in sorted(sample_list):
        for cellType_name in sorted(cellType_list):
            if cellType_name not in covDict[sample_name]["CpG"].keys():
                covDict[sample_name]["CpG"][cellType_name] = set()
            if cellType_name not in covDict[sample_name]["DMR"].keys():
                covDict[sample_name]["DMR"][cellType_name] = set()
    #1) Print basic statitistics for number of CpG/DMR covered for each sample
    f_stats.write("---------- Raw Statistics (Total Num CpG, Num CpG (DMR) per cell type) ----------\n")
    if DMR_bool is True: #NOTE: If DMR is set, these are only sum of CpG counts found within DMR of all cell types
        f_stats.write("NOTE: Because DMR option is True, we are only considering CpG within DMRs across all cell types\n")
    f_stats.write("Sample,Num Total CpG," + ','.join(i + " CpG (DMR)" for i in sorted(cellType_list)) + "\n")
    for sample_name in sorted(sample_list):
        f_stats.write(sample_name + ',' + str(len(covDict[sample_name]["all_CpG"])) + ',') #We want to print out the total number of CpG across all cell types for sample
        for cellType_name in sorted(cellType_list):
            f_stats.write(str(len(covDict[sample_name]["CpG"][cellType_name])))
            f_stats.write('(' + str(len(covDict[sample_name]["DMR"][cellType_name])) + '),') #Print number of CpG and DMR covered per cell type per sample
        f_stats.write("\n")
    #2) Print number of shared CpG for each pairwise sample comparison
    f_stats.write("\n---------- Num Shared CpG between each sample----------\n")
    f_stats.write(',' + ','.join(sorted(sample_list)) + "\n")
    for sample_name1 in sorted(sample_list):
        f_stats.write(sample_name1 + ',')
        for sample_name2 in sorted(sample_list):
            f_stats.write(str(len(covDict[sample_name1]["all_CpG"].intersection(covDict[sample_name2]["all_CpG"]))) + ',')
        f_stats.write("\n")
    #3) Print number of DMR shared
    for cellType_name in sorted(cellType_list):
        f_stats.write("\n---------- Num DMR shared between each sample (cellType:" + cellType_name + ") ----------\n")
        f_stats.write(',' + ','.join(sorted(sample_list)) + "\n")
        for sample_name1 in sorted(sample_list):
            f_stats.write(sample_name1 + ',')
            for sample_name2 in sorted(sample_list):
                f_stats.write(str(len(covDict[sample_name1]["DMR"][cellType_name].intersection(covDict[sample_name2]["DMR"][cellType_name]))) + ',')
            f_stats.write("\n")
    f_stats.close()

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--sample', action="store", dest="sample_pkl", help="Pre-computed methDict containing sample information")
    parser.add_argument('--cellType_info', action="store", dest="cellType_info", help="Tab-delimited file containing the cell type tsv files/locations along with cell type names")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specifies prefix for output files")
    parser.add_argument('-DMR', action="store_true", default=False, help="[Optional] Constrain bases of interest in reference file to DMR regions")
    parser.add_argument('-covStats', action="store_true", default=True, help="Calculate CpG and DMR coverage statistics for each sample (pairwise between samples and cell types)")
    parser.add_argument('--seqDepth', action="store", dest="seqDepth", default=1, help="[Optional] Minimum number reads covering CpG for calculating pd (Default: 1)")
    parser.add_argument('--min_shared', action="store", dest="min_shared", default=100, help="[Optional] Minimum number CpG shared between each sample and cell type comparison (Default: 100)")
    parser.add_argument('--min_rate', action="store", dest="min_rate", default=0.5, help="[Optional] Minimum methylated or unmethylated rate for cellType reference for filtering confident methyl calls (Default: 0.5 = no filtering; 1 = only perfect call)")
    args = parser.parse_args()

    #Import cell type file information into cellDict
    print("Import typeDict")
    typeDict = {}
    typeDict["cellType"] = {}
    with open(args.cellType_info) as f_info:
        for cellType in f_info:
            (file_loc, cellType_name) = cellType.split()
            typeDict["cellType"][cellType_name] = {}
            typeDict["cellType"][cellType_name]["file_loc"] = file_loc

    #Import pre-computed sample methDict
    print("Import sampleDict")
    with open(args.sample_pkl, 'rb') as sample_file:
        sampleDict = pickle.load(sample_file)

    #Calculate pairwise dissimilarity matrix
    print("Calculating PD")
    covDict = calcPD(sampleDict, typeDict, args.DMR, int(args.seqDepth), int(args.min_shared), float(args.min_rate), args.prefix)

    #Print statistics if specified
    if args.covStats is True:
        printStats(covDict, args.DMR, args.prefix)

#%%
if __name__ == "__main__":
    main()
