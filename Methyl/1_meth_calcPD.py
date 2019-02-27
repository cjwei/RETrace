#!/usr/bin/env python3
import argparse
import pickle

'''
Usage: python script.py --sample ... --cellType ... --prefix ... -stats --ref_CGI ...
This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039> and Paiwise Dissimilarity calculations from Hui 2018 Stem Cell Reports paper <https://doi.org/10.1016/j.stemcr.2018.07.003>
To run this script, we must include the following input:
    1) sample_methDict = pre-computed methDict containing single cell sample information.  The structure of methDict is as follows:
    methDict
            file_name
                "base"
                    base_loc
                        "Total" = total number reads
                        "Meth" = number of methylated reads
                        "Type" = C type (only save CpG [indicated as "CGN"], ignore CHH and CHG)
    2) cellType_methDict = pre-computed methDict containing reference cell type data
    3) ref_CGI = reference bed file containing CGI locations across the desired reference genome
    4) prefix = prefix for output files
The script will then go through and output the following statistics/files:
    1) prefix.PD.txt = summarizes pairwise dissimilarity scores between samples and given reference cell-types
    2) prefix.covStats.txt = contains basic statistics of read coverage (i.e. num CpG covered, avg read coverage).  It will also include CpG island coverage statistics
    3) prefix.PD.pkl = exported dictionary containing pairwise-dissimilarity between samples and cell types
'''

def calcCovStats(methDict, file_name, cutoff):
    (num_CpG, total_cov) = (0, 0)
    for base_loc in sorted(methDict[file_name]["base"].keys()):
        if (methDict[file_name]["base"][base_loc]["Total"] >= cutoff) and methDict[file_name]["base"][base_loc]["Type"] == "CpG":
            num_CpG += 1
            total_cov += methDict[file_name]["base"][base_loc]["Total"]
    if(num_CpG == 0):
        return (num_CpG, "NA")
    mean_cov = round(total_cov / num_CpG)
    return (num_CpG, mean_cov)

def calcStats(methDict, prefix):
    '''
    This funciton will calculate the following statistics per file:
        1) Total number of CpG positions covered
        2) Avg number of reads covering each position
        3) Number CGI's covered
    '''
    import pybedtools
    stats_output = open(prefix + ".covStats.txt", 'a+')
    stats_output.write("File\t"
        + "Unique CpG (1x)\tMean Coverage (1x)\t"
        + "Unique CpG (5x)\tMean Coverage (5x)\t"
        + "Unique CpG (10x)\tMean Coverage (10x)\t"
        + "Number CGI\tMean Coverage\n")
    for file_name in sorted(methDict.keys()):
        stats_output.write(file_name + "\t")
        #Summarize CGI stats
        num_CGI = len(methDict[file_name]["CGI"].keys())
        if (num_CGI > 0):
            mean_CGI_cov = round(sum(methDict[file_name]["CGI"].values())/num_CGI)
        else:
            mean_CGI_cov = "NA"
        #Determine CpG coverage
        for cutoff in [1,5,10]:
            (num_CpG, mean_cov) = calcCovStats(methDict, file_name, cutoff)
            stats_output.write(str(num_CpG) + "\t" + str(mean_cov) + "\t")
        stats_output.write(str(num_CGI) + "\t" + str(mean_CGI_cov) + "\n")
    stats_output.close()
    return

def calcPD(sampleDict, typeDict, typeDict_total, seqDepth, prefix):
    PDdict = {} #We want to save all paiwise_dis and p-values as dictionary where we have ordered lists for each file analyzed (ordered alphabetically by filename)
    PDdict["PD"] = {}
    PD_output = open(prefix + ".PD.txt", 'w')
    PD_output.write("Sample\tCell Type\tPairwise Dissimilarity\tNum Shared\n")
    for sample_name in sorted(sampleDict.keys()):
        PDdict["PD"][sample_name] = []
        for type_name in sorted(typeDict.keys()):
            dis_sum = 0
            num_shared = 0
            for shared_base in sampleDict[sample_name]["base"].keys() & typeDict[type_name]["base"].keys():
                if sampleDict[sample_name]["base"][shared_base]["Type"] == "CGN" and typeDict[type_name]["base"][shared_base]["Type"] == "CGN":
                    if sampleDict[sample_name]["base"][shared_base]["Total"] >= seqDepth and typeDict[type_name]["base"][shared_base]["Total"] >= seqDepth:
                        sample_methRate = float(sampleDict[sample_name]["base"][shared_base]["Meth"]/sampleDict[sample_name]["base"][shared_base]["Total"])
                        type_methRate = float(typeDict[type_name]["base"][shared_base]["Meth"]/typeDict[type_name]["base"][shared_base]["Total"])
                        if sample_methRate.is_integer() and type_methRate.is_integer():
                            num_shared += 1
                            if sample_methRate == type_methRate:
                                dis_sum += 0
                            else:
                                dis_sum += 100
            pairwise_dis = float(dis_sum/num_shared)
            PD_output.write(sample_name + "\t" + type_name + "\t" + str(round(pairwise_dis,4)) + "\t" + str(num_shared) + "\n")
            PDdict["PD"][sample_name].append(pairwise_dis)
    PD_output.close()

    #Export PDdict to file using pickle
    PDdict["index"] = sorted(typeDict.keys()) #This contains cell type names used for pandas dataframe index

    with open(prefix + ".PD.pkl", 'rb') as PDdict_file:
        pickle.dump(PDdict, PDdict_file, protocol=pickle.HIGHEST_PROTOCOL)

    return

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--sample', action="store", dest="sample_pkl", help="Pre-computed methDict containing sample information")
    parser.add_argument('--cellType', action="store", dest="cellType_pkl", help="Pre-computed methDict containing cell type information")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specifies prefix for output files")
    parser.add_argument('--seqDepth', action="store", dest="seqDepth", default=1, help="[Optional] Minimum number reads covering CpG for calculating pd")
    parser.add_argument('-stats', action="store_true", help="[Optional] Output statistics including unique CpG/CGI counts")
    parser.add_argument('--ref_CGI', action="store", dest="ref_CGI", help="Bed file containing reference genome CGI locations")
    args = parser.parse_args()

    #Import sample_methDict and cellType_methDict
    with open(args.sample_pkl, 'rb') as sample_file:
        sampleDict = pickle.load(sample_file)
    with open(args.cellType_pkl, 'rb') as cellType_file:
        typeDict = pickle.load(cellType_file)

    #Calculate methylation coverage statistics
    if args.stats is True:
        if args.ref_CGI:
            sampleDict = CGIstats(sampleDict, args.ref_CGI)
        else:
            parser.error('Cannot calculate statistics without specifying ref_CGI file')
        calcStats(sampleDict, args.prefix)

    #Calculate pairwise dissimilarity matrix
    calcPD(sampleDict, typeDict, int(args.seqDepth), args.prefix)

#%%
if __name__ == "__main__":
    main()
