#!/usr/bin/env python3
import argparse
import os
import subprocess
import pickle
from Bio import SeqIO #Remove when running on TSCC (does not have Bio module)
from Bio.Seq import Seq #Remove when running on TSCC (does not have Bio module)
import pybedtools

'''
Usage: python script.py --sample_info ... --methDict ...
This script will import all CpG methylation signal from sample tsv files into a methDict dictionary to speed up processing of pairwise dissimilarity calculation.
The structure of the outputted methDict is as follows:
methDict
    "sample"
        sample_name
            file_loc
    "base"
        base_loc (chr:pos;chain)
            sample_name
                (num_meth, total_reads)
'''

def CGIstats(methDict, ref_CGI):
    for sample_name in sorted(methDict["sample"].keys()):
        methDict["sample"][sample_name]["CGI"] = {}
        sample_bed = pybedtools.BedTool(methDict["sample"][sample_name]["bam_file"])
        ref_bed = pybedtools.BedTool(ref_CGI)
        CGI_intersect = ref_bed.intersect(sample_bed, bed=True, wa=True, wb=True)
        for CGI in CGI_intersect:
            CGI_loc = CGI.fields[0] + ":" + CGI.fields[1] + "-" + CGI.fields[2]
            read_name = CGI.fields[7]
            if CGI_loc not in methDict["sample"][sample_name]["CGI"].keys():
                methDict["sample"][sample_name]["CGI"][CGI_loc] = 1
            else:
                methDict["sample"][sample_name]["CGI"][CGI_loc] += 1
    return methDict

def calcStats(methDict, prefix):
    '''
    This funciton will calculate the following statistics per file:
        1) Total number of CpG positions covered
        2) Avg number of reads covering each position
        3) Number CGI's covered
    '''
    stats_output = open(prefix + ".covStats.txt", 'w')
    stats_output.write("Sample\t"
        + "Unique CpG (1x)\tMean Coverage (1x)\t"
        + "Unique CpG (5x)\tMean Coverage (5x)\t"
        + "Unique CpG (10x)\tMean Coverage (10x)\t"
        + "Number CGI\tMean Coverage\n")
    #Parse saved methDict data to calculate read coverage stats
    for base_loc in sorted(methDict["base"].keys()):
        for sample_name in sorted(methDict["base"][base_loc].keys()):
            if not "Cov" in methDict["sample"][sample_name].keys():
                methDict["sample"][sample_name]["Cov"] = [] #This contains a list of all CpG coverage values for each CpG site
            methDict["sample"][sample_name]["Cov"].append(methDict["base"][base_loc][sample_name][1])
    for sample_name in sorted(methDict["sample"].keys()):
        stats_output.write(sample_name + "\t")
        #Summarize CGI stats
        num_CGI = len(methDict["sample"][sample_name]["CGI"].keys())
        if (num_CGI > 0):
            mean_CGI_cov = round(sum(methDict["sample"][sample_name]["CGI"].values())/num_CGI, 2)
        else:
            mean_CGI_cov = "NA"
        #Determine CpG coverage for various cutoff values
        for cutoff in [1,5,10]:
            filtered_CpG = [CpG for CpG in methDict["sample"][sample_name]["Cov"] if CpG >= cutoff]
            num_CpG = len(filtered_CpG)
            if num_CpG > 0:
                mean_CpG_cov = round(sum(filtered_CpG)/len(filtered_CpG), 2)
            else:
                mean_CpG_cov = "NA"
            stats_output.write(str(num_CpG) + "\t" + str(mean_CpG_cov) + "\t")
        stats_output.write(str(num_CGI) + "\t" + str(mean_CGI_cov) + "\n")
    return

def parseMethCall(methDict):
    for sample_name in sorted(methDict["sample"].keys()):
        print("Importing:\t" + sample_name)
        f_methCall = open(methDict["sample"][sample_name]["tsv_file"])
        for line in f_methCall:
            if line.startswith('#'):
                continue
            else:
                (chr, pos, chain, ctype, meth, total_reads) = line.split()[0:6]
                if ctype.startswith('CG'): #We only care about CGN methylated cytosine calls
                    base_loc = chr.replace("chr","") + ":" + pos + ";" + chain
                    if base_loc not in methDict["base"].keys():
                        methDict["base"][base_loc] = {}
                    methDict["base"][base_loc][sample_name] = (int(meth), int(total_reads))
    return methDict

def main():
    parser = argparse.ArgumentParser(description="Import methylation information for cell types or single cells and export methDict")
    parser.add_argument('--sample_info', action="store", dest="sample_info", help="Tab-delimited file containing the sample/cell type tsv files/locations along with sample/cell type names")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Prefix name for exporting methDict and coverage statistics information")
    parser.add_argument('-stats', action="store_true", help="[Optional] Output statistics including unique CpG and CGI counts")
    parser.add_argument('--ref_CGI', action="store", dest="ref_CGI", help="Bed file containing reference genome CGI locations")
    args = parser.parse_args()

    #Import methylation calls into methDict
    methDict = {}
    methDict["sample"] = {}
    methDict["base"] = {}
    with open(args.sample_info) as f_info:
        for sample in f_info:
            (bam_file, tsv_file, sample_name) = sample.split()
            methDict["sample"][sample_name] = {}
            methDict["sample"][sample_name]["bam_file"] = bam_file
            methDict["sample"][sample_name]["tsv_file"] = tsv_file
    methDict = parseMethCall(methDict)

    #Export methDict prior to calculating stats file
    with open(args.prefix + ".methDict.pkl", 'wb') as methDict_file:
        pickle.dump(methDict, methDict_file, protocol=pickle.HIGHEST_PROTOCOL)

    #We also want to calculate methylation coverage stats
    if args.stats is True:
        if args.ref_CGI:
            methDict = CGIstats(methDict, args.ref_CGI)
        else:
            parser.error('Cannot calculate statistics without specifying ref_CGI file')
        calcStats(methDict, args.prefix)

#%%
if __name__ == "__main__":
    main()
