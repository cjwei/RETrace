#!/usr/bin/env python3
import argparse
import os
import subprocess
import pickle
from Bio import SeqIO #Remove when running on TSCC (does not have Bio module)
from Bio.Seq import Seq #Remove when running on TSCC (does not have Bio module)

'''
Usage: python script.py --cell_files ... --ref_fa ... --ref_CGI ... --methDict ...
This script will import methlylation information for all given cell types or single cells and export methDict.
The structure of methDict is as follows:
    methDict
        file_name
            "base"
                base_loc
                    "Total" = total number reads
                    "Meth" = number of methylated reads
                    "Type" = C type (only save CpG, ignore CHH and CHG)
'''

def parseMethCall(methDict):
    for sample_name in sorted(methDict["sample"].keys()):
        print("Importing:\t" + sample_name)
        f_methCall = open(methDict["sample"][sample_name]["file_loc"])
        for line in f_methCall:
            if line.startswith('#'):
                continue
            else:
                (chr, pos, chain, ctype, meth, total_bases) = line.split()[0:6]
                if ctype.startswith('CG'): #We only care about CGN methylated cytosine calls
                    base_loc = chr.replace("chr","") + ":" + pos + ";" + chain
                    if base_loc not in methDict["base"].keys():
                        methDict["base"][base_loc] = {}
                    methDict["base"][base_loc][sample_name] = (int(meth), int(total_bases))
    return methDict

def main():
    parser = argparse.ArgumentParser(description="Import methylation information for cell types or single cells and export methDict")
    parser.add_argument('--sample_info', action="store", dest="sample_info", help="Tab-delimited file containing the sample/cell type tsv files/locations along with sample/cell type names")
    parser.add_argument('--methDict', action="store", dest="methDict_file", help="File name containing exported methDict (*.pkl)")
    args = parser.parse_args()

    #Import methylation calls into methDict
    methDict = {}
    methDict["sample"] = {}
    methDict["base"] = {}
    with open(args.sample_info) as f_info:
        for sample in f_info:
            (file_loc, sample_name) = sample.split()
            methDict["sample"][sample_name] = {}
            methDict["sample"][sample_name]["file_loc"] = file_loc
    methDict = parseMethCall(methDict)

    with open(args.methDict_file, 'wb') as methDict_file:
        pickle.dump(methDict, methDict_file, protocol=pickle.HIGHEST_PROTOCOL)

#%%
if __name__ == "__main__":
    main()
