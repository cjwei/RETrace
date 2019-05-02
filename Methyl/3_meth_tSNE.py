#!/usr/bin/env python3
import argparse
import os
import pickle
from tqdm import tqdm
import re

'''
Usage: python script.py --sample ... --gff ... --prefix ...
The script will calculate methRate across Ensembl Regulatory Build windows in order to perform NMF and tSNE plotting for single cells
To run this script, we must include the following input:
    1) sample_methDict = pre-computed methDict containing single cell sample information.  The structure of methDict is as follows:
    methDict
        "sample"
            sample_name
                "file_loc" = file_loc
        "base"
            base_loc (chr:pos;chain)
                sample_name = (num_meth, total_reads)
    2) gff = file containing location sof Ensembl Regulatory Build windows
    3) prefix = prefix for output files
The script will then go through and output the following files:
    1) prefix.methRate.csv = matrix containing methRate across Ensembl Regulatory Build
'''

def importReg(Ensembl_gff): #High memory requirement (~40gb for hg19 Ensembl Reg Build)
    regDict = {}
    #We want to import the Ensembl Regulatory Build windows from the gff file
    with open(Ensembl_gff) as f_gff:
        for line in f_gff:
            (chr, source, reg_type, start, end, score, strand, frame, attribute) = line.split("\t") #Based on gff format on Ensembl <https://uswest.ensembl.org/info/website/upload/gff.html>
            reg_name = reg_type + ':' + chr + ':' + start + '-' + end
            if chr not in regDict.keys():
                if chr.startswith('GL'): #We need to make a special condition for chromosomes that start with "GL" in order to match against sampleDict
                    chr = chr.replace('GL', 'gl').replace('.1', '')
                regDict[chr] = {}
            for pos in range(int(start), int(end) + 1): #We want to iterate through all base positions within reg window range
                regDict[chr][int(pos)] = reg_name
    return regDict

def regMeth(sampleDict, regDict, min_CpG, prefix):
    methDict = {}
    #We want to calculate the methRate across Ensembl Regulatory Build windows
    for base_loc in tqdm(sorted(sampleDict["base"].keys())):
        (chr, pos, chain) = re.split(':|;', base_loc)
        if 'gl' in chr: #We need tom ake a special condition for chromosomes that contain 'gl'
            chr = chr.split('_')[1]
        if chr in regDict.keys():
            if int(pos) in regDict[chr].keys():
                reg_name = regDict[chr][int(pos)]
                if reg_name not in methDict.keys():
                    methDict[reg_name] = {}
                for sample_name in sorted(sampleDict["base"][base_loc].keys()):
                    meth_reads = sampleDict["base"][base_loc][sample_name][0]
                    unmeth_reads = sampleDict["base"][base_loc][sample_name][1] - meth_reads
                    if sample_name not in methDict[reg_name].keys():
                        methDict[reg_name][sample_name] = {}
                        methDict[reg_name][sample_name]["meth_reads"] = 0
                        methDict[reg_name][sample_name]["unmeth_reads"] = 0
                        methDict[reg_name][sample_name]["CpG_sites"] = 0
                    methDict[reg_name][sample_name]["meth_reads"] += meth_reads
                    methDict[reg_name][sample_name]["unmeth_reads"] += unmeth_reads
                    methDict[reg_name][sample_name]["CpG_sites"] += 1
    #We then want calculate and output the methRate for each reg window and sample
    f_methRate = open(prefix + '.methRate.csv', 'w')
    f_methRate.write('reg_name,' + ','.join(sorted(sampleDict["sample"].keys())) + "\n")
    for reg_name in sorted(methDict.keys()):
        methRate_list = []
        for sample_name in sorted(sampleDict["sample"].keys()):
            if sample_name in methDict[reg_name].keys():
                meth_total = methDict[reg_name][sample_name]["meth_reads"]
                unmeth_total = methDict[reg_name][sample_name]["unmeth_reads"]
            else:
                meth_total = 0
                unmeth_total = 0
            methRate = float((meth_total + 1) / (meth_total + unmeth_total + 2))
            if methRate != 0.5 and methDict[reg_name][sample_name]["CpG_sites"] >= min_CpG: #We want to ensure the number of CpG sites captured per reg window is >= min_CpG
                methRate_list.append(methRate) #Because of the way we calculate, even if sample does not have reg_name, it'll be assigned a methRate of 0.5
            else:
                methRate_list.append('')
        f_methRate.write(reg_name + ',' + ','.join(str(x) for x in methRate_list) + "\n")
    f_methRate.close()
    return

def main():
    parser = argparse.ArgumentParser(description="Calculate methRate across Ensembl Regulatory Build windows")
    parser.add_argument('--sample', action="store", dest="sample_pkl", help="Pre-computed methDict containing sample information")
    parser.add_argument('--gff', action="store", dest="Ensembl_gff", help="gff file containing Ensembl Regulatory Build windows")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Output prefix for methRate.csv file")
    parser.add_argument('--min_CpG', action="store", dest="min_CpG", default=5, help="Minimum number of CpG sites within Ensembl Reg Build window")
    args = parser.parse_args()

    if os.path.exists(args.prefix + '.methRate.csv') is False:
        #Import pre-computed sample methDict
        print("Import sampleDict")
        with open(args.sample_pkl, 'rb') as sample_file:
            sampleDict = pickle.load(sample_file)

        #Import Ensemble Regulatory Build windows
        print("Import regDict")
        regDict = importReg(args.Ensembl_gff)

        #Calc methRate per sample_name in each reg window
        print("Calculating methRate across Ensembl Regulatory Build windows")
        regMeth(sampleDict, regDict, int(args.min_CpG), args.prefix)
    else:
        print(args.prefix + '.methRate.csv already exists')

if __name__ == "__main__":
    main()
