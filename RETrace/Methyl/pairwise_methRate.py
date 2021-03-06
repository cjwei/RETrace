#!/usr/bin/env python3
import pickle
from tqdm import tqdm
import re
import gzip

def importReg(Ensemble_gff): #High memory requirement (~40gb for hg19 Ensembl Reg Build)
    regDict = {}
    #We want to import the Ensembl Regulatory Build windows from the gff file
    with gzip.GzipFile(Ensemble_gff, 'rb') as f_gff:
        for line in f_gff:
            (chr, source, reg_type, start, end, score, strand, frame, attribute) = line.decode('utf-8').split("\t") #Based on gff format on Ensembl <https://uswest.ensembl.org/info/website/upload/gff.html>
            reg_name = reg_type + ':' + chr + ':' + start + '-' + end
            if chr not in regDict.keys():
                if chr.startswith('GL'): #We need to make a special condition for chromosomes that start with "GL" in order to match against sampleDict
                    chr = chr.replace('GL', 'gl').replace('.1', '')
                regDict[chr] = {}
            for pos in range(int(start), int(end) + 1): #We want to iterate through all base positions within reg window range
                regDict[chr][int(pos)] = reg_name
    return regDict

def calc_methRate(filtered_samples, sampleDict, regDict, prefix):
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
                    if sample_name in filtered_samples:
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
    f_methRate.write('reg_name,' + ','.join(sorted(filtered_samples)) + "\n")
    for reg_name in sorted(methDict.keys()):
        methRate_list = []
        for sample_name in sorted(filtered_samples):
            if sample_name in methDict[reg_name].keys():
                meth_total = methDict[reg_name][sample_name]["meth_reads"]
                unmeth_total = methDict[reg_name][sample_name]["unmeth_reads"]
            else:
                meth_total = 0
                unmeth_total = 0
            methRate = float((meth_total + 1) / (meth_total + unmeth_total + 2))
            # if methRate != 0.5 and methDict[reg_name][sample_name]["CpG_sites"] >= min_CpG: #We want to ensure the number of CpG sites captured per reg window is >= min_CpG
            #     methRate_list.append(methRate) #Because of the way we calculate, even if sample does not have reg_name, it'll be assigned a methRate of 0.5
            # else:
            #     methRate_list.append('')
            methRate_list.append(methRate) #Include if we do not want any empty calls (if no reads mapped, will have methRate of 0.5)
        f_methRate.write(reg_name + ',' + ','.join(str(x) for x in methRate_list) + "\n")
    f_methRate.close()

def pairwise_methRate(sample_list, sample_methDict, Ensemble_gff, prefix):
    '''
    The script will calculate methRate across Ensembl Regulatory Build windows in order to perform NMF and UMAP/tSNE plotting for single cells
    It will output a prefix.methRate.csv file summarizing the matrix containing methRate across Ensembl Regulatory Build windows for each pairwise comparison of cells
    '''

    #Import the samples that we want to keep for analysis
    filtered_samples = open(sample_list, 'r').read().splitlines()

    #Import pre-computed sample methDict
    print("Import sampleDict")
    with open(sample_methDict, 'rb') as sample_methDict_file:
        sampleDict = pickle.load(sample_methDict_file)

    #Import Ensemble Regulatory Build windows
    print("Import regDict (containing Ensemble Regulatory Build Windows)")
    regDict = importReg(Ensemble_gff)

    #Calculate pairwise methRate across Ensemble Regulatory Build windows
    print("Calculating pairwise methRate")
    calc_methRate(filtered_samples, sampleDict, regDict, prefix)
