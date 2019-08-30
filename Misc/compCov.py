#!/usr/bin/env python3
import argparse
import pickle
from tqdm import tqdm
import re
import gzip

def importBed(bed_file, gff_flag):
    bedDict = {}
    #We want to impor the bed_file into a dictonary
    with gzip.open(bed_file) as f_bed:
        for line in f_bed:
            line = line.decode("utf-8")
            if gff_flag is True:
                chrom = line.split()[0]
                start = line.split()[3]
                end = line.split()[4]
            else:
                (chrom, start, end) = line.split()[0:3]
            region_name = chrom + ":" + start + "-" + end
            if chrom.replace('chr','') not in bedDict.keys():
                bedDict[chrom.replace('chr','')] = {}
            for pos in range(int(start), int(end) + 1): #We want to iterate through all base positions
                bedDict[chrom.replace('chr','')][int(pos)] = (region_name)
    return bedDict

def compCov(infoDict, sampleDict, bedDict, interest_group, prefix):
    covDict = {}
    covDict["interest_group"] = set()
    covDict["background"] = set()
    for base_loc in tqdm(sorted(sampleDict["base"].keys())):
        #Only want to consider base_loc that are found in bedDict
        (chrom, pos, chain) = re.split(':|;', base_loc)
        if chrom.replace('chrom','') in bedDict.keys():
            if int(pos) in bedDict[chrom].keys():
                for sample in sampleDict["base"][base_loc].keys():
                    if sample in infoDict["group"][interest_group]:
                        covDict["interest_group"].add(bedDict[chrom][int(pos)])
                        covDict["background"].add(bedDict[chrom][int(pos)])
                        if sample not in covDict.keys():
                            covDict[sample] = set() #We want to keep track of all individual regions covered for each individual sample within interest_group
                        covDict[sample].add(bedDict[chrom][int(pos)])
                    elif sample in infoDict["group"]["background"]:
                        covDict["background"].add(bedDict[chrom][int(pos)])
    for group_name in covDict.keys():
        f_output = open(prefix + "." + group_name + ".bed", 'w')
        for region_name in covDict[group_name]:
            (chrom, start, end) = re.split(':|-', region_name)
            f_output.write("chr" + chrom.replace('chr','') + "\t" + start + "\t" + end + "\n")
        f_output.close()

def main():
    '''
    This script will take as input a sample_info file that gives details on the samples we want to compare:
        sampleName1    Outgroup
        sampleName2
        sampleName3
        ...
    The script will then merge all of the single cells labeled as the group we want to analyze.
    It'll then merge the single cells into a pseudo-bulk containing all the CpGs (or regions if a bed file is provided) that were covered.
    The script will then calculate the CGI (or any other region that is provided in the bed file) that are covered between the group of interest only, background only, or all cells.  Output format is 3 bed files:
        chrom   chromStart  chromEnd    region
    '''
    parser = argparse.ArgumentParser(description="Compare methyl signal from single cells against pseudo-bulk")
    parser.add_argument("--sample_info", action="store", dest="sample_info", help="Specify sample_info that contains the single cells and bulk to analyze")
    parser.add_argument("--methDict", action="store", dest="sample_methDict", help="methDict containing all meth calls for analysis")
    parser.add_argument("--bed", action="store", dest="bed_file", help="Bed file containing regions of interest")
    parser.add_argument("--group", action="store", dest="interest_group", help="Group label of single cells of interest")
    parser.add_argument("--prefix", action="store", dest="prefix", help="Output prefix")
    parser.add_argument("-gff", action="store_true", default=False, help="Flag for specifying gff (not bed) format especially used for Ensembl regulatory windows")
    args = parser.parse_args()

    infoDict = {}
    infoDict["group"] = {}
    infoDict["group"]["background"] = [] # We want to set background as all single cells that are unlabeled
    with open(args.sample_info) as f_info:
        for line in f_info:
            if len(line.split()) == 2:
                (sample_name, sample_group) = line.split()
                if sample_group not in infoDict["group"].keys():
                    infoDict["group"][sample_group] = []
                infoDict["group"][sample_group].append(sample_name)
            else:
                sample_name = line.rstrip()
                infoDict["group"]["background"].append(sample_name)

    #Import pre-computed sample methDict
    print("Import sampleDict")
    with open(args.sample_methDict, 'rb') as sample_methDict_file:
        sampleDict = pickle.load(sample_methDict_file)

    #Import bed file containing regions of interest
    print("Import bed containing region info")
    bedDict = importBed(args.bed_file, args.gff)

    #Import methyl calls and asign to appropriate sample_group/sample_name
    print("Comparing coverage across regions")
    compCov(infoDict, sampleDict, bedDict, args.interest_group, args.prefix)

if __name__ == "__main__":
    main()
