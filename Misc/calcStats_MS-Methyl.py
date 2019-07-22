#!/usr/bin/env python3
import pysam
import argparse
from tqdm import tqdm
import gzip
import itertools
import random
import os
import pickle

def parseMS(sampleDict, f_target):
    '''
    This will go through the bam file and count number of reads that completely cover the MS locus
    Simpler analysis to as the analysis found in <calcDepth.py>
    '''
    print("Calculating MS Statistics")
    for sample_name in sorted(sampleDict.keys()):
        print("\tImporting:\t" + sample_name)
        samfile = pysam.AlignmentFile(sampleDict[sample_name]["MS"]["bam"])
        with open(f_target) as f:
            for line in f:
                num_complete = 0
                (chrom, chromStart, chromEnd) = line.split()[0:3]
                target_iter = samfile.fetch(chrom, int(chromStart), int(chromEnd))
                for fastq_read in target_iter:
                    (read_start, read_end) = (min(fastq_read.get_reference_positions()), max(fastq_read.get_reference_positions()))
                    if read_start <= int(chromStart) and read_end >= int(chromEnd):
                        num_complete += 1
                target_id = chrom + ":" + chromStart + "-" + chromEnd
                sampleDict[sample_name]["MS"]["depthDict"][target_id] = num_complete
        samfile.close()
    return sampleDict

def parseMethyl(sampleDict):
    print("Calculating Methyl Statistics")
    for sample_name in sorted(sampleDict.keys()):
        print("\tImporting:\t" + sample_name)
        if "gz" in sampleDict[sample_name]["Methyl"]["tsv"]:
            f_methCall = gzip.open(sampleDict[sample_name]["Methyl"]["tsv"], 'rb')
        else:
            f_methCall = open(sampleDict[sample_name]["Methyl"]["tsv"])
        for line in f_methCall:
            if "gz"in sampleDict[sample_name]["Methyl"]["tsv"]:
                line = line.decode('utf-8')
            if line.startswith('#'):
                continue
            else:
                (chr, pos, chain, ctype, meth, total_reads) = line.split()[0:6]
                base_loc = chr.replace("chr","") + ":" + pos + ";" + chain
                if ctype.startswith('CG'): #We only care about CGN methylated cytosine calls
                    sampleDict[sample_name]["Methyl"]["depthDict"][base_loc] = int(total_reads)
    return sampleDict

def printStats(sampleDict, prefix, min_MS, min_Methyl):
    #We want to first print statistics for each single cell (or bulk sample) separately
    print("Analyzing single cell stats")
    f_out_single = open(prefix + '.singleStats.txt', 'w')
    f_out_single.write("sample_name\tnum_MS\tnum_Methyl\n")
    for sample_name in sorted(sampleDict.keys()):
        num_MS = 0
        for target_id in sampleDict[sample_name]["MS"]["depthDict"].keys():
            if sampleDict[sample_name]["MS"]["depthDict"][target_id] >= min_MS:
                num_MS += 1
        num_Methyl = 0
        for base_loc in sampleDict[sample_name]["Methyl"]["depthDict"].keys():
            if sampleDict[sample_name]["Methyl"]["depthDict"][base_loc] >= min_Methyl:
                num_Methyl += 1
        f_out_single.write(sample_name + "\t" + str(num_MS) + "\t" + str(num_Methyl) + "\n")
    f_out_single.close()

    #We also want to determine the number of target_id and base_loc that are present from merging 2-4 single cells together
    f_out_merge = open(prefix + '.mergeStats.txt', 'w')
    f_out_merge.write("sample_merge\tnum_merge\tnum_MS\tnum_Methyl\n")
    for n_merge in [2, 10]:
        SC_samples = [sample_name for sample_name in sampleDict.keys() if "SC" in sample_name]

        # #We want to perform reservoir downsampling of all possible combinations of single cell samples to produce a random sampling of 1000 merged samples
        # merge_list = []
        # for indx, sample_merge in tqdm(enumerate(itertools.combinations(SC_samples, n_merge))):
        #     merge_name = ';'.join(sample_merge)
        #     if indx < 1000:
        #         merge_list.append(merge_name)
        #     elif indx >= 1000 and random.random() < 1000/float(indx + 1):
        #         replace_indx = random.randint(0, len(merge_list) - 1)
        #         merge_list[replace_indx] = merge_name

        #We want to randomly build merge_list from the ground up
        print("Analyzing merged single cells:\t" + str(n_merge))
        merge_set = set()
        while len(merge_set) < 1000:
            merge_name = ';'.join(sorted(random.sample(SC_samples, n_merge)))
            merge_set.add(merge_name)

        for merge_name in sorted(merge_set):
            sample_merge = merge_name.split(';')
            MS_analyzed = set()
            num_MS = 0
            for sample_name in sample_merge:
                for target_id in sampleDict[sample_name]["MS"]["depthDict"].keys():
                    if target_id not in MS_analyzed:
                        if sampleDict[sample_name]["MS"]["depthDict"][target_id] >= min_MS:
                            MS_analyzed.add(target_id)
                            num_MS += 1
            Methyl_analyzed = set()
            num_Methyl = 0
            for sample_name in sample_merge:
                for base_loc in sampleDict[sample_name]["Methyl"]["depthDict"].keys():
                    if base_loc not in Methyl_analyzed:
                        if sampleDict[sample_name]["Methyl"]["depthDict"][base_loc] >= min_Methyl:
                            Methyl_analyzed.add(base_loc)
                            num_Methyl += 1
            f_out_merge.write(','.join(sample_merge) + "\t" + str(n_merge) + "\t" + str(num_MS) + "\t" + str(num_Methyl) + "\n")
    return

def calcStats():
    '''
    This script will calculate the following coverage statistics for MS/Methyl data:
        sample  MS_loci CpG_sites
    We will take as input a sample_info file that contains given samples (SC or bulk) that we want to calculate coverage stats:
        sample_name MS_bam  Methyl_tsv
    '''
    parser = argparse.ArgumentParser(description="Calculate number of MS loci and CpG sites captured for given samples")
    parser.add_argument("--sample_info", action="store", dest="sample_info", help="Specify sample_info file containing MS/methyl file loc")
    parser.add_argument("--target_bed", action="store", dest="f_target", help="Bed file containing target loci locations")
    parser.add_argument("--prefix", action="store", dest="prefix", help="Output file containing information for plotting")
    parser.add_argument("--min_MS", action="store", dest="min_MS", default=30, help="Minimum number of reads for MS loci")
    parser.add_argument("--min_Methyl", action="store", dest="min_Methyl", default=1, help="Minimum number of reads for Methyl CpG sites")
    args = parser.parse_args()

    #Import sample_info
    if not os.path.isfile(args.prefix + ".sampleDict.pkl"):
        sampleDict = {}
        with open(args.sample_info, 'r') as f_info:
            for line in f_info:
                (sample_name, MS_bam, Methyl_tsv) = line.rstrip().split()
                sampleDict[sample_name] = {}
                sampleDict[sample_name]["MS"] = {}
                sampleDict[sample_name]["MS"]["bam"] = MS_bam
                sampleDict[sample_name]["MS"]["depthDict"] = {} #Contains target_id as keys and num reads as values
                sampleDict[sample_name]["Methyl"] = {}
                sampleDict[sample_name]["Methyl"]["tsv"] = Methyl_tsv
                sampleDict[sample_name]["Methyl"]["depthDict"] = {} #Contains CpG base_loc as keys and num reads as values
        #Import MS statistics
        sampleDict = parseMS(sampleDict, args.f_target)
        #Import Methyl statistics
        sampleDict = parseMethyl(sampleDict)
        pickle.dump(sampleDict, open(args.prefix + ".sampleDict.pkl", "wb"))
    else:
        sampleDict = pickle.load(open(args.prefix + ".sampleDict.pkl", "rb"))

    #Calculate number of MS loci and Methyl CpG base_loc per single cell or merge of cells
    printStats(sampleDict, args.prefix, args.min_MS, args.min_Methyl)

if __name__ == "__main__":
    calcStats()
