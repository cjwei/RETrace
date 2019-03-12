#!/usr/bin/env python3
import argparse
from importer import parseProbes
from counter import counter
from count_visualizer import plotCounts
import os
import pysam
import numpy as np
from tqdm import tqdm
import pickle

'''
Usage: python script.py --input sample_list.txt --probes probes.info.txt --prefix output_prefix
This script will take as input a sample_list.txt file that gives a summary of all the following single cell data:
    1) File locaiton of sorted sample bam
    2) Sample name
    3) [Optional] Sample type (i.e. known clone in ex vivo tree)
The script will then go through and output a list of microsatellite subunit counts per target per sample.  It will then output the sampleDict that contains all msCounts in a pickle file
'''

#%%
def msCount():
    parser = argparse.ArgumentParser(description="Peform msCounts from mapped reads")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--probes', action="store", dest="probe_file", help="File location of probe info file summarizing targets captured")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specify output file containign pairwise distance calculations")
    parser.add_argument('--nproc', action="store", dest="nproc", default=10, help="Number of processes for msCounting")
    parser.add_argument('--counter', action="store", dest="count_type", default="aln", help="[simple/aln] Method for microsatellite subunit counting (default: aln)")
    parser.add_argument('--min_reads', action="store", dest="min_reads", default=10, help="Integer cutoff of number of reads required per target (default: 10)")
    parser.add_argument('-plot', action="store_true", help="Flag for indicating whether we want to output a plot file visualzing msCounts per targetID")
    args = parser.parse_args()

    #Parse sample_info file
    sampleDict = {}
    with open(args.sample_info) as f:
        for line in f:
            if len(line.split()) == 3: #If clone is specified
                (bam, sample, clone) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
                sampleDict[sample]["clone"] = clone
                sampleDict[sample]["msCount"] = {} #Contains list of all msCounts found within sample/target_id
            elif len(line.split()) == 2: #If clone is not specified
                (bam, sample) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
            else:
                print("Incorrect formatting for line (bam, sample, [optional] clone):\n" + "\t".join(line.split()))
                return

    #Perform msCount for each target_id within targetDict across all samples simultaneously (in order to decrease computational time)
    if not os.path.isfile(args.prefix + '.sampleDict.pkl'):
        print("Running msCount")
        #Import targetDict from probe_file
        targetDict = parseProbes(args.probe_file)

        for target_id in tqdm(sorted(targetDict.keys())):
            read_list = [] #Keep list of unique read sequences found across all samples
            for sample in sorted(sampleDict.keys()): #Iterate through all sample bams in order to extract all of the reads for given target
                samfile = pysam.AlignmentFile(sampleDict[sample]["bam"], "rb")
                for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
                    if read.query not in read_list:
                        read_list.append(read.query)
            msCount_list = counter(str(args.count_type), read_list, targetDict, target_id, int(args.nproc))
            print(target_id + "\t" + str(len(read_list)) + "\t" + ','.join(str(x) for x in msCount_list))
            #We want to re-iterate through samples in order to assign the msCount to the reads aligned for given sample
            for sample in sorted(sampleDict.keys()): #Iterate through all sample bams in order to assign msCounts
                sampleDict[sample]["msCount"][target_id] = []
                samfile = pysam.AlignmentFile(sampleDict[sample]["bam"], "rb")
                for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
                    try:
                        count_indx = read_list.index(read.query)
                    except:
                        print("Error: Unable to find read within sample")
                        return
                    if msCount_list[count_indx] > 0: #We only want to save valid msCounts for each sample
                        sampleDict[sample]["msCount"][target_id].append(msCount_list[count_indx])
                samfile.close()
                # print(target_id + "\t" + sample + "\t" + ','.join(str(x) for x in sampleDict[sample]["msCount"][target_id]))
                #We want to also save the information in targetDict for each sample (instead of sample as keys, we have target_id as keys)
                if len(sampleDict[sample]["msCount"][target_id]) > args.min_reads:
                    targetDict[target_id]["sample_msCount"][sample] = {}
                    targetDict[target_id]["sample_msCount"][sample]["count_list"] = sorted(set(sampleDict[sample]["msCount"][target_id]))
                    targetDict[target_id]["sample_msCount"][sample]["count_freq"] = [sampleDict[sample]["msCount"][target_id].count(i)/len(sampleDict[sample]["msCount"][target_id]) for i in sorted(set(sampleDict[sample]["msCount"][target_id]))]
        print("Pickling sampleDict.pkl and targetDict.pkl")
        pickle.dump(sampleDict, open(args.prefix + '.sampleDict.pkl', 'wb'))
        pickle.dump(targetDict, open(args.prefix + '.targetDict.pkl', 'wb'))
    else:
        print("Importing msCounts from " + args.prefix + ".sampleDict.pkl")
        sampleDict = pickle.load(open(args.prefix + '.sampleDict.pkl', 'rb'))
        targetDict = pickle.load(open(args.prefix + '.targetDict.pkl', 'rb'))

    #Visualize msCounts
    if args.plot is True:
        print("Plotting msCounts for each sample")
        plotCounts(args.prefix, targetDict)

#%%
if __name__ == "__main__":
    msCount()
