#!/usr/bin/env python3
import argparse
from importer import import_targetDict, import_sampleDict
from counter import counter
from count_visualizer import plotCounts
import os
import pysam
import numpy as np
from tqdm import tqdm
import pickle
from numba import jit
import multiprocessing
manager = multiprocessing.Manager()

'''
Usage: python script.py --input sample_list.txt --probes probes.info.txt --prefix output_prefix
This script will take as input a sample_list.txt file that gives a summary of all the following single cell data:
    1) File locaiton of sorted sample bam
    2) Sample name
    3) Sex (mostly used as information for allelotyping)
    4) [Optional] Sample type (i.e. known clone in ex vivo tree)
The script will then go through and output a list of microsatellite subunit counts per target per sample.  It will then output the targetDict that contains all msCounts in a pickle file with the following structure:
    targetDict
        target_id
            "sample"
                sample
                    msCount_list = list containing all msCounts for given sample at target_id
'''

#%%
def multi_counter(count_type, target_group, targetDict, sampleDict, counterDict, prefix): #This module allow for multiprocessing of the counter by allowing target_id's to be processed in parallel
    for target_id in sorted(target_group):
        read_list = [] #Keep list of unique read sequences found across all samples
        for sample in sorted(sampleDict.keys()): #Iterate through all sample bams in order to extract all of the reads for given target
            samfile = pysam.AlignmentFile(sampleDict[sample]["bam"], "rb")
            for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
                if read.query not in read_list:
                    read_list.append(read.query)
        counterDict[target_id] = counter(count_type, read_list, targetDict, target_id, prefix)
    return

def msCount():
    parser = argparse.ArgumentParser(description="Peform msCounts from mapped reads")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--probes', action="store", dest="probe_file", help="File location of probe info file summarizing targets captured")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specify output prefix (for targetDict and alleleDict, along with any stats or plot files)")
    parser.add_argument('--nproc', action="store", dest="nproc", type=int, default=10, help="Number of processes for msCounting")
    parser.add_argument('--counter', action="store", dest="count_type", default="aln", help="[simple/aln] Method for microsatellite subunit counting (default: aln)")
    parser.add_argument('--min_reads', action="store", dest="min_reads", type=int, default=10, help="Integer cutoff of number of reads required per target (default: 10)")
    parser.add_argument('-plot', action="store_true", help="Flag for indicating whether we want to output a plot file visualzing msCounts per targetID")
    args = parser.parse_args()

    #Parse sample_info file
    sampleDict = import_sampleDict(args.sample_info)
    if not sampleDict:
        return

    #Perform msCount for each target_id within targetDict across all samples simultaneously (in order to decrease computational time)
    if not os.path.isfile(args.prefix + '.targetDict.pkl'):
        print("Running msCount")
        #Import targetDict from probe_file
        targetDict = import_targetDict(args.probe_file)

        #We want to perform multiprocessing using the args.nproc specified.  In order to do this, we can split target_id's into chunks of sized args.nproc
        #This is based off of: <https://further-reading.net/2017/01/quick-tutorial-python-multiprocessing/>, <https://stackoverflow.com/questions/10415028/how-can-i-recover-the-return-value-of-a-function-passed-to-multiprocessing-proce>
        counterDict = manager.dict() #This allows for parallel processing of msCount for multiple target_id at once
        jobs = []
        for target_group in [sorted(targetDict.keys())[i::args.nproc] for i in range(args.nproc)]:
            p = multiprocessing.Process(target = multi_counter, args = (args.count_type, target_group, targetDict, sampleDict, counterDict, args.prefix))
            jobs.append(p)
            p.start()

        #Join counterDict
        for p in jobs:
            p.join()

        #Re-iterate through samples in order to assign msCount to the reads aligned for given sample (save in targetDict)
        for target_id in sorted(counterDict.keys()):
            #We want to extract all read and msCount info from counterDict
            read_list = [str(count_info.split('=')[0]) for count_info in counterDict[target_id]]
            msCount_list = [int(count_info.split('=')[1]) for count_info in counterDict[target_id]]
            #Assign appropriate msCounts to each sample and save in targetDict
            targetDict[target_id]["sample"] = {}
            for sample in sorted(sampleDict.keys()):
                targetDict[target_id]["sample"][sample] = []
                samfile = pysam.AlignmentFile(sampleDict[sample]["bam"], "rb")
                for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
                    try:
                        count_indx = read_list.index(read.query)
                    except:
                        print("Error: Unable to find read within sample")
                        return
                    if msCount_list[count_indx] > 0: #We only want to save valid msCounts for each sample
                        targetDict[target_id]["sample"][sample].append(msCount_list[count_indx])
                samfile.close()
                if len(targetDict[target_id]["sample"][sample]) < args.min_reads: #We want to remove sample from targetDict if there are fewer than min_reads
                    targetDict[target_id]["sample"].pop(sample, None)
        print("Pickling targetDict.pkl")
        pickle.dump(targetDict, open(args.prefix + '.targetDict.pkl', 'wb'))

    # Visualize msCounts
    if args.plot is True:
        print("Plotting msCounts for each sample")
        plotCounts(args.prefix, targetDict)

#%%
if __name__ == "__main__":
    msCount()
