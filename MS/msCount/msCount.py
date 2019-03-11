#!/usr/bin/env python3
import argparse
from importer import parseProbes
from counter import counter
from count_visualizer import plotCounts
import pysam
import numpy as np
from tqdm import tqdm

'''
Usage: python script.py --probes probes.info.txt --bam input.sorted.bam
This script will take as input:
    1) probes.info.txt = contains all the probes designed per target_id
    2) input.sorted.bam = sorted bam file containing mapping results against hg19
The script will then go through and output a list of microsatellite subunit counts per target
'''

#%%
def msCount():
    parser = argparse.ArgumentParser(description="Perform msCounts from mapped reads")
    parser.add_argument('--probes', action="store", dest="probe_file")
    parser.add_argument('--bam', action="store", dest="input_bam")
    parser.add_argument('--num_reads', action="store", dest="num_reads", default=10, help="Integer cutoff of number of reads required per target (default: 10)")
    parser.add_argument('--counter', action="store", dest="count_type", default="aln", help="[simple/aln] Method for microsatellite subunit counting (default: aln)")
    parser.add_argument('--plot_file', action="store", dest="plot_file", default=None, help="Specify output plot_file if you want to plot msCounts for each target")
    parser.add_argument('--nproc', action="store", dest="nproc", default=5, help="Number of proc for msCounting")
    args = parser.parse_args()

    #Import targetDict from probe_file
    targetDict = parseProbes(args.probe_file)

    #Perform msCount for each target_id within targetDict
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    for target_id in tqdm(sorted(targetDict.keys())):
        read_list = []
        for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
            read_list.append(read.query)
        if len(read_list) >= args.num_reads:
            targetDict = counter(str(args.count_type), read_list, targetDict, target_id, int(args.nproc))
            print(target_id + "\t" + str(len(read_list)) + "\t " + ",".join(str(int(count)) for count in targetDict[target_id]["count_list"]))
    samfile.close()

    #Visualize msCounts
    if args.plot_file is not None:
        plotCounts(args.plot_file, targetDict)

#%%
if __name__ == "__main__":
    msCount()
