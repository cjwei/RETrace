#!/usr/bin/env python3
import argparse
import pysam
import numpy as np
import random
import os
from tqdm import tqdm
import mmap

'''
Usage: python script.py --CA_bed CA_order.info.bed --bam input.bam --prefix ...
This scritp will take as input:
    1) CA_order.info.bed = bed file containing target locations of microsatellite capture
    2) input.bam = bam file containing reads mapped against hg19
The script will then determine the number of reads that mapped to microsatellite regions
'''

def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def num_targetReads():
    parser = argparse.ArgumentParser(description="Determine number of reads mapping to targeted microsatellite region")
    parser.add_argument("--CA_bed", action="store", dest="CA_bed", help="Bed file containing targetted microsatellite regions")
    # parser.add_argument("--ref_bed", action="store", dest="ref_bed", help="Bed file containing reference fragments")
    parser.add_argument("--bam", action="store", dest="bam_file", help="Bam file containing reads mapped to hg19")
    parser.add_argument("--prefix", action="store", dest="prefix", help="Prefix for depth data")
    args = parser.parse_args()

    depthDict = {}
    #We want to initialize depthDict for all readDepths <=100 reads (in order to make sure we cover all of the possible readDepths)
    for i in range(101):
        depthDict[i] = 0
    samfile = pysam.AlignmentFile(args.bam_file, "rb")
    with open(args.CA_bed) as f:
        for line in tqdm(f, total=get_num_lines(args.CA_bed)):
            num_complete = 0 #We want to keep track of reads that completely span the microsatellite sequence
            num_incomplete = 0 #We also want to keep track of the number of reads that fail to span the microsatellite
            (chrom, start, end) = line.split()[0:3]
            target_iter = samfile.fetch(chrom, int(start), int(end))
            for fastq_read in target_iter:
                (read_start, read_end) = (min(fastq_read.get_reference_positions()), max(fastq_read.get_reference_positions()))
                if read_start <= int(start) and read_end >= int(end):
                    num_complete += 1
                else:
                    num_incomplete += 1
            # print(chrom + ":" + start + "-" + end + "\t" + str(num_complete) + "\t" + str(num_incomplete))
            if num_complete not in depthDict.keys():
                if len([k for k in depthDict.keys() if k > num_complete]) > 0:
                    depthDict[num_complete] = depthDict[min(k for k in depthDict.keys() if k > num_complete)] #We want to use the number of targets for next greatest read_depth as current value
                else:
                    depthDict[num_complete] = 0
            for read_depth in sorted(depthDict.keys()): #Update all read_depth <= num_complete
                if read_depth <= num_complete:
                    depthDict[read_depth] += 1
    samfile.close()

    #We want to downsample the number of mapped reads in the sam file to calculate target coverage
    depth_out = open(args.prefix + ".depthDict.txt", 'w') #We want to print the number of targets captured per number of input reads based on the given
    depth_out.write("target_minReads\tnum_targets\n")
    for depth in sorted(depthDict.keys()):
        depth_out.write(args.prefix + "\t" + str(depth) + "\t" + str(depthDict[depth]) + "\n")
    depth_out.close()
    return

if __name__ == "__main__":
    num_targetReads()
