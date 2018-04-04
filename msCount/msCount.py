#!/usr/bin/env python
import argparse
from importer import parseProbes
import pysam

'''
Usage: python script.py --probes probes.info.txt --input input.sorted.bam
This script will take as input:
    1) probes.info.txt = contains all the probes designed per targetID
    2) input.sorted.bam = sorted bam file containing mapping results against hg19
The script will then go through and output a list of microsatellite subunit counts per target
'''

#%%
def msCount():
    parser = argparse.ArgumentParser(description="Perform msCounts from mapped reads")
    parser.add_argument('--probes', action="store", dest="probe_file")
    parser.add_argument('--bam', action="store", dest="input_bam")
    parser.add_argument('--num_reads', action="store", dest="num_reads", default=10)
    args = parser.parse_args()
    
    #Import targetDict from probe_file
    targetDict = parseProbes(args.probe_file)
    
    #Perform msCount for each targetID within targetDict
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    for targetID in sorted(targetDict.keys()):
        for read in samfile.fetch(targetDict[targetID]["chrom"], int(targetDict[targetID]["chromStart"]), int(targetDict[targetID]["chromEnd"])):
            print(read)
    samfile.close()
    
#%%
if __name__ == "__main__":
    msCount()