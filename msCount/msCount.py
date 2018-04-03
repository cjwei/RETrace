#!/usr/bin/env python
from __future__ import print_function
import argparse
from importer import parseProbes

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
    parser.add_argument('--input', action="store", dest="input_bam")
    args = parser.parse_args()
    
    targetDict = parseProbes(args.probe_file)
    
#%%
if __name__ == "__main__":
    msCount()