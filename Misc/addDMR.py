#!/usr/bin/env python3
import argparse
from tqdm import tqdm

'''
Usage: python script.py --tsv ... --DMR ... --output
This script will take as input the reference tsv file created from methylpy:
    8       11048   -       CGN     3       3       1
    8       11065   -       CGA     3       3       1
    8       11102   -       CGG     6       6       1
    8       11788   +       CGT     1       1       1
    8       11789   -       CGT     7       7       1
    ...
It'll then add DMR information to bases that occur in the DMR bed file as identified by Luo et al, Science 2017:
    chr8    1121061 1121131
    chr8    1178799 1178862
    chr8    1809865 1809878
    chr8    2747427 2747449
    chr8    3005866 3006040
    ...
'''

def importDMR(DMR_file):
    DMRdict = {}
    #We want to import the DMR bed file into a DMRDict
    with open(DMR_file) as f_DMR:
        for line in f_DMR:
            (chr, start, end) = line.split()
            DMR_name = "DMR:" + chr + ':' + start + '-' + end
            if chr.replace('chr','') not in DMRdict.keys():
                DMRdict[chr.replace('chr','')] = {}
            for pos in range(int(start), int(end) + 1): #We want to iterate through all base positions within DMR range
                DMRdict[chr.replace('chr','')][pos] = DMR_name
    return DMRdict

def filterDMR(tsv_file, DMRdict, output_file):
    tsv_output = open(output_file, 'w')
    with open(tsv_file) as f_tsv:
        for i, line in enumerate(tqdm(f_tsv)):
            (chr, pos) = line.split()[0:2]
            if int(pos) in DMRdict[chr.replace('chr','')].keys():
                tsv_output.write(line.rstrip() + "\t" + DMRdict[chr.replace('chr','')][int(pos)] + "\n")
            else:
                tsv_output.write(line) #We want to keep the CpG information line in the file, but without label of DMR.  This allows us to have options for downstream analysis
    return

def main():
    parser = argparse.ArgumentParser(description="Filter methlyation signal found in DMR regions")
    parser.add_argument('--tsv', action="store", dest="tsv_file", help="Original methylpy tsv file for all CGN")
    parser.add_argument('--DMR', action="store", dest="DMR_file", help="Bed-formatted file containing genomic regions identified as DRMs")
    parser.add_argument('--output', action="store", dest="output_file", help="Filename for output filtered tsv file")
    args = parser.parse_args()

    DMRdict = importDMR(args.DMR_file)
    filterDMR(args.tsv_file, DMRdict, args.output_file)

if __name__ == "__main__":
    main()
