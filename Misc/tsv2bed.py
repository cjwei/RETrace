#!/usr/bin/env python3
import argparse
import pybedtools
import gzip
import os

def main():
    '''
    This simple script will take as input methylpy tsv files in the following format:
        chrom   pos     strand  C_context   num_meth    num_reads   methyl_call (ignore)
    The script will then convert to bed files compatible with metilene with the following format:
        chrom   start   end     methylation_rate
    This will be utilized for metilene script in order to call DMRs
    '''
    parser = argparse.ArgumentParser(description="Convert methylpy tsv file to metilene-compatible bed file")
    parser.add_argument("--tsv", action="store", nargs='+', dest="tsv_file", help="Specify input methylpy tsv file (if specify multiple tsv, will merge prior to conversion to bed)")
    parser.add_argument("--bed", action="store", dest="bed_file", help="Specify bed file output")
    args = parser.parse_args()

    methDict = {} #We want to merge all shared CpG across single cells in given tsv files
    for tsv_file in args.tsv_file:
        if "gz" in tsv_file:
            f_tsv = gzip.open(tsv_file, 'rb')
        else:
            f_tsv = open(tsv_file, 'r')
        for line in f_tsv:
            if "gz" in tsv_file:
                line = line.decode('utf-8')
            (chrom, pos, strand, ctype, num_meth, num_reads) = line.split()[0:6]
            if ctype.startswith('CG'): #We only care about CGN methylated cytosines
                base_loc = chrom + ":" + pos
                if base_loc not in methDict.keys():
                    methDict[base_loc] = {}
                    methDict[base_loc]["num_meth"] = 0
                    methDict[base_loc]["num_reads"] = 0
                methDict[base_loc]["num_meth"] += int(num_meth)
                methDict[base_loc]["num_reads"] += int(num_reads)
        f_tsv.close()

    f_temp = open("temp.bed", 'w')
    for base_loc in sorted(methDict.keys()):
        (chrom, pos) = base_loc.split(':')
        meth_rate = float(methDict[base_loc]["num_meth"] / methDict[base_loc]["num_reads"])
        f_temp.write("chr" + chrom.replace('chr','') + "\t" + pos + "\t" + str(int(pos) + 1) + "\t" + str(meth_rate) + "\n")
    f_temp.close()

    #We want to sort the temp file
    with open(args.bed_file, 'w') as f_out:
        f_out.write(str(pybedtools.BedTool("temp.bed").sort()))
    os.remove("temp.bed")

if __name__ == "__main__":
    main()
