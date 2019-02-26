#!/usr/bin/env python3
import argparse
import os
import subprocess
import pickle
from Bio import SeqIO #Remove when running on TSCC (does not have Bio module)
from Bio.Seq import Seq #Remove when running on TSCC (does not have Bio module)

'''
Usage: python script.py --cell_files ... --ref_fa ... --ref_CGI ... --methDict ...
This script will import methlylation information for all given cell types or single cells and export methDict.
The structure of typeDict is as follows:
    methDict
        file_name
            "base"
                base_loc
                    "Total" = total number reads
                    "Meth" = number of methylated reads
                    "Unmeth" = number of unmethylated reads
                    "Type" = C type (only save CpG, ignore CHH and CHG)
'''

def parseRef(ref_fa, ref_pkl):
    refDict = {}
    with open(ref_fa, "r") as ref_file:
        for chr in SeqIO.parse(ref_file, "fasta"):
            refDict[str(chr.id)] = str(chr.seq)
    with open(ref_pkl,'wb') as refDict_file:
        pickle.dump(refDict, refDict_file, protocol=pickle.HIGHEST_PROTOCOL)
    return refDict

def findctype(ref_bases):
    H_BASES = ("A", "T", "C")
    if "N" in ref_bases:
        ctype = "NA"
        return ctype
    if (len(ref_bases) == 2):
        if (ref_bases[1] == "G"):
            ctype = "CGN"
    else:
        if (ref_bases[1] in H_BASES and ref_bases[2] == "G"):
            ctype = "CHG"
        elif (ref_bases[1] in H_BASES and ref_bases[2] in H_BASES):
            ctype = "CHH"
        elif (ref_bases.startswith("CG")):
            ctype = "CGN" #Indicates CpG context
    return ctype

def parseMethCall(methDict, ref_fa, refDict):
    for file_name in sorted(methDict.keys()):
        print(file_name)
        methDict[file_name]["base"] = {}
        pileup_name = file_name + ".pileup"
        methCall_name = file_name + ".methCall.tsv"
        #Check whether we have the methCall.tsv file within directory
        if not os.path.isfile('./' + methCall_name):
            #Check whether we already have pileup within directory
            if not os.path.isfile('./' + pileup_name):
                mpileup_command = "samtools mpileup -f" + ref_fa + " " + methDict[file_name]["file_loc"] + ">" + pileup_name
                subprocess.call(mpileup_command, shell=True)
            #We will parse the pieup file to determient he relevant bases (i.e. C/G's within CpG context) alogn with read count and methylation levels
            f_methCall = open(methCall_name, 'w')
            f_methCall.write("#Chr\tPos\tChain\tType\tMeth\tTotal\n")
            f_pileup = open(pileup_name)
            for line in f_pileup:
                if len(line.split) < 6:
                    continue
                (chr, pos, ref, depth, bases, quality) = line.split()
                if (ref.upper() not in ("C", "G") or
                    "random" in chr or
                    chr == "chrM" or
                    int(depth) == 0):
                    continue
                tmp_pos = int(pos) - 1
                if ("C" in ref.upper()): #We want to make case-insensitive string comparison
                    chain = "+" #We are on the '+' strand if looking at C's
                    meth = bases.upper().count('.')
                    unmeth = bases.count('T') #Methylated C's will remain C, unmethylated will be BS converted to T
                    total_bases = meth + unmeth
                    if (total_bases == 0):
                        continue
                    #We next want to determine the C context (i.e. in CpG context) depending on reference bases at given C position
                    ref_bases = refDict[chr][tmp_pos:tmp_pos+3].upper()
                    ctype = findctype(ref_bases)
                elif ("G" in ref.upper()):
                    chain = "-" #We are on the '-' strand if looking at G's
                    meth = bases.upper().count(',')
                    unmeth = bases.count('a') #Methylated C's (G on opposite strand) will be BS converted to A
                    total_bases = meth + unmeth
                    if (total_bases == 0):
                        continue
                    #We want to determine the C context
                    if (tmp_pos == 0):
                        continue
                    elif (tmp_pos == 1):
                        ref_bases = refDict[chr][0:2].upper()
                    elif (tmp_pos > 1):
                        tmp_pos = tmp_pos - 2
                        ref_bases = refDict[chr][tmp_pos:tmp_pos + 3].upper()
                    #In order to determien C context, we must reverse complement the ref_bases
                    ref_bases = str(Seq(ref_bases).reverse_complement())
                    ctype = findctype(ref_bases)
                if (ctype == "CGN"): #Only save CpG (ignore CHH and CHG)
                    f_methCall.write(chr.replace("chr","") + "\t" + pos + "\t" + chain + "\t" + ctype + "\t" + str(meth) + "\t" + str(total_bases) + "\n")
            f_methCall.close()
            f_pileup.close()
        f_methCall = open(methCall_name)
        for line in f_methCall:
            if line.startswith('#'):
                continue
            else:
                (chr, pos, chain, ctype, meth, total_bases) = line.split()
                base_loc = chr + ":" + pos + ";" + chain
                methDict[file_name]["base"][base_loc] = {}
                methDict[file_name]["base"][base_loc]["Total"] = int(total_bases)
                methDict[file_name]["base"][base_loc]["Meth"] = int(meth)
                methDict[file_name]["base"][base_loc]["Type"] = ctype
        f_methCall.close()
    return methDict

def main():
    parser = argparse.ArgumentParser(description="Import methylation information for cell types or single cells and export methDict")
    parser.add_argument('--cell_files', action="store", dest="samples", help="Tab-delimited file containing the sample/cell type bam name/locations along with sample/cell type names")
    parser.add_argument('--ref_fa', action="store", dest="ref_fa", help="[Optional: if not precomputed dict] Genome reference fasta location")
    parser.add_argument('--refDict', action="store", dest="ref_pkl", help="Genome reference dictionary")
    parser.add_argument('--methDict', action="store", dest="methDict_file", help="File name containing exported methDict (*.pkl)")
    args = parser.parse_args()

    #Import reference fasta file
    if args.ref_fa:
        refDict = parseRef(args.ref_fa, args.ref_pkl)
    elif args.ref_pkl:
        with open(args.ref_pkl, 'rb') as refDict_file:
            refDict = pickle.load(refDict_file)
    else:
        print("Error: Specify ref_fa or ref_pkl for computing refDict")
        return

    #Import methylation calls into methDict
    methDict = {}
    with open(args.samples) as f_sample:
        for sample in f_sample:
            (file_loc, sample_name) = sample.split()
            methDict[sample_name] = {}
            methDict[sample_name]["file_loc"] = file_loc
    methDict = parseMethCall(methDict, args.ref_fa, refDict)

    with open(args.methDict_file, 'wb') as methDict_file:
        pickle.dump(methDict, methDict_file, protocol=pickle.HIGHEST_PROTOCOL)

#%%
if __name__ == "__main__":
    main()
