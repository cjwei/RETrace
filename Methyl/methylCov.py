#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

'''
Usage: python script.py input1.mpileup input2.mpileup ... --ref ref.fa
This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039>
To run this script, we must include the following input:
    1) input.mpileup = list of files containing mpileup conversion of the original mapped reads.  Generate this file using the following command:
        samtools mpileup -f <genome.fa> <read1_val_1.sort.bam> > <read1_val_1.pileup>
    2) ref.fa = reference fasta file containing sequences to which reads were mapped
The script will then go through and output the following statistics/files:
    1) covStats.txt = contains basic statistics of read coverage (i.e. num CpG/CHG/CHH read, avg read coverage)
    2) methylLevel.png = contains plots of frequency of different methylation levels in dataset (similar to Supp Fig 2 in Guo 2015 Nat Prot paper)
'''

def parseRef(ref_fa):
    refDict = {}
    with open(ref_fa, "r") as ref_file:
        for chr in SeqIO.parse(ref_file, "fasta"):
            refDict[str(chr.id)] = str(chr.seq)
    return refDict

def findctype(ref_bases):
    H_BASES = ("A", "T", "C")
    if (len(ref_bases) == 2):
        if (ref_bases[1] == "G"):
            ctype = "CpG"
    else:
        if (ref_bases[1] in H_BASES and ref_bases[2] == "G"):
            ctype = "CHG"
        elif (ref_bases[1] in H_BASES and ref_bases[2] in H_BASES):
            ctype = "CHH"
        elif (ref_bases.startswith("CG")):
            ctype = "CpG"
    return ctype

def parsePileup(file_list, refDict):
    methDict = {}
    for f in file_list:
        name_list = f.name.split('/')[-1].split('.')[:-1]
        file_name = ".".join(name_list)
        methDict[file_name] = {}
        output = open(file_name + '.methCov.txt', 'w')
        output.write("#Chr\tPos\tRef\tChain\tTotal\tMeth\tUnMeth\tMethRate\tRef_context\tType\n")
        for line in f:
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
            meth_rate = meth / total_bases
            if (ctype):
                output.write(chr + "\t" + pos + "\t" + ref + "\t" + chain + "\t" + str(total_bases) + "\t" + str(meth) + "\t" + str(unmeth) + "\t" + str(meth_rate) + "\t" + ref_bases + "\t" + ctype + "\n")
                #We also want to save all of these stats to methDict"
                base_loc = chr + ":" + pos + ";" + ref + ";" + chain
                methDict[file_name][base_loc] = {}
                methDict[file_name][base_loc] = {}
                methDict[file_name][base_loc]["Total"] = total_bases
                methDict[file_name][base_loc]["Meth"] = meth
                methDict[file_name][base_loc]["UnMeth"] = unmeth
                methDict[file_name][base_loc]["Ref_context"] = ref_bases
                methDict[file_name][base_loc]["Type"] = ctype
        output.close()
    return methDict

def calcCovStats(methDict, file_name, cutoff):
    (num_CpG, total_cov) = (0, 0)
    for base_loc in sorted(methDict[file_name].keys()):
        if (methDict[file_name][base_loc]["Total"] >= cutoff) and methDict[file_name][base_loc]["Type"] == "CpG":
            num_CpG += 1
            total_cov += methDict[file_name][base_loc]["Total"]
    if(num_CpG == 0):
        return (num_CpG, "NA")
    mean_cov = round(total_cov / num_CpG)
    return (num_CpG, mean_cov)

def calcStats(methDict, prefix):
    '''
    This funciton will calculate the following statistics per file:
        1) Total number of CpG/CHG/CHH positions covered
        2) Avg number of reads covering each position
    '''
    stats_output = open(prefix + ".covStats.txt", 'a+')
    stats_output.write("File\tTotal CpG/CHG/CHH\t"
        + "Unique CpG (1x)\tMean Coverage (1x)\t"
        + "Unique CpG (5x)\tMean Coverage (5x)\t"
        + "Unique CpG (10x)\tMean Coverage (10x)\n")
    for file_name in sorted(methDict.keys()):
        #We first want to calculate the total read/ctype statistics
        (total_CpG, total_CHG, total_CHH) = (0,0,0)
        for base_loc in sorted(methDict[file_name].keys()):
            if (methDict[file_name][base_loc]["Type"] == "CpG"):
                total_CpG += 1
            elif (methDict[file_name][base_loc]["Type"] == "CHG"):
                total_CHG += 1
            elif (methDict[file_name][base_loc]["Type"] == "CHH"):
                total_CHH += 1
        stats_output.write(file_name + "\t" + str(total_CpG) + "/" + str(total_CHG) + "/" + str(total_CHH) + "\t")
        for cutoff in [1,5,10]:
            (num_CpG, mean_cov) = calcCovStats(methDict, file_name, cutoff)
            stats_output.write(str(num_CpG) + "\t" + str(mean_cov) + "\t")
        stats_output.write("\n")
    stats_output.close()

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--ref', action="store", dest="ref_fa", help="Genome reference fasta location")
    parser.add_argument('--prefix', action="store", dest="prefix", default="methCov", help="Specifies prefix for output files")
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+', help="List of mpileup files")
    args = parser.parse_args()

    #Import reference fasta file
    refDict = parseRef(args.ref_fa)

    #Import pileup file
    methDict = parsePileup(args.file, refDict)

    #Calculate raw statistics for all files
    calcStats(methDict, args.prefix)

#%%
if __name__ == "__main__":
    main()
