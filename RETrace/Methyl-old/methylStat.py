#!/usr/bin/env python3
import argparse
import subprocess
import pybedtools
import itertools
from Bio import SeqIO
from Bio.Seq import Seq

'''
Usage: python script.py input1.mpileup input2.mpileup ... --ref ref.fa
This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039>
To run this script, we must include the following input:
    1) input.bam = list of bam files containing Bismark mapped reads.  This will be used to either generate the mpileup file (to calculate raw number of CpGs covered) or imported in bed format for CGI counting
        samtools mpileup -f <genome.fa> <read1_val_1.sort.bam> > <read1_val_1.pileup>
    2) ref.fa = reference fasta file containing sequences to which reads were mapped
    3) cpgIslandExt.bed = reference bed file containing CGI locations across the desired reference genome
The script will then go through and output the following statistics/files:
    1) covStats.txt = contains basic statistics of read coverage (i.e. num CpG/CHG/CHH read, avg read coverage).  It will also include CpG island covereage statistics.
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

def parsePileup(methDict, ref_fa, refDict):
    for file_name in sorted(methDict.keys()):
        #We first want to convert each bam file into mpileup format
        bam_name = file_name + ".bam"
        pileup_name = file_name + ".pileup"
        mpileup_command = "samtools mpileup -f"+ref_fa+" "+bam_name+" >"+pileup_name
        subprocess.call(mpileup_command, shell=True)

        #Next, we will parse the pileup fill to determine the relevant bases (i.e. C/G's within CpG/CHG/CHH pairs) along with read count and methylation levels
        output = open(file_name + '.methCov.txt', 'w')
        output.write("#Chr\tPos\tRef\tChain\tTotal\tMeth\tUnMeth\tMethRate\tRef_context\tType\n")
        f_pileup = open(pileup_name)
        for line in f_pileup:
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
        f_pileup.close()
    return methDict

def CGIstats(methDict, ref_CGI):
    for file_name in sorted(methDict.keys()):
        methDict[file_name]["CGI"] = {}
        sample_bed = pybedtools.BedTool(file_name+".bam")
        ref_bed = pybedtools.BedTool(ref_CGI)
        CGI_intersect = ref_bed.intersect(sample_bed, bed=True, wa=True, wb=True)
        for CGI in CGI_intersect:
            CGI_loc = CGI.fields[0] + ":" + CGI.fields[1] + "-" + CGI.fields[2]
            read_name = CGI.fields[7]
            if CGI_loc not in methDict[file_name]["CGI"].keys():
                methDict[file_name]["CGI"][CGI_loc] = [read_name]
            else:
                methDict[file_name]["CGI"][CGI_loc].append(read_name)
    return methDict

def calcCovStats(methDict, file_name, cutoff):
    (num_CpG, total_cov) = (0, 0)
    for base_loc in sorted(methDict[file_name].keys()):
        if base_loc is "CGI":
            continue
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
        + "Unique CpG (10x)\tMean Coverage (10x)\t"
        + "Number CGI\tMean Coverage\n")
    for file_name in sorted(methDict.keys()):
        #We first want to calculate the total read/ctype statistics
        (total_CpG, total_CHG, total_CHH) = (0,0,0)
        for base_loc in sorted(methDict[file_name].keys()):
            if base_loc is "CGI": #Calculage CGI stats
                num_CGI = len(methDict[file_name]["CGI"].keys())
                if (num_CGI > 0):
                    mean_CGI_cov = round(len(list(itertools.chain.from_iterable(methDict[file_name]["CGI"].values())))/num_CGI)
                else:
                    mean_CGI_cov = "NA"
                continue
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
        stats_output.write(str(num_CGI) + "\t" + str(mean_CGI_cov) + "\n")
    stats_output.close()

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--ref', action="store", dest="ref_fa", help="Genome reference fasta location")
    parser.add_argument('--prefix', action="store", dest="prefix", default="methCov", help="Specifies prefix for output files")
    parser.add_argument('--ref_CGI', action="store", dest="ref_CGI", help="Bed file containing reference genome CGI locations")
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+', help="List of mpileup files")
    args = parser.parse_args()

    methDict = {}
    for f in args.file:
        name_list = f.name.split('/')[-1].split('.')[:-1]
        file_name = ".".join(name_list)
        methDict[file_name] = {}

    #Import reference fasta file
    refDict = parseRef(args.ref_fa)

    #Import pileup file
    methDict = parsePileup(methDict, args.ref_fa, refDict)

    #Determine number of CGI's covered within reads
    methDict = CGIstats(methDict,args.ref_CGI)

    #Calculate raw statistics for all files
    calcStats(methDict, args.prefix)

#%%
if __name__ == "__main__":
    main()
