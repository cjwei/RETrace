#!/usr/bin/env python3
import argparse
import subprocess
import pybedtools
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
import os
import pandas
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

'''
Usage: python script.py --sample_bams ... --cell_bams ... --ref_fa ... --prefix ...
This script is based off of SingleC_MetLevel.pl from Guo 2015 Nat Prot paper <https://doi.org/10.1038/nprot.2015.039> and Paiwise Dissimilarity calculations from Hui 2018 Stem Cell Reports paper <https://doi.org/10.1016/j.stemcr.2018.07.003>
To run this script, we must include the following input:
    1) sample_bams = tab-delimited file containing sample file locations and names
    2) cell_bams = tab-delimited file containing reference cell-type file locations and names
    2) ref.fa = reference fasta file containing sequences to which reads were mapped
    3) cpgIslandExt.bed = reference bed file containing CGI locations across the desired reference genome
The script will then go through and output the following statistics/files:
    1) input.methCov.txt = for each CpG position covered, output the following information:
        Chr, Pos, Ref, Chain, Total Reads, Meth, UnMeth, MethRate, Ref_context, Type (eg. CpG, CHH, CHG)
    2) prefix.covStats.txt = contains basic statistics of read coverage (i.e. num CpG/CHG/CHH read, avg read coverage).  It will also include CpG island covereage statistics.
    3) prefix.PD.txt = summarizes pairwise dissimilarity scores between samples and given reference cell-types
'''

def parseRef(ref_fa):
    refDict = {}
    with open(ref_fa, "r") as ref_file:
        for chr in SeqIO.parse(ref_file, "fasta"):
            refDict[str(chr.id)] = str(chr.seq)
    return refDict

def findctype(ref_bases):
    H_BASES = ("A", "T", "C")
    if "N" in ref_bases:
        ctype = "NA"
        return ctype
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
        methDict[file_name]["base"] = {}
        pileup_name = file_name + ".pileup"
        #Check whether we already have pileup file within directory
        if not os.path.isfile('./' + pileup_name):
            mpileup_command = "samtools mpileup -f" + ref_fa + " " + methDict[file_name]["file_loc"] + " >" + pileup_name
            subprocess.call(mpileup_command, shell=True)
        #We will parse the pileup fill to determine the relevant bases (i.e. C/G's within CpG/CHG/CHH pairs) along with read count and methylation levels
        output = open(file_name + '.methCall.txt', 'w')
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
                methDict[file_name]["base"][base_loc] = {}
                methDict[file_name]["base"][base_loc]["Total"] = total_bases
                methDict[file_name]["base"][base_loc]["Meth"] = meth
                methDict[file_name]["base"][base_loc]["UnMeth"] = unmeth
                methDict[file_name]["base"][base_loc]["Ref_context"] = ref_bases
                methDict[file_name]["base"][base_loc]["Type"] = ctype
        output.close()
        f_pileup.close()
    return methDict

def CGIstats(methDict, ref_CGI):
    for file_name in sorted(methDict.keys()):
        methDict[file_name]["CGI"] = {}
        sample_bed = pybedtools.BedTool(methDict[file_name]["file_loc"])
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
    for base_loc in sorted(methDict[file_name]["base"].keys()):
        if (methDict[file_name]["base"][base_loc]["Total"] >= cutoff) and methDict[file_name]["base"][base_loc]["Type"] == "CpG":
            num_CpG += 1
            total_cov += methDict[file_name]["base"][base_loc]["Total"]
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
        #Summarize CGI stats
        num_CGI = len(methDict[file_name]["CGI"].keys())
        if (num_CGI > 0):
            mean_CGI_cov = round(len(list(itertools.chain.from_iterable(methDict[file_name]["CGI"].values())))/num_CGI)
        else:
            mean_CGI_cov = "NA"
        #Determine C coverage
        for base_loc in sorted(methDict[file_name]["base"].keys()):
            if (methDict[file_name]["base"][base_loc]["Type"] == "CpG"):
                total_CpG += 1
            elif (methDict[file_name]["base"][base_loc]["Type"] == "CHG"):
                total_CHG += 1
            elif (methDict[file_name]["base"][base_loc]["Type"] == "CHH"):
                total_CHH += 1
        stats_output.write(file_name + "\t" + str(total_CpG) + "/" + str(total_CHG) + "/" + str(total_CHH) + "\t")
        for cutoff in [1,5,10]:
            (num_CpG, mean_cov) = calcCovStats(methDict, file_name, cutoff)
            stats_output.write(str(num_CpG) + "\t" + str(mean_cov) + "\t")
        stats_output.write(str(num_CGI) + "\t" + str(mean_CGI_cov) + "\n")
    stats_output.close()
    return

def calcPD(sampleDict, typeDict, seqDepth, prefix):
    PDdict = {} #We want to save all paiwise_dis as dictionary where we have ordered lists for each file analyzed (ordered alphabetically by filename)
    PD_output = open(prefix + ".PD.txt", 'w')
    PD_output.write("Sample\tCell Type\tPairwise Dissimilarity\tNum Shared\n")
    for sample_name in sorted(sampleDict.keys()):
        PDdict[sample_name] = []
        for type_name in sorted(typeDict.keys()):
            dis_sum = 0
            num_shared = 0
            for shared_base in sampleDict[sample_name]["base"].keys() & typeDict[type_name]["base"].keys():
                if sampleDict[sample_name]["base"][shared_base]["Type"] == "CpG" and typeDict[type_name]["base"][shared_base]["Type"] == "CpG":
                    if sampleDict[sample_name]["base"][shared_base]["Total"] >= seqDepth and typeDict[type_name]["base"][shared_base]["Total"] >= seqDepth:
                        methRate1 = float(sampleDict[sample_name]["base"][shared_base]["Meth"]/sampleDict[sample_name]["base"][shared_base]["Total"])
                        methRate2 = float(typeDict[type_name]["base"][shared_base]["Meth"]/typeDict[type_name]["base"][shared_base]["Total"])
                        if methRate1.is_integer() and methRate2.is_integer():
                            num_shared += 1
                            if methRate1 == methRate2:
                                dis_sum += 0
                            else:
                                dis_sum += 100
            pairwise_dis = float(dis_sum/num_shared)
            PD_output.write(sample_name + "\t" + type_name + "\t" + str(round(pairwise_dis,4)) + "\t" + str(num_shared) + "\n")
            PDdict[sample_name].append(pairwise_dis)
    PD_output.close()

    #Print clustermap using the calculated pairwise disimilarity
    PDdf = pandas.DataFrame(PDdict, index=sorted(typeDict.keys()))
    PD_clustermap = sns.clustermap(PDdf)
    PD_clustermap.savefig(prefix + ".png")

def main():
    parser = argparse.ArgumentParser(description="Calculate methylation coverage")
    parser.add_argument('--sample_bams', action="store", dest="samples", help="Tab-delimited file containing the sample bam name/locations along with sample names")
    parser.add_argument('--cell_bams', action="store", dest="cell_types", help="Tab-delimited file containing the reference cell-type specific bam files along with the known cell-type")
    parser.add_argument('--ref_fa', action="store", dest="ref_fa", help="Genome reference fasta location")
    parser.add_argument('--prefix', action="store", dest="prefix", default="methCov", help="Specifies prefix for output files")
    parser.add_argument('-stats', action="store_true", help="Optional: Output statistics including unique CpG/CGI counts")
    parser.add_argument('--ref_CGI', action="store", dest="ref_CGI", help="Bed file containing reference genome CGI locations")
    parser.add_argument('--seqDepth', action="store", dest="seqDepth", default=1, help="Minimum number reads covering CpG for calculating pd")
    args = parser.parse_args()

    #Import reference fasta file
    refDict = parseRef(args.ref_fa)

    #Import methylation calls for sample files
    sampleDict = {}
    with open(args.samples) as f_sample:
        for sample in f_sample:
            (file_loc, sample_name) = sample.split()
            sampleDict[sample_name] = {}
            sampleDict[sample_name]["file_loc"] = file_loc
    sampleDict = parsePileup(sampleDict, args.ref_fa, refDict)

    #Import methylation calls for reference cell-type files
    typeDict = {}
    with open(args.cell_types) as f_cell:
        for cell_type in f_cell:
            (file_loc, type_name) = cell_type.split()
            typeDict[type_name] = {}
            typeDict[type_name]["file_loc"] = file_loc
    typeDict = parsePileup(typeDict, args.ref_fa, refDict)

    refDict.clear() #Remove refDict in order to clear up memory

    #Caculate methylation coverage statistics
    if args.stats is True:
        if args.ref_CGI:
            sampleDict = CGIstats(sampleDict, args.ref_CGI)
        else:
            parser.error('Cannot calculate statistics without specifying ref_CGI file')

    #Calculate raw statistics for all files
    if args.stats is True:
        calcStats(sampleDict, args.prefix)

    #Calculate pairwise dissimilarity matrix
    calcPD(sampleDict, typeDict, int(args.seqDepth), args.prefix)

#%%
if __name__ == "__main__":
    main()
