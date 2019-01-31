#!/usr/bin/env python3
import argparse
from skbio import DistanceMatrix
from skbio.tree import nj

'''
Usage: python script.py input1.vcf input2.vcf... --dist [Abs/NormAbs] --output distMatrix.txt
This script will take as input multiple vcf files: (delimited by spaces)
    input1.vcf input2.vcf ...
Each of these vcf files will have the following format for each chromosomal target captured:
    #CHROM  POS_ID  REF ALT QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 ...
It will then parse these vcf files and extract the allelotype/msCounts from each of the samples contained within each vcf in order to perform a pairwise comparison of distance between each sample.
The distance metrics used can be either: (Chapal-Ilani et al, 2013; Biezuner et al, 2016)
    1) Abs [Default]: the distance between two samples, i and j, is the average absolute difference between their number of repeats in all alleles which were analyzed in both samples:
        D(i,j) = 1/L sum(abs(Ai,l - Aj,l)) for all l in L
    2) NormAbs: the distance between two samples, i and j, is the average normalized absolute distance between their number of repeats in all alleles which were analyzed in both samples:
        D(i,j) = 1/L sum(Ai,l/sum(abs(Ai,l)) - Aj,l/sum(abs(Aj,l))) for all l in L
'''

def parseVCF(file_list):
    alleleDict = {}
    for f in file_list:
        for line in f:
            if not line.startswith('##'):
                if line.startswith('#'): #This contains the vcf fields info (comes before any allelotyping calls within vcf)
                    info_line = line.split()
                    sample_names = info_line[9:] #Contains sample names from HipSTR analysis
                else:
                    vcf_line = line.split()
                    target_id = vcf_line[2] #This contains targetID (eg. "chr1:14290254-14290298_22xTG")
                    sample_info = vcf_line[9:]
                    for indx, val in enumerate(sample_info):
                        sample_id = sample_names[indx]
                        if val is not '.':
                            allelotype = val.split(':')[1]
                            if sample_id not in alleleDict.keys():
                                alleleDict[sample_id] = {}
                            alleleDict[sample_id][target_id] = [int(i) for i in allelotype.split('|')]
    return alleleDict

def calcDist(alleleDict,dist_metric,output_file):
    output = open(output_file, 'w')
    distMatrix = []
    if dist_metric == "Abs":
        for sample1 in sorted(alleleDict.keys()):
            sample1_dist = []
            for sample2 in sorted(alleleDict.keys()):
                sum_diff = 0
                target_list = set(alleleDict[sample1].keys()).intersection(set(alleleDict[sample2].keys()))
                for target_id in sorted(target_list):
                    sum_diff += abs(alleleDict[sample1][target_id][0] - alleleDict[sample2][target_id][0]) + abs(alleleDict[sample1][target_id][1] - alleleDict[sample2][target_id][1])
                abs_diff = float(sum_diff/len(target_list))
                sample1_dist.append(abs_diff)
                output.write(sample1 + "\t" + sample2 + "\t" + str(abs_diff) + "\t" + str(len(target_list)) + "\n")
            distMatrix.append(sample1_dist)
    distObj = DistanceMatrix(distMatrix,sorted(alleleDict.keys()))
    NJTree = nj(distObj)
    NJNewick = nj(distObj, result_constructor=str)

    output.write(NJTree.ascii_art() + "\n")
    output.write(NJNewick)
    output.close()

    return

def main():
    parser = argparse.ArgumentParser(description="Parse HipSTR vcf output and calculate pairwise distance of all samples")
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    parser.add_argument('--dist', action="store", dest="dist_metric", help="Specify distance metric for pairwise comparisons", default="Abs")
    parser.add_argument('--output', action="store", dest="output", help="Specify output file containing pairwise distance calculations")
    args = parser.parse_args()

    #Import allelotype for all samples found within vcf files
    alleleDict = parseVCF(args.file)

    #Calculate distance
    calcDist(alleleDict, args.dist_metric, args.output)

if __name__ == "__main__":
    main()
