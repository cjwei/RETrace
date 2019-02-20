#!/usr/bin/env python3
import argparse
import os
import subprocess
from skbio import DistanceMatrix
from skbio.tree import nj

'''Usage: python script.py --input sample_list.txt --vcf output.vcf
This script will take as input a sample list file that contains the following information:
    1) File location of sample bam
    2) Sample name
    3) [Optional] Sample type (i.e. known clone in ex vivo tree)
It will then perform read group assignment and run HipSTR on the given files.  We will then parse the HipSTR vcf file output and extract the allelotype/msCounts from each given sample.
The HipSTR vcf file have the followign format for each microsatellite target captured:
    #CHROM  POS_ID  REF ALT QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 ...
'''

def run_HipSTR(sampleDict, vcf_output, prefix):
    '''Wrapper command to add readGroup to bam files and run HipSTR'''
    readGroup_list = []
    for sample in sorted(sampleDict.keys()):
        bam_readGroup = '.'.join(sampleDict[sample]["bam"].split('.')[:-1]) + ".readGroup.bam"
        readGroup_list.append(bam_readGroup)
        sampleDict[sample]["bam_readGroup"] = bam_readGroup
        if not os.path.isfile(bam_readGroup):
            readGroup_command = "java -jar /home/cjwei/software/picard.jar AddOrReplaceReadGroups " + \
                "I=" + sampleDict[sample]["bam"] + " " + "O=" + bam_readGroup + " " + \
                "RGLB=" + sample + " " + "RGPL=illumina RGPU=unit1 RGSM=" + sample
            index_command = "samtools index " + bam_readGroup
            subprocess.call(readGroup_command, shell=True)
            subprocess.call(index_command, shell=True)
    hipSTR_command = "/home/cjwei/software/HipSTR/HipSTR --bams " + ",".join(sorted(readGroup_list)) + " " + \
        "--fasta /media/6TB_slot4/GenomeDB/hg19/raw_fasta/hg19_reference.fa " + \
        "--regions /home/cjwei/software/RETrace/MS/HipSTR/20171211_probes.bed " + \
        "--str-vcf " + vcf_output + ".gz --log " + prefix + ".log --min-reads 15 --use-unpaired --no-rmdup"
    print(hipSTR_command)
    gunzip_command = "gunzip " + vcf_output
    subprocess.call(hipSTR_command, shell=True)
    subprocess.call(gunzip_command, shell=True)
    return

def parseVCF(sampleDict, vcf_output, min_qual, min_reads, max_stutter):
    with open(vcf_output) as vcf:
        for line in vcf:
            if not line.startswith('##'):
                if line.startswith('#'): #This contains the vcf fields info (comes before any allelotyping calls within vcf)
                    info_line = line.split()
                    sample_names = info_line[9:] #Contains sample names from HipSTR analysis
                    for sample in sample_names:
                        if sample in sampleDict.keys():
                            sampleDict[sample]["HipSTR"] = {}
                        else:
                            print("Mismatched sample names in vcf file (rerun HipSTR)")
                            return
                else:
                    vcf_line = line.split()
                    target_id = vcf_line[2] #This contains targetID (eg. "chr1:14290254-14290298_22xTG")
                    sample_info = vcf_line[9:]
                    for indx, val in enumerate(sample_info):
                        sample = sample_names[indx]
                        if val is not '.':
                            allelotype = val.split(':')[1] #GB (base pair differences of genotype from reference)
                            prob_genotype = val.split(':')[2] #Q (posteriior probability of unphased genotype)
                            num_reads = val.split(':')[4] #DP (number of reads used for sample's genotype)
                            num_stutter = val.split(':')[6] #DSTUTTER (number of reads with a stutter indel in the STR region)
                            if float(prob_genotype) >= min_qual and int(num_reads) >= min_reads and float(int(num_stutter)/int(num_reads)) <= max_stutter:
                                sampleDict[sample]["HipSTR"][target_id] = {}
                                sampleDict[sample]["HipSTR"][target_id]["allelotype"] = [int(i) for i in allelotype.split('|')]
                                sampleDict[sample]["HipSTR"][target_id]["stats"] = [float(prob_genotype),int(num_reads),int(num_stutter)] #This will contain important stats for HipSTR genotype call
    return sampleDict

def calcDist(sampleDict, target_file, dist_metric, verbose, prefix):
    distDict = {}
    distDict["sampleComp"] = {}
    distDict["targetComp"] = {}
    #We want to create a list of all targets to be used for distance calculation (specified by user)
    if target_file != "All":
        with open(target_file) as f:
            specified_targets = f.read().splitlines()
    if dist_metric == "Abs":
        for sample1 in sorted(sampleDict.keys()):
            distDict["sampleComp"][sample1] = {}
            for sample2 in sorted(sampleDict.keys()):
                distDict["sampleComp"][sample1][sample2] = {}
                sum_diff = 0
                target_list = set(sampleDict[sample1]["HipSTR"].keys()).intersection(set(sampleDict[sample2]["HipSTR"].keys()))
                num_targets = 0
                for target_id in sorted(target_list):
                    try:
                        if target_id not in specified_targets:
                            continue
                    except:
                        pass
                    num_targets += 1
                    genotype_diff = abs(sampleDict[sample1]["HipSTR"][target_id]["allelotype"][0] - sampleDict[sample2]["HipSTR"][target_id]["allelotype"][0]) + abs(sampleDict[sample1]["HipSTR"][target_id]["allelotype"][1] - sampleDict[sample2]["HipSTR"][target_id]["allelotype"][1])
                    sum_diff += genotype_diff
                    if target_id not in distDict["targetComp"].keys():
                        distDict["targetComp"][target_id] = {}
                    distDict["targetComp"][target_id][tuple(sorted([sample1,sample2]))] = genotype_diff
                abs_diff = float(sum_diff/len(target_list))
                distDict["sampleComp"][sample1][sample2]["dist"] = abs_diff
                distDict["sampleComp"][sample1][sample2]["num_targets"] = num_targets
    if verbose is True: #We want to determine useful targets
        targetOutput = open(prefix + ".stats.out", 'w')
        targetOutput.write("targetID\tIntra-clone Dist\tNum Intra-clone Pairs\tInter-clone Dist\tNum Inter-clone Pairs\tDist Bool\t" + "\t".join(sorted(sampleDict.keys())) + "\n")
        for target_id in sorted(distDict["targetComp"].keys()):
            #Calculate average intra (within) and inter (across) clone distances
            (intra_dist, inter_dist, num_intra, num_inter, avg_intra_dist, avg_inter_dist) = (0,0,0,0,0,0) #Declare value of -1 for undefined
            print("----------" + target_id + "----------")
            for (sample1,sample2) in distDict["targetComp"][target_id].keys():
                print("\t".join([sample1,sample2]) + "\t" + str(distDict["targetComp"][target_id][(sample1,sample2)]) + "\t" + '|'.join(str(i) for i in sampleDict[sample1]["HipSTR"][target_id]["allelotype"]) + "\t" + '|'.join(str(i) for i in sampleDict[sample2]["HipSTR"][target_id]["allelotype"]))
                if sampleDict[sample1]["clone"] == sampleDict[sample2]["clone"] and sample1 != sample2:
                    intra_dist += distDict["targetComp"][target_id][(sample1,sample2)]
                    num_intra += 1
                elif sampleDict[sample1]["clone"] != sampleDict[sample2]["clone"] and sampleDict[sample1]["clone"] != sampleDict[sample2]["clone"]:
                    inter_dist += distDict["targetComp"][target_id][(sample1,sample2)]
                    num_inter += 1
            if num_intra > 0 and num_inter > 0:
                avg_intra_dist = float(intra_dist/num_intra) #Average allelotype difference within group
                avg_inter_dist = float(inter_dist/num_inter) #Average allelotype difference across groups
                if avg_intra_dist > avg_inter_dist:
                    diff_bool = "intra"
                elif avg_intra_dist < avg_inter_dist:
                    diff_bool = "inter" #We want inter_dist (across groups) to be larger than intra_dist
                elif avg_intra_dist == avg_inter_dist:
                    diff_bool = "eq"
            else:
                diff_bool = "NA"
            #We want to save the allotypes for all samples
            allelotype_list = []
            for sample in sorted(sampleDict.keys()):
                if target_id in sampleDict[sample]["HipSTR"].keys():
                    allelotype_info = '|'.join(str(i) for i in sampleDict[sample]["HipSTR"][target_id]["allelotype"]) + "(" + ','.join(str(j) for j in sampleDict[sample]["HipSTR"][target_id]["stats"]) + ")"
                    allelotype_list.append(allelotype_info)
                else:
                    allelotype_list.append('.')
            targetOutput.write(target_id + "\t" + str(round(avg_intra_dist,2)) + "\t" + str(num_intra) + "\t" + str(round(avg_inter_dist,2)) + "\t" + str(num_inter) + "\t" + diff_bool + "\t" + "\t".join(allelotype_list) + "\n")
        targetOutput.close()
    return distDict

def drawTree(distDict, sampleDict, prefix):
    tree_output = open(prefix + ".tree.out", 'w')
    distMatrix = []
    targetMatrix = []
    for sample1 in sorted(sampleDict.keys()):
        sample1_dist = []
        sample1_targets = []
        for sample2 in sorted(sampleDict.keys()):
            sample1_dist.append(distDict["sampleComp"][sample1][sample2]["dist"])
            sample1_targets.append(distDict["sampleComp"][sample1][sample2]["num_targets"])
        distMatrix.append(sample1_dist)
        targetMatrix.append(sample1_targets)
    for dist_indx,dist_list in enumerate(distMatrix): #Print matrix containing distances
        tree_output.write(sorted(sampleDict.keys())[dist_indx] + "," + ",".join(str(round(i,3)) for i in dist_list) + "\n")
    for target_indx,target_list in enumerate(targetMatrix): #Print matrix containing number targets shared between each pair
        tree_output.write(sorted(sampleDict.keys())[target_indx] + "," + ",".join(str(j) for j in target_list) + "\n")
    distObj = DistanceMatrix(distMatrix,sorted(sampleDict.keys()))
    NJTree = nj(distObj)
    NJNewick = nj(distObj, result_constructor=str)
    tree_output.write(NJTree.ascii_art() + "\n")
    tree_output.write(NJNewick + "\n")
    tree_output.close()
    return

def main():
    parser = argparse.ArgumentParser(description="Run and analyze HipSTR output to determine distance between single cells")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specify output file containing pairwise distance calculations")
    parser.add_argument('--targets', action="store", dest="target_file", default="All", help="[Optional] Specify list of targets used for distance calculation")
    parser.add_argument('--dist', action="store", dest="dist_metric", help="Specify distance metric for pairwise comparisons", default="Abs")
    parser.add_argument('--vcf', action="store", dest="vcf_output")
    parser.add_argument('--min-call-qual', action="store", dest="min_qual", default=0, help="Specify the minimum posterior probability of genotype for filtering")
    parser.add_argument('--min-reads', action="store", dest="min_reads", default=1, help="Cutoff for minimum number of reads required for calling allelotype")
    parser.add_argument('--max-stutter', action="store", dest="max_stutter", default=1, help="Define maximum number of reads that can be classified as stutter")
    parser.add_argument('-v', action="store_true", help="Flag for determining whether we want to output all statistics for shared targets in output")
    args = parser.parse_args()

    #We want to first parse the sample_info file
    sampleDict = {}
    with open(args.sample_info) as f:
        for line in f:
            if len(line.split()) == 3: #If clone is specified
                (bam, sample, clone) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
                sampleDict[sample]["clone"] = clone
            elif len(line.split()) == 2: #If clone is not specified
                (bam, sampleName) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
            else:
                print("Incorrect formmating for line (bam, sampleName, [optional] clone):\n" + "\t".join(line.split))
                return

    #Run HipSTR if you haven't done so already
    if not os.path.isfile('./' + args.vcf_output):
        run_HipSTR(sampleDict, args.vcf_output, args.prefix)

    #Parse vcf file to import HipSTR allelotypes
    sampleDict = parseVCF(sampleDict, args.vcf_output, float(args.min_qual), int(args.min_reads), float(args.max_stutter))

    #Calculate pairwise distance
    distDict = calcDist(sampleDict, args.target_file, args.dist_metric, args.v, args.prefix)

    #Draw neighbor joining tree
    drawTree(distDict, sampleDict, args.prefix)

if __name__ == "__main__":
    main()
