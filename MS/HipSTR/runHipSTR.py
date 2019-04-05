#!/usr/bin/env python3
import argparse
import os
import subprocess
import collections
import more_itertools
from skbio import DistanceMatrix
from skbio.tree import nj
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
from itertools import cycle
from tqdm import tqdm
from ete3 import Tree #Call ETE toolkit <http://etetoolkit.org/docs/latest/tutorial/index.html>
import numpy as np

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
        "--str-vcf " + vcf_output + ".gz --log " + prefix + ".log --use-unpaired --no-rmdup"
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
                    #We want to parse target_info into each specific field (i.e. INFRAME_PGEOM, INFRAME_UP, DP, DSTUTTER, etc)
                    info_name = []
                    info_val = []
                    for info in vcf_line[7].split(';'):
                        info_name.append(info.split('=')[0])
                        info_val.append(info.split('=')[1])
                    sample_format = vcf_line[9:]
                    for indx, val in enumerate(sample_format):
                        sample = sample_names[indx]
                        if val is not '.':
                            allelotype = val.split(':')[1] #GB (base pair differences of genotype from reference)
                            prob_genotype = val.split(':')[2] #Q (posteriior probability of unphased genotype)
                            num_reads = val.split(':')[4] #DP (number of reads used for sample's genotype)
                            num_stutter = val.split(':')[6] #DSTUTTER (number of reads with a stutter indel in the STR region)
                            msCounts = val.split(':')[-2] #ALLREADS (bp difference observed in each read's alignment)
                            if float(prob_genotype) >= min_qual and int(num_reads) >= min_reads and float(int(num_stutter)/int(num_reads)) <= max_stutter:
                                sampleDict[sample]["HipSTR"][target_id] = {}
                                sampleDict[sample]["HipSTR"][target_id]["allelotype"] = [int(i) for i in allelotype.split('|')]
                                sampleDict[sample]["HipSTR"][target_id]["stats"] = [float(prob_genotype),int(num_reads),int(num_stutter)] #This will contain important stats for HipSTR genotype call
                                sampleDict[sample]["HipSTR"][target_id]["msCounts"] = msCounts.split(';')
                                #We want to define "info" sub-dictionary containing target info names/values
                                info_zip = zip(info_name, info_val)
                                sampleDict[sample]["HipSTR"][target_id]["info"] = dict(info_zip)
    return sampleDict

def calcBulk(sampleDict, min_freq):
    alleleDict = {}
    for sample in sorted(sampleDict.keys()):
        for target_id in sorted(sampleDict[sample]["HipSTR"].keys()):
            if target_id not in alleleDict.keys():
                alleleDict[target_id] = {}
                alleleDict[target_id]["rawCounts"] = [] #Contains merged list of all single cell allelotypes
            alleleDict[target_id]["rawCounts"].extend(sampleDict[sample]["HipSTR"][target_id]["allelotype"])
    #Determine allele frequency from merged single cell allelotypes
    for target_id in sorted(alleleDict.keys()):
        # print("-----\t" + target_id + "\t-----")
        sub_len = len(target_id.split('x')[-1])
        alleleDict[target_id]["countFreq"] = collections.Counter(alleleDict[target_id]["rawCounts"])
        filtered_alleles = []
        for allele, count in sorted(alleleDict[target_id]["countFreq"].items()):
            if float(count/sum(alleleDict[target_id]["countFreq"].values())) >= min_freq: #We want to filter out and keep only alleles that have frequency of >=0.1 across all allelotypes
                filtered_alleles.append(int(allele/sub_len)) #We only care about whole subunit mutations
        #Group all consecutive filtered_alleles (this would indicate allele1/allele2 found in merged files)
        alleleDict[target_id]["allele_groups"] = {}
        # alleleDict[target_id]["allele_groups"]["All"] = [allele * sub_len for allele in list(set(filtered_alleles))] #For no allele grouping condition
        for indx, group in enumerate(more_itertools.consecutive_groups(sorted(set(filtered_alleles)))): #Need to use "set(filtered_alleles)" in order to prevent repeated int in filtered_alleles
            alleleDict[target_id]["allele_groups"][indx] = [allele * sub_len for allele in list(group)] #Save allele groups in terms of raw number of bases difference from ref (same as HipSTR output)
        # print(alleleDict[target_id]["allele_groups"].values())
    return alleleDict

def plotVCF(sampleDict, prefix, num_plot):
    #Reformat sampleDict to use target_id as first key and sample as second key
    targetDict = {} #Dictionary for plotting VCF
    for sample in sampleDict.keys():
        for target_id in sampleDict[sample]["HipSTR"].keys():
            if not target_id in targetDict.keys():
                targetDict[target_id] = {}
            msCount_list = []
            msCount_freq = []
            for msCount_info in sampleDict[sample]["HipSTR"][target_id]["msCounts"]:
                if msCount_info is not '.':
                    [msCount,freq] = msCount_info.split('|')
                    msCount_list.append(int(msCount))
                    msCount_freq.append(int(freq))
            targetDict[target_id][sample] = {}
            targetDict[target_id][sample]["msCount_list"] = msCount_list
            targetDict[target_id][sample]["msCount_freq"] = msCount_freq
            targetDict[target_id][sample]["allelotype"] = sampleDict[sample]["HipSTR"][target_id]["allelotype"]

    #Plot MS counts
    plots = PdfPages(prefix + ".VCFplot.pdf")
    color_options = ['xkcd:pale green','xkcd:pale blue','xkcd:light grey','xkcd:pale pink']
    color_cycle = cycle(color_options)
    next_color = next(color_cycle)

    for target_id in tqdm(sorted(targetDict.keys())[0:min(num_plot,len(targetDict.keys()))]):
    # for target_id in tqdm(sorted(targetDict.keys())):
        if len(targetDict[target_id].keys()) >= 2:
            select_color, next_color = next_color, next(color_cycle)
        else:
            select_color = 'xkcd:white'
        for sample in sorted(targetDict[target_id].keys()):
            if len(targetDict[target_id][sample]["msCount_list"]) < 1:
                continue
            #Convert msCount and allelelotype to actual bases (not just base differences from reference microsatellite)
            sub_info = target_id.split('_')[-1]
            [num_sub, sub_seq] = sub_info.split('x')
            ref_len = int(num_sub) * int(len(sub_seq))
            count_list = [(int(i) + ref_len) for i in sorted(set(targetDict[target_id][sample]["msCount_list"]))]
            count_freq = [float(i/sum(targetDict[target_id][sample]["msCount_freq"])) for i in targetDict[target_id][sample]["msCount_freq"]]
            allele_list = [(int(i) + ref_len) for i in sorted(targetDict[target_id][sample]["allelotype"])]
            #Plot raw alleles/msCounts (blue) along with final allelotype (red)
            ax = plt.gca()
            ax.set_ylim([0,1])
            ax.set_xlim([min(count_list) - 5, max(count_list) + 5])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.plot(allele_list, [0.5]*len(allele_list), 'ro')
            plt.vlines(count_list, [0], count_freq, linestyle="dashed", color="b")
            plt.title(target_id + ", " + sample)
            ax.set_facecolor(select_color)
            plt.savefig(plots, format='pdf')
            plt.clf()
    plots.close()
    return

def calcDist(sampleDict, alleleDict, distDict, sample1, sample2, final_target_list, dist_metric):
    total_dist = 0
    num_alleles = 0
    # print("----------" + sample1 + "\t" + sample2 + "\t" + str(len(final_target_list)) + "----------")
    for target_id in final_target_list:
        target_dist = 0
        # print(sample1 + "\t" + sample2 + "\t" + target_id + "\t" + ','.join(str(x) for x in sampleDict[sample1]["HipSTR"][target_id]["allelotype"]) + "\t" + ','.join(str(y) for y in sampleDict[sample2]["HipSTR"][target_id]["allelotype"]))
        for allele_indx, allele_group in alleleDict[target_id]["allele_groups"].items():
            allelotype1 = list(set(sampleDict[sample1]["HipSTR"][target_id]["allelotype"]).intersection(set(allele_group)))
            allelotype2 = list(set(sampleDict[sample2]["HipSTR"][target_id]["allelotype"]).intersection(set(allele_group)))
            if len(allelotype1) == 2 and len(allelotype2) == 2:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]) + abs(allelotype1[1] - allelotype2[1]), abs(allelotype1[1] - allelotype2[0]) + abs(allelotype1[0]- allelotype2[1]))
                elif dist_metric == "EqorNot":
                    target_dist += int(len(set(allelotype1).symmetric_difference(set(allelotype2)))/2)
                num_alleles += 2
            elif len(allelotype1) == 1 and len(allelotype2) == 2:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]), abs(allelotype1[0] - allelotype2[1]))
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0] and allelotype1[0] != allelotype2[1]:
                        target_dist += 1
                num_alleles += 1
            elif len(allelotype1) == 2 and len(allelotype2) == 1:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]), abs(allelotype1[1] - allelotype2[0]))
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0] and allelotype1[1] != allelotype2[0]:
                        target_dist += 1
                num_alleles += 1
            elif len(allelotype1) == 1 and len(allelotype2) == 1:
                if dist_metric == "Abs":
                    target_dist += abs(allelotype1[0] - allelotype2[0])
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0]:
                        target_dist += 1
                num_alleles += 1
            # if target_dist > 0:
            #     print(sample1 + "\t" + sample2 + "\t" + target_id + "\t" + ','.join(str(i) for i in allele_group) + "\t" + ','.join(str(j) for j in allelotype1) + "\t" + ','.join(str(k) for k in allelotype2) + "\t" + str(target_dist))
        if target_id not in distDict["targetComp"].keys():
            distDict["targetComp"][target_id] = {}
        distDict["targetComp"][target_id][tuple(sorted([sample1,sample2]))] = target_dist #Save contribution to distance from each target_id
        total_dist += target_dist
    # print("Total Dist:\t" + str(total_dist) + "\tNum Alleles:\t" + str(num_alleles))
    distDict["sampleComp"][sample1][sample2]["dist"] = float(total_dist/num_alleles)
    distDict["sampleComp"][sample1][sample2]["num_targets"] = len(final_target_list)
    return distDict

def makeDistMatrix(sampleDict, alleleDict, target_list, dist_metric, verbose, prefix, bootstrap):
    distDict = {}
    distDict["sampleComp"] = {}
    distDict["targetComp"] = {}
    #Calculate pairwise distance for all samples depending on shared targets (follow order specified in target_list [esp for bootstrapping, which may have duplicates due to sampling w/ replacement])
    for sample1 in sorted(sampleDict.keys()):
        distDict["sampleComp"][sample1] = {}
        for sample2 in sorted(sampleDict.keys()):
            distDict["sampleComp"][sample1][sample2] = {}
            shared_targets = set(sampleDict[sample1]["HipSTR"].keys()).intersection(set(sampleDict[sample2]["HipSTR"].keys()))
            final_target_list = []
            for target_id in sorted(target_list): #Loop through all targets specified in target_list
                try:
                    if target_id not in shared_targets:
                        continue
                except:
                    pass
                final_target_list.append(target_id)
            #Use calcDist function to calculate genotype_dist between samples and contribution of each target to distance
            distDict = calcDist(sampleDict, alleleDict, distDict, sample1, sample2, final_target_list, dist_metric)
    if verbose is True: #We want to determine useful targets
        info_fields = ["INFRAME_PGEOM","INFRAME_UP","INFRAME_DOWN","OUTFRAME_PGEOM","OUTFRAME_UP","OUTFRAME_DOWN","START","END","PERIOD","NSKIP","NFILT","BPDIFF","DP","DSNP","DSTUTTER","DFLANKINDEL","AN","REFAC","AC"]
        statsOutput = open(prefix + ".stats.out", 'w')
        statsOutput.write("targetID\tIntra-clone Dist\tNum Intra-clone Pairs\tInter-clone Dist\tNum Inter-clone Pairs\tDist Bool\tTotal Dist\tNum Total Pairs\t" + \
            "\t".join(info_fields) + "\t" + "\t".join(sorted(sampleDict.keys())) + "\talleleCount Freq\n")
        for target_id in sorted(distDict["targetComp"].keys()):
            infoDict = {} #Contains overall statistics for each target_id (this will be repeated as same for all samples per target_id)
            msCount = {} #Contains all msCounts for each sample for target_id
            alleleCount = {} #Contains all alleles called for each sample for target_id
            #Calculate average intra (within) and inter (across) clone distances.  Also want to track average pairwise distance across all samples (total_dist)
            (intra_dist, inter_dist, total_dist, num_intra, num_inter, num_total, avg_intra_dist, avg_inter_dist, avg_total_dist) = (0,0,0,0,0,0,0,0,0) #Declare value of -1 for undefined
            for (sample1,sample2) in distDict["targetComp"][target_id].keys():
                infoDict = sampleDict[sample1]["HipSTR"][target_id]["info"]
                if sampleDict[sample1]["clone"] == sampleDict[sample2]["clone"] and sample1 != sample2:
                    intra_dist += distDict["targetComp"][target_id][(sample1,sample2)]
                    num_intra += 1
                elif sampleDict[sample1]["clone"] != sampleDict[sample2]["clone"] and sampleDict[sample1]["clone"] != sampleDict[sample2]["clone"]:
                    inter_dist += distDict["targetComp"][target_id][(sample1,sample2)]
                    num_inter += 1
                total_dist += distDict["targetComp"][target_id][(sample1,sample2)]
                num_total += 1
            avg_total_dist = float(total_dist/num_total) #Average allelotype difference across all samples
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
            #We want to save the msCounts and allotypes for all samples
            HipSTR_call = []
            for sample in sorted(sampleDict.keys()):
                if target_id in sampleDict[sample]["HipSTR"].keys():
                    allelotype_info = '|'.join(str(i) for i in sampleDict[sample]["HipSTR"][target_id]["allelotype"]) + "(" + ','.join(str(j) for j in sampleDict[sample]["HipSTR"][target_id]["stats"]) + ")"
                    HipSTR_call.append(allelotype_info)
                else:
                    HipSTR_call.append('.')
            statsOutput.write(target_id + "\t" + str(round(avg_intra_dist,2)) + "\t" + str(num_intra) + "\t" + str(round(avg_inter_dist,2)) + "\t" + str(num_inter) + "\t" + diff_bool + \
                "\t" + str(round(avg_total_dist,2)) + "\t" + str(num_total) + "\t")
            for info_name in info_fields: #Print overall stutter model info for target_id
                if info_name in infoDict.keys():
                    statsOutput.write(infoDict[info_name] + "\t")
                else:
                    statsOutput.write("NA\t")
            statsOutput.write("\t".join(HipSTR_call) + "\t")
            for allele in sorted(alleleDict[target_id]["countFreq"], key=alleleDict[target_id]["countFreq"].get, reverse=True): #Sort dictionary by value
                statsOutput.write(str(allele) + "," + str(round(float(alleleDict[target_id]["countFreq"][allele]/sum(alleleDict[target_id]["countFreq"].values())), 2)) + ";")
            statsOutput.write("\n")
        statsOutput.close()
    return distDict

def drawTree(distDict, sampleDict, outgroup, prefix, bootstrap):
    # tree_output = open(prefix + ".tree.out", 'w')
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
    if bootstrap is False: #Only output statistics for distance and number targets shared if for original tree (don't output for bootstrap resampling)
        statsOutput = open(prefix + ".stats.out", 'a')
        for dist_indx,dist_list in enumerate(distMatrix): #Print matrix containing distances
            statsOutput.write(sorted(sampleDict.keys())[dist_indx] + "," + ",".join(str(round(i,3)) for i in dist_list) + "\n")
        for target_indx,target_list in enumerate(targetMatrix): #Print matrix containing number targets shared between each pair
            statsOutput.write(sorted(sampleDict.keys())[target_indx] + "," + ",".join(str(j) for j in target_list) + "\n")
        statsOutput.close()
    distObj = DistanceMatrix(distMatrix,sorted(sampleDict.keys()))
    NJTree = nj(distObj, result_constructor=str)
    tree_original = Tree(NJTree) #We use skbio to first make a tree from distance matrix then convert to ete tree
    if outgroup is "NA":
        # tree_output.write(tree_unrooted.write(format=0)) + "\n")
        return tree_original
    else:
        if outgroup == "Midpoint":
            tree_midpoint = tree_original.get_midpoint_outgroup()
            tree_original.set_outgroup(tree_midpoint)
        else:
            tree_original.set_outgroup(outgroup)
    # tree_output.write(tree_original.write(format=0) + "\n")
    return tree_original

def bootstrapTree(nodeDict, treeTemp):
    for temp_node in treeTemp.search_nodes():
        leaf_list = []
        for leaf in temp_node:
            leaf_list.append(leaf.name)
        if tuple(sorted(leaf_list)) in nodeDict.keys():
            nodeDict[tuple(sorted(leaf_list))]["Bootstrap"] += 1
    return nodeDict

def main():
    parser = argparse.ArgumentParser(description="Run and analyze HipSTR output to determine distance between single cells")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specify output file containing pairwise distance calculations")
    parser.add_argument('--targets', action="store", dest="target_file", default="All", help="[Optional] Specify list of targets used for distance calculation")
    parser.add_argument('--dist', action="store", dest="dist_metric", help="Specify distance metric for pairwise comparisons [Abs, EqorNot]", default="EqorNot")
    parser.add_argument('--vcf', action="store", dest="vcf_output")
    parser.add_argument('--min-freq', action="store", dest="min_freq", default=0.1, help="Specify required minimum frequency across all single cell to keep allelotype")
    parser.add_argument('--min-call-qual', action="store", dest="min_qual", default=0, help="Specify the minimum posterior probability of genotype for filtering")
    parser.add_argument('--min-reads', action="store", dest="min_reads", default=1, help="Cutoff for minimum number of reads required for calling allelotype")
    parser.add_argument('--max-stutter', action="store", dest="max_stutter", default=1, help="Define maximum number of reads that can be classified as stutter")
    parser.add_argument('--outgroup', action="store", dest="outgroup", default="NA", help="[Optional] Specify outgroup for rooted NJ tree (if use midpoint, specify 'Midpoint')")
    parser.add_argument('-v', action="store_true", help="Flag for determining whether we want to output all statistics for shared targets in output")
    parser.add_argument('-plot', action="store_true", help="Flag for indicating whether we want to output a plot file visualizing msCounts per targetID")
    parser.add_argument('--num_plot', action="store", dest="num_plot", default=1784, help="If plotting allelotype, specify number of targets to plot")
    parser.add_argument('-bootstrap', action="store_true", help="Flag for indicating whether we want to bootstrap the tree to determine node support (random sampling shared targets per cell pair)")
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
                (bam, sample) = line.split()
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

    #Calculate "bulk" allelotype counts by merging all single cell HipSTR calls
    alleleDict = calcBulk(sampleDict, float(args.min_freq))

    #Plot msCounts and allelotype if -plot flag is used
    if args.plot is True:
        plotVCF(sampleDict, args.prefix, int(args.num_plot))

    #Calculate original tree using all targets found within alleleDict (or subset as specified in args.target_file)
    if args.target_file is "All":
        target_list = sorted(alleleDict.keys())
    else:
        with open(target_file) as f:
            target_list = f.read().splitlines()
    #Calculate pairwise distance
    distDict_original = makeDistMatrix(sampleDict, alleleDict, target_list, args.dist_metric, args.v, args.prefix, False)
    #Draw neighbor joining tree
    tree_original = drawTree(distDict_original, sampleDict, args.outgroup, args.prefix, False) #We want to declare Fase for args.bootstrap because we want to output stats file for original tree

    #Bootstrap resample to create new distance matrices/trees and add support values to internal nodes of original tree
    if args.bootstrap is True:
        #Determine dictionary of tree nodes from tree_original
        nodeDict = {}
        for node in tree_original.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name) #We need to compare leaves in a node cluster without regard to order
            nodeDict[tuple(sorted(leaf_list))] = {}
            nodeDict[tuple(sorted(leaf_list))]["Bootstrap"] = 0 #Contains values for number of times random bootstrap tree (tree_temp) contains given node
            nodeDict[tuple(sorted(leaf_list))]["NodeID"] = node.write(format = 9)
        for i in tqdm(range(10)): #Bootstrap resample 1000 times
            #Random downsample from pool of available targets (target_list) to use for distance calculation
            bootstrap_targets = list(np.random.choice(target_list, len(target_list), replace=True))
            distDict_temp = makeDistMatrix(sampleDict, alleleDict, bootstrap_targets, args.dist_metric, False, args.prefix, args.bootstrap)
            tree_temp = drawTree(distDict_temp, sampleDict, args.outgroup, args.prefix, args.bootstrap)
            nodeDict = bootstrapTree(nodeDict, tree_temp) #Determine whether each node in original tree is found in tree_temp
        #Add support information to original tree
        for node in tree_original.search_nodes():
            leaf_list = []
            for leaf in node:
                leaf_list.append(leaf.name)
            node.add_features(support = round(nodeDict[tuple(sorted(leaf_list))]["Bootstrap"]/10,2))

    #Output tree (if bootstrapped, output with support values)
    tree_output = open(args.prefix + ".tree.out", 'w')
    tree_output.write(tree_original.write(format = 0) + "\n")
    tree_output.close()

if __name__ == "__main__":
    main()
