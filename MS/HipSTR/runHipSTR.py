#!/usr/bin/env python3
import argparse
import os
import subprocess

'''Usage: python script.py --sample_list input.txt
This script will take as input a sample list file that contains the following information:
    1) File location of sample bam
    2) Sample name
    3) [Optional] Sample type (i.e. known clone in ex vivo tree)
It will then perform read group assignment and run HipSTR on the given files.  We will then parse the HipSTR vcf file output and extract the allelotype/msCounts from each given sample.
The HipSTR vcf file have the followign format for each microsatellite target captured:
    #CHROM  POS_ID  REF ALT QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 ...
'''

def run_HipSTR(sampleDict, software_dir):
    for sample in sampleDict.keys():
        bam_readGroup = sampleDict[sample]["bam"].split('.')[:-1] + ".readGroup.bam"
        sampleDict["bam_readGroup"] = bam_readGroup
        if not os.path.isfile('./' + bam_bname + ".readGroup.bam"):
            readGroup_command = "java -jar " + software_dir + "picard.jar AddOrReplaceReadGroups " +
                "I=" + sampleDict[sample]["bam"] + " " + "O=" + bam_readGroup + " " +
                "RGLB=" + sample + " " + "RGPL=illumina RGPU=unit1 RGSM" + sample
            index_command = "samtools index " + bam_readGroup
            subprocess.call(readGroup_command, shell=True)
            subprocess.call(index_command, shell=True)
    hipSTR_command = software_dir + "HipSTR/HipSTR --bams " + ",".join(sorted(sampleDict["bam_readGroup"])) +
        "--fasta /media/6TB_slot4/GenomeDB/hg19/raw_fasta/hg19_reference.fa " +
        "--regions /home/cjwei/software/RETrace/MS/HipSTR/20171211_probes.bed " +
        "--str-vcf " +

def main():
    parser = argparse.ArgumentParser(description="Run and analyze HipSTR output to determine distance between single cells")
    parser.add_argument('--software_dir', action="store", dest="software_dir", default='~/software/', help="Sofware directory containing picard.jar and HipSTR install")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--vcf', action="store", dest="vcf_output")

    #We want to first parse the sample_info file
    sampleDict = {}
    with open(args.sample_info) as f:
        for line in f:
            if len(line.split) == 3: #If clone is specified
                (bam, sample, clone) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
                sampleDict[sample]["clone"] = clone
            elif len(line.split) == 2: #If clone is not specified
                (bam, sampleName) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
            else:
                print("Incorrect formmating for line (bam, sampleName, [optional] clone):\n" + "\t".join(line.split))
                return

    #Run HipSTR if you haven't done so already
    if not os.path.isfile('./' + args.vcf_output):
        run_HipSTR(sampleDict, args.software_dir, args.vcf_output)

if __name__ == "__main__":
    main()
