#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict, import_targetDict
import os
import subprocess
import pickle

def run_HipSTR(sampleDict, HipSTR_vcf, fasta_loc, picard_loc, HipSTR_loc, target_bed):
    '''
    Wrapper command to add readgroup to bam files and run HipSTR
    '''
    print("Running HipSTR")
    readGroup_list = []
    for sample in sorted(sampleDict.keys()):
        bam_readGroup = '.'.join(sampleDict[sample]["bam"].split('.')[:-1]) + ".readGroup.bam"
        readGroup_list.append(bam_readGroup)
        sampleDict[sample]["bam_readGroup"] = bam_readGroup
        if not os.path.isfile(bam_readGroup):
            readGroup_command = "java -jar " + picard_loc + "/picard.jar AddOrReplaceReadGroups " + \
                "I=" + sampleDict[sample]["bam"] + " " + "O=" + bam_readGroup + " " + \
                "RGLB=" + sample + " " + "RGPL=illumina RGPU=unit1 RGSM=" + sample
            index_command = "samtools index " + bam_readGroup
            subprocess.call(readGroup_command, shell=True)
            subprocess.call(index_command, shell=True)
    HipSTR_command = HipSTR_loc + "/HipSTR --bams " + ",".join(sorted(readGroup_list)) + " " + \
        "--fasta " + fasta_loc + " " + \
        "--regions " + target_bed + " " + \
        "--str-vcf " + HipSTR_vcf + ".gz --log " + HipSTR_vcf.replace('.vcf.gz', '') + ".log --use-unpaired --no-rmdup"
    print(HipSTR_command)
    gunzip_command = "gunzip " + HipSTR_vcf + ".gz"
    subprocess.call(HipSTR_command, shell=True)
    subprocess.call(gunzip_command, shell=True)
    return

def parseVCF(sampleDict, targetDict, HipSTR_vcf, alleleDict_file, alleleOutput_file, min_qual, min_reads, max_stutter):
    '''
    Parse HipSTR vcf output into alleleDict with the following structure:
        alleleDict
            target_id (from targetDict)
                "sample"
                    sample (from sampleDict, which is already defined when labeling readGroups prior to HipSTR)
                        "msCount"
                            list of msCounts
                        "allelotype"
                            list of alleles (2 alleles)
    '''
    print("Extracting allelotype from HipSTR output")
    alleleDict = {}
    with open(HipSTR_vcf) as vcf:
        for line in vcf:
            if not line.startswith('##'):
                if line.startswith('#'): #This contains the vcf fields info (comes before any allelotyping calls within vcf)
                    info_line = line.split()
                    sample_names = info_line[9:] #Contains sample names from HipSTR analysis
                else:
                    vcf_line = line.split()
                    target_id = '_'.join(vcf_line[2].split('_')[0:-1]) #This contains target_id that matches with targetDict (eg. "chr1:14290254-14290298"); also takes into account some chromosomes with '_' in the middle (such as "chrUn_*")
                    if target_id not in targetDict.keys():
                        print(target_id + " not found in targetDict.  Aborting...")
                        return
                    else:
                        alleleDict[target_id] = {}
                        alleleDict[target_id]["sample"] = {}
                    sample_format = vcf_line[9:]
                    for indx, val in enumerate(sample_format):
                        sample = sample_names[indx]
                        if val is not '.':
                            allelotype = val.split(':')[1] #GB (base pair differences of genotype from reference)
                            prob_genotype = val.split(':')[2] #Q (posteriior probability of unphased genotype)
                            num_reads = val.split(':')[4] #DP (number of reads used for sample's genotype)
                            num_stutter = val.split(':')[6] #DSTUTTER (number of reads with a stutter indel in the STR region)
                            msCounts = val.split(':')[-2] #ALLREADS (bp difference observed in each read's alignment)
                            if float(prob_genotype) >= min_qual and int(num_reads) >= min_reads and float(int(num_stutter)/int(num_reads)) <= max_stutter and msCounts is not '.':
                                ref_MS_len = len(targetDict[target_id]["sub_seq"]) * int(targetDict[target_id]["num_sub"])
                                alleleDict[target_id]["sample"][sample] = {}
                                alleleDict[target_id]["sample"][sample]["msCount"] = []
                                for msCount_info in msCounts.split(';'):
                                    [msCount, freq] = msCount_info.split('|')
                                    alleleDict[target_id]["sample"][sample]["msCount"].extend([int(msCount)] * int(freq))
                                alleleDict[target_id]["sample"][sample]["allelotype"] = [int(allele) + ref_MS_len for allele in allelotype.split('|')]
    #We want to output the allelotype calls into a tabulated file for publication
    if alleleOutput_file:
        f_output = open(alleleOutput_file, 'w')
        f_output.write("targetID\t" + "\t".join(sorted(sample_names)) + "\n")
        for target_id in sorted(alleleDict.keys()):
            f_output.write(target_id + "\t")
            for sample in sorted(sample_names):
                if sample in alleleDict[target_id]["sample"].keys():
                    f_output.write('|'.join(str(allele) for allele in alleleDict[target_id]["sample"][sample]["allelotype"]) + "\t")
                else:
                    f_output.write('.|.' + "\t")
            f_output.write("\n")
        f_output.close()
    return alleleDict

def HipSTR_allelotype(sample_info, HipSTR_vcf, alleleDict_file, alleleOutput_file,
    fasta_loc, picard_loc, HipSTR_loc,
    target_bed, target_info,
    min_qual, min_reads, max_stutter):
    '''
    This function serves as a wrapper for running HipSTR and importing the allelotype information.  The following are required inputs:
        sample_info = tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)
        picard_loc = directory location containing picard command
        HipSTR_loc = directory location containing HipSTR command
        target_bed = location of probe bed file
        target_info = location of probe info file
    We also want to output the allelotype calls into a tab-delimited file for publication
    '''
    sampleDict = import_sampleDict(sample_info)

    targetDict = import_targetDict(target_info)

    #Run HipSTR if you haven't done so already
    if not os.path.isfile('./' + HipSTR_vcf):
        run_HipSTR(sampleDict, HipSTR_vcf, fasta_loc, picard_loc, HipSTR_loc, target_bed)

    #Import HipSTR to alleleDict
    if not os.path.isfile('./' + alleleDict_file):
        alleleDict = parseVCF(sampleDict, targetDict, HipSTR_vcf, alleleDict_file, alleleOutput_file, min_qual, min_reads, max_stutter)
        pickle.dump(alleleDict, open(alleleDict_file, "wb"))
    else:
        print("Allelotype already available at:\t" + alleleDict_file)
