#!/usr/bin/env python3
import numpy as np
import multiprocessing
import pysam
from tqdm import tqdm
import random
manager = multiprocessing.Manager()

def construct(read, flank, bool): #bool=0 for up_seq, bool=1 for down_seq, bool=2 for pseudo_ref (ms).  We want to make it so that the score for match decreases as you get towards the microsatellite region for both up/down-seq.  That way, the local alignment program would more likely choose to skip rather than psuh for a match in those locations.
    read_length = len(read)
    flank_length = len(flank)
    read_list = list(read)
    flank_list = list(flank)

    #We want to start backtrack from max s[i][j].  Consequently, we want to determine the maximum value of score and the location within the read (max_i) and flanking sequence (max_j).  This will then be returned from the construct subroutine
    max_s = float("-inf") #With i being the read position and j being the flanking region
    s = np.empty((read_length + 1, flank_length + 1), dtype=int)
    backtrack = np.empty((read_length + 1, flank_length + 1), dtype=object)
    for i in range(read_length + 1):
        s[i][0] = 0
        backtrack[i][0] = 'down'
    for j in range(flank_length + 1):
        s[0][j] = 0
        backtrack[0][j] = 'right'
    for i in range(1, read_length + 1):
        for j in range(1, flank_length + 1):
            if read_list[i-1] == flank_list[j-1]: #If nucleotides match
                comp1 = float("-inf")
                if bool == 0: #for up_seq local alignment
                    if j <= 0.8*flank_length: #We want to keep the advantage of a match to be constant +1 up until you are 0.8 of the length of up_length.  Once you are past that, you want to taper off drastically and decrease the advantage of a match
                        comp2 = s[i-1][j-1] + 1
                    else:
                        comp2 = s[i-1][j-1] + 1/j
                elif bool == 1: #for down-seq local alignment
                    if j > 0.2*flank_length:
                        comp2 = s[i-1][j-1] + 1
                    else:
                        comp2 = s[i-1][j-1] + 1/(1 + (flank_length - j))
                else:
                    comp2 = s[i-1][j-1] + 1
            else: #If there is a mismatch
                comp1 = s[i-1][j-1] - 1
                comp2 = float("-inf")
            #Penalty for indel
            comp3 = s[i-1][j] - 1 #If there is an insertion inf the read compared to the flank_seq
            comp4 = s[i][j-1] - 1 #If ther eis a deletion int eh read compared ot the flank_seq
            comp5 = 0 #If there is a skip in the local alignment
            s[i][j] = max(comp1, comp2, comp3, comp4, comp5)
            #We want to define backtrack
            if s[i][j] == comp1 or s[i][j] == comp2:
                backtrack[i][j] = 'diag'
            elif s[i][j] == comp3:
                backtrack[i][j] = 'down'
            elif s[i][j] == comp4:
                backtrack[i][j] = 'right'
            elif s[i][j] == comp5:
                backtrack[i][j] = 'skip'
            #We now want to determine and adjust max_s, max_i, and max_j depending ont he value of s[i][j]
            if s[i][j] > max_s:
                max_s = s[i][j]
                max_i = i
                max_j = j
    return backtrack, max_i, max_j, max_s


def outputLCS(backtrack, i, j): #Adjustment variable may or may not be used depending on whether it is defined or not originally in the initial call of outputLCS (whether it is a saved in the output of the initial call of outputLCS)
    while i > 0 and j > 0:
        if backtrack[i][j] == 'down': #If there is an insertion in the read compared to the flank_seq
            i += -1
        elif backtrack[i][j] == 'right': #If there is a deletion i nthe read compared to flank_seq
            j += -1
        elif backtrack[i][j] == 'diag':
            i += -1
            j += -1
        elif backtrack[i][j] == 'skip':
            return (i, j)
    return (i, j)


def aln_counter(read, targetDict, target_id):
    '''We incoporate "fuzzy" zones in the alignment where the regions at the ends of the suspected MS site and the regions of the up/down-seqs nearest the MS site (nearest the middle) will expreince less of an advantage to try to find a match.  This will push the program to just skip these regions in the local alignment, which in the end makes the MS call more accurate by getting rid of artifacts that occur more often in these more non-unique regions of the read'''
    #1) We first want to use local alignment to substract out the up-stream sequence from teh read by determining the last base of the highest-socring match of the up-seq within the read
    (backtrack_up, max_i_up, max_j_up, max_s_up) = construct(read, targetDict[target_id]["up_seq"], 0) #max_i signifies position within read containing greatest score, max_j position within flank containing greatest score
    if max_s_up <= 0.5*len(targetDict[target_id]["up_seq"]): #Remove reads with low alignment score
        return -1
    ms_start  = max_i_up + (len(targetDict[target_id]["up_seq"]) - max_j_up) #This will show the true position of the microsatellite start by taking into account the flanking region that was not mapped

    #2) We want to use local alignment to subtract out the down-stream sequence form the read by determining the first base of the highest-socring match of the down-seq within the read
    (backtrack_down, max_i_down, max_j_down, max_s_down) = construct(read, targetDict[target_id]["down_seq"], 1)
    if max_s_down <= 0.5 * len(targetDict[target_id]["down_seq"]):
        return -1
    (min_i_down, min_j_down) = outputLCS(backtrack_down, max_i_down, max_j_down)
    ms_end = min_i_down - min_j_down #This will show the true position of the micorsatellite end by taking into account the flanking regions that was not mapped

    #3) We finally want to determine the number of subunits by subtracking ms_end and ms_start and dividding tby the length of the microsatellite subunit
    num_bases = ms_end - ms_start #We base all of our calculations off of the absolute number of bases rather than the number of subunits

    return num_bases

def simple_counter(read, targetDict, target_id):
    '''This is a simple counter to determine microsatellite subunits with exact match of up_seq (last 10 bases) and down_seq'''
    up_start = read.index(targetDict[target_id]["up_seq"][-10:])
    down_start = read.index(targetDict[target_id]["down_seq"][:10])
    num_bases = down_start-(up_start+len(targetDict[target_id]["up_seq"][-10:]))

    return num_bases

def counter(read_list, targetDict, target_id):
    #Run microsatellite counting in parallel
    print("Analyzing:\t" + target_id)
    msCount_list = []
    for read in read_list:
        if targetDict[target_id]["up_seq"][-10:] in read and targetDict[target_id]["down_seq"][:10] in read: #Speed up processing by simple alignment if exact match to portion of up/down-seq
            msCount_list.append(read + "=" + str(simple_counter(read, targetDict, target_id)))
        else:
            msCount_list.append(read + "=" + str(aln_counter(read, targetDict, target_id))) #Save both read sequence and msCount
    return msCount_list

def multi_counter(target_group, targetDict, sampleDict, counterDict):
    '''
    This module allos for multiprocessing of the counter by allowing target_id's to be processed in parallel
    '''
    for target_id in tqdm(sorted(target_group)):
        read_list = [] #Keep list of unique read sequences found across all samples
        for sample in sorted(sampleDict.keys()): #Iterate through all sample bams in order to extract all of the reads for given target
            samfile = pysam.AlignmentFile(sampleDict[sample]["bam"], "rb")
            for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
                if read.query not in read_list:
                    read_list.append(read.query)
        counterDict[target_id] = counter(read_list, targetDict, target_id)
    return

def msCount(sampleDict, prefix, targetDict, nproc):
    '''
    This is a wrapper script for running multiple instances of msCount based on nproc.  In order to do this, we can split target_id's into chunks
    This is based off of: <https://further-reading.net/2017/01/quick-tutorial-python-multiprocessing/>, <https://stackoverflow.com/questions/10415028/how-can-i-recover-the-return-value-of-a-function-passed-to-multiprocessing-proce>
    Save all of the msCount in alleleDict with the following structure:
        alleleDict
            target_id (from targetDict)
                "sample"
                    sample (from sampleDict, which is already defined when labeling readGroups prior to HipSTR)
                        "msCount"
                            list of msCounts
                        "allelotype"
                            list of alleles (2 alleles)
    '''
    print("Performing msCount from raw reads")
    counterDict = manager.dict() #This allows for parallel processing of msCount for multiple target_id at once
    jobs = []
    target_list = list(targetDict.keys())
    random.shuffle(target_list) #We want to randomize in order to even processing time
    for target_group in [target_list[i::nproc] for i in range(nproc)]:
        p = multiprocessing.Process(target = multi_counter, args = (target_group, targetDict, sampleDict, counterDict))
        jobs.append(p)
        p.start()

    #Join counterDict
    for p in jobs:
        p.join()

    #Re-iterate through samples in order to assign msCount to teh reads aligned for given sample (save in alleleDict)
    msCount_output = open(prefix + ".Custom.msCount.csv", 'w')
    alleleDict = {}
    for target_id in sorted(counterDict.keys()):
        alleleDict[target_id] = {}
        #We want to extract all read and msCount info from counterDict
        read_list = [str(count_info.split('=')[0]) for count_info in counterDict[target_id]]
        msCount_list = [int(count_info.split('=')[1]) for count_info in counterDict[target_id]]
        #Assign appropriate msCounts to each sample and save in alleleDict
        alleleDict[target_id]["sample"] = {}
        for sample in sorted(sampleDict.keys()):
            samfile = pysam.AlignmentFile(sampleDict[sample]["bam"], "rb")
            for read in samfile.fetch(targetDict[target_id]["chrom"], int(targetDict[target_id]["chromStart"]), int(targetDict[target_id]["chromEnd"])):
                try:
                    count_indx = read_list.index(read.query)
                except:
                    print("Error: Unable to find read within sample")
                    return
                if msCount_list[count_indx] > 0: #We only want to save valid msCounts for each sample
                    if sample not in alleleDict[target_id]["sample"].keys():
                        alleleDict[target_id]["sample"][sample] = {}
                        alleleDict[target_id]["sample"][sample]["msCount"] = []
                    alleleDict[target_id]["sample"][sample]["msCount"].append(msCount_list[count_indx])
            samfile.close()
            #We want to print the msCount info for each target_id and sample if exists
            if sample in alleleDict[target_id]["sample"].keys():
                msCount_output.write(target_id + ',' + sample + ',' + ','.join(str(msCount) for msCount in alleleDict[target_id]["sample"][sample]["msCount"]) + "\n")
    return alleleDict
