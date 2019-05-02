import numpy as np
# import itertools
# from joblib import Parallel, delayed
# import multiprocessing
# from random import sample
# from tqdm import tqdm

def simple_counter(read, targetDict, target_id):
    '''This is a simple counter to determine microsatellite subunits with exact match of up_seq and down_seq'''
    num_bases = -1
    if targetDict[target_id]["up_seq"] in read and targetDict[target_id]["down_seq"] in read:
        up_start = read.index(targetDict[target_id]["up_seq"])
        down_start = read.index(targetDict[target_id]["down_seq"])
        num_bases = down_start-(up_start+len(targetDict[target_id]["up_seq"]))
    return num_bases


def construct(read, flank, bool): #bool=0 for up_seq, bool=1 for down_seq, bool=2 for pseudo_ref (ms).  We want to make it so that the score for match decreases as you get towards the microsatellite region for both up/down-seq.  That way, the local alignment program would more likely choose to skip rather than psuh for a match in those locations.
    read_length = len(read)
    flank_length = len(flank)
    read_list = list(read)
    flank_list = list(flank)

    #We want to start backtrack from max s[i][j].  Consequently, we want to determine the maximum value of score and the location within the read (max_i) and flanking sequence (max_j).  This will then be returned from the construct subroutine
    max_s = float("-inf") #With i being the read position and j being the flanking region
    #Initialize scoring and backtrack matrices.  We want to avoid using numpy arrays because it would cause considerable slowdown
    # s = [0] * (read_length + 1)
    # backtrack = ['NA'] * (read_length + 1)
    # for i in range(read_length + 1):
    #     s[i] = [0] * (flank_length + 1)
    #     backtrack[i] = ['NA'] * (flank_length + 1)
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
    # print(str(len(read)) + "\t" + read)
    '''We incoporate "fuzzy" zones in the alignment where the regions at the ends of the suspected MS site and the regions of the up/down-seqs nearest the MS site (nearest the middle) will expreince less of an advantage to try to find a match.  This will push the program to just skip these regions in the local alignment, which in the end makes the MS call more accurate by getting rid of artifacts that occur more often in these more non-unique regions of the read'''
    #1) We firs twant to use local alignment to substract out the up-stream sequence from teh read by determining the last base of the highest-socring match of the up-seq within the read
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

def counter(count_type, read_list, targetDict, target_id, prefix):
    #Run microsatellite counting in parallel
    print(prefix + "\tAnalyzing:\t" + target_id)
    msCount_list = []
    for read in read_list:
        if count_type == "aln":
            msCount_list.append(read + "=" + str(aln_counter(read, targetDict, target_id))) #Save both read sequence and msCount
        elif count_type == "simple":
            msCount_list.append(read + "=" + str(simple_counter(read, targetDict, target_id))) #Save both read sequence and msCount
    return msCount_list
