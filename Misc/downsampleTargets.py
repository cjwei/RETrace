#!/usr/bin/env python3
import argparse
import pickle
import random

def import_targetDict(probe_file):
    '''Import probe file into targetDict'''
    targetDict = {}
    with open(probe_file) as f:
        for line in f:
            target_info = line.split(',')[0]
            info_list = target_info.split('_')
            #Reformat target_id and save to targetDict
            chrom = "_".join(info_list[0:len(info_list)-3])
            chromStart = info_list[-3]
            chromEnd = info_list[-2]
            target_id = chrom + ":" + chromStart + "-" + chromEnd

            #Save relevant information to targetDict
            targetDict[target_id] = {}
            (targetDict[target_id]["chrom"], targetDict[target_id]["chromStart"], targetDict[target_id]["chromEnd"]) = (chrom, chromStart, chromEnd)
            (targetDict[target_id]["num_sub"], targetDict[target_id]["sub_seq"]) = info_list[-1].split('x')
            #Also save the expected up/down_seq around the microsatellite
            MS_frag = line.split()[1]
            targetDict[target_id]["MS_frag"] = MS_frag
            frag_seq = targetDict[target_id]["sub_seq"]*int(targetDict[target_id]["num_sub"])
            frag_split = MS_frag.split(frag_seq)
            #There are a few cases (i.e. chr6:55179847-55179875) in which microsatellite frag_seq is found more than once in reference fragment.  Thus, we need to choose the optimal up/down-seq based on which pair has the greater length
            (up_seq, down_seq) = ('', '')
            for i in range(len(frag_split) - 1):
                if len(frag_split[i]) > len(up_seq) and len(frag_split[i + 1]) > len(down_seq):
                    (up_seq, down_seq) = frag_split[i:i+2]
            (targetDict[target_id]["up_seq"], targetDict[target_id]["down_seq"]) = (up_seq, down_seq)
            # targetDict[target_id]["sample_msCount"] = {} #Place holder for sample msCounts in targetDict
    return targetDict

def downsampleTarget(alleleDict_file, target_info, numTargets, output_alleleDict):
    '''
    This script will downsample the alleleDict file such that samples have the user-specified numTargets.
    If samples have < numTargets, remove from alleleDict["samples"]
    '''

    alleleDict = pickle.load(open(alleleDict_file, 'rb'))
    targetDict = import_targetDict(target_info)

    #We want to run through the alleleDict to obtain initial list of target_ids per sample
    sampleTargets = {} #Keep track of the targets that were captured per sample
    for target_id in sorted(alleleDict.keys()):
        if target_id in targetDict.keys(): #We only want to analyze target_id specified in targetDict
            for sample in sorted(alleleDict[target_id]["sample"].keys()):
                if sample not in sampleTargets.keys():
                    sampleTargets[sample] = {}
                    sampleTargets[sample]["original"] = set()
                sampleTargets[sample]["original"].add(target_id)

    #Downsample target_ids in sampleTargets
    for sample in sampleTargets.keys():
        if len(sampleTargets[sample]["original"]) >= numTargets:
            sampleTargets[sample]["downsample"] = random.sample(sampleTargets[sample]["original"], numTargets)
        else:
            sampleTargets[sample]["downsample"] = set() #Create empty set if sample does not satisfy numTargets requirement

    #Run through alleleDict again and remove samples from calls
    for target_id in sorted(alleleDict.keys()):
        original_numSamples = len(alleleDict[target_id]["sample"].keys())
        for sample in sorted(alleleDict[target_id]["sample"].keys()):
            if target_id not in sampleTargets[sample]["downsample"]:
                del alleleDict[target_id]["sample"][sample]

    #Save alleleDict into output file
    pickle.dump(alleleDict, open(output_alleleDict, 'wb'))

def main():
    parser = argparse.ArgumentParser(description="Randomly downsample targets from alleleDict")
    parser.add_argument('--alleleDict', action="store", dest="alleleDict_file", help="Original alleleDict file")
    parser.add_argument('--target_info', action="store", dest="target_info", help="Target_info file location")
    parser.add_argument('--numTargets', action="store", dest="numTargets", type=int, help="Number of targets for random downsampling")
    parser.add_argument('--output', action="store", dest="output_alleleDict", help="Output alleleDict file containing only randomly downsampled target_id")
    args = parser.parse_args()

    downsampleTarget(args.alleleDict_file, args.target_info, args.numTargets, args.output_alleleDict)

if __name__ == "__main__":
    main()
