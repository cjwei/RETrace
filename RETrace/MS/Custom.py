#!/usr/bin/env python3
from RETrace.MS.utilities import import_sampleDict, import_targetDict
from RETrace.MS.Custom_msCount import msCount
from RETrace.MS.Custom_allelotype import allelotype
import os
import pickle

def import_msCount(prefix):
    '''
    We want to import alleleDict if msCount has already been performed previously.
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
    alleleDict = {}
    with open(prefix + '.Custom.msCount.csv', 'r') as msCount_output:
        for line in msCount_output:
            msCount_line = line.split(',')
            target_id = msCount_line[0]
            sample = msCount_line[1]
            if target_id not in alleleDict.keys():
                alleleDict[target_id] = {}
                alleleDict[target_id]["sample"] = {}
            if sample not in alleleDict[target_id]["sample"].keys():
                alleleDict[target_id]["sample"][sample] = {}
            alleleDict[target_id]["sample"][sample]["msCount"] = [int(msCount) for msCount in msCount_line[2:]]
    return alleleDict


def Custom_allelotype(sample_info, prefix, target_info, nproc, min_cov, min_ratio):
    '''
    This funciton will perform custom microsatellite counting with "fuzzy" ends and calculate allelotype based off of a stutter model similar to LobSTR.  The following are required inputs:
        sample_info = tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)
        target_info = location of probe info file
    '''

    sampleDict = import_sampleDict(sample_info)

    targetDict = import_targetDict(target_info)

    if not os.path.isfile(prefix + '.Custom.alleleDict.pkl'):
        if not os.path.isfile(prefix + '.Custom.msCount.csv'):
            alleleDict = msCount(sampleDict, prefix, targetDict, nproc) #We want to perform msCount if prefix.Custom_msCount.csv is not found within directory
        else:
            alleleDict = import_msCount(prefix) #Import msCount information for
        #Once we have msCount imported into alleleDict, we next want to calculate allelotype using a LobSTR-based stutter model
        alleleDict = allelotype(sampleDict, prefix, targetDict, alleleDict, min_cov, min_ratio)
        pickle.dump(alleleDict, open(prefix + ".Custom.alleleDict.pkl", "wb"))
    else:
        print("Allelotype already available at:\t" + prefix + ".Custom.alleleDict.pkl")
