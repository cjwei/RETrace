#!/usr/bin/env python3
import pickle
from tqdm import tqdm

def merge_allelotype(file_list, output):
    '''
    Merge dictionaries containing allelotypes in order to combine multiple datasets together.
    The structure of alleleDict is as follows:
        alleleDict
            target_id (from targetDict)
                "sample"
                    sample (from sampleDict, which is already defined when labeling readGroups prior to HipSTR)
                        "msCount"
                            list of msCounts
                        "allelotype"
                            list of alleles (2 alleles)
    The purpose of this script is that we need to allelotype for each experiment separately because of variations in stutter rate especially during capture.
    We then can merge alleleDicts in order to build phylogeny using samples across multiple experiments
    '''

    alleleDict = {}
    for f in file_list:
        alleleDict_temp = pickle.load(open(f, 'rb'))
        print("Merging:\t" + f)
        for target_id in tqdm(sorted(alleleDict_temp.keys())):
            for sample_name in alleleDict_temp[target_id]["sample"]:
                if "allelotype" in alleleDict_temp[target_id]["sample"][sample_name].keys():
                    if target_id not in alleleDict.keys():
                        alleleDict[target_id] = {}
                        alleleDict[target_id]["sample"] = {}
                    if sample_name not in alleleDict[target_id]["sample"]:
                        alleleDict[target_id]["sample"][sample_name] = {}
                    alleleDict[target_id]["sample"][sample_name]["msCount"] = alleleDict_temp[target_id]["sample"][sample_name]["msCount"]
                    alleleDict[target_id]["sample"][sample_name]["allelotype"] = alleleDict_temp[target_id]["sample"][sample_name]["allelotype"]
    pickle.dump(alleleDict, open(output, "wb"))
