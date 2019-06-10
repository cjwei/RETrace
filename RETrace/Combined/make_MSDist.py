#!/usr/bin/env python3
import more_itertools
import pandas as pd
from tqdm import tqdm

def calcBulk(alleleDict, targetDict):
    '''
    Cluster all alleles belonging to given target across all sample.  Use these as possible "allele_groups"
    '''
    for target_id in sorted(alleleDict.keys()):
        if target_id in targetDict.keys(): #We only want to analyze target_id specified in targetDict
            sub_len = len(targetDict[target_id]["sub_seq"])
            all_alleles = set()
            for sample in sorted(alleleDict[target_id]["sample"].keys()):
                try:
                    all_alleles.update([int(allele/sub_len) for allele in alleleDict[target_id]["sample"][sample]["allelotype"]]) #Convert alleles to number of subunits (because we want to group all consecutive integer number subunit alleles)
                except:
                    continue
            alleleDict[target_id]["allele_groups"] = {}
            for indx, group in enumerate(more_itertools.consecutive_groups(sorted(set(all_alleles)))): #Need to use "set(filtered_alleles)" in order to prevent repeated int in filtered_alleles
                alleleDict[target_id]["allele_groups"][indx] = [allele * sub_len for allele in list(group)] #Save allele groups in terms of raw number of bases difference from ref
    return alleleDict

def calcDist_MS(sample1, sample2, alleleDict, sharedDict, dist_metric):
    '''
    Calculate distance between two samples given their shared_targets and dist_metric
    '''
    total_dist = 0
    num_alleles = 0
    sample_pair = tuple(sorted([sample1, sample2]))
    for target_id in sharedDict[sample_pair]:
        target_dist = 0
        target_alleles = 0
        for allele_indx, allele_group in alleleDict[target_id]["allele_groups"].items():
            allelotype1 = list(set(alleleDict[target_id]["sample"][sample1]["allelotype"]).intersection(set(allele_group)))
            allelotype2 = list(set(alleleDict[target_id]["sample"][sample2]["allelotype"]).intersection(set(allele_group)))
            if len(allelotype1) == 2 and len(allelotype2) == 2:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]) + abs(allelotype1[1] - allelotype2[1]), abs(allelotype1[1] - allelotype2[0]) + abs(allelotype1[0]- allelotype2[1]))
                elif dist_metric == "EqorNot":
                    target_dist += int(len(set(allelotype1).symmetric_difference(set(allelotype2)))/2)
                target_alleles = 2
            elif len(allelotype1) == 1 and len(allelotype2) == 2:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]), abs(allelotype1[0] - allelotype2[1]))
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0] and allelotype1[0] != allelotype2[1]:
                        target_dist += 1
                target_alleles = 1
            elif len(allelotype1) == 2 and len(allelotype2) == 1:
                if dist_metric == "Abs":
                    target_dist += min(abs(allelotype1[0] - allelotype2[0]), abs(allelotype1[1] - allelotype2[0]))
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0] and allelotype1[1] != allelotype2[0]:
                        target_dist += 1
                target_alleles = 1
            elif len(allelotype1) == 1 and len(allelotype2) == 1:
                if dist_metric == "Abs":
                    target_dist += abs(allelotype1[0] - allelotype2[0])
                elif dist_metric == "EqorNot":
                    if allelotype1[0] != allelotype2[0]:
                        target_dist += 1
                target_alleles = 1
        num_alleles += target_alleles
        total_dist += target_dist
    if num_alleles > 0:
        MS_dist = float(total_dist / num_alleles)
    else:
        MS_dist = "NA"
    return MS_dist

def make_MSDist(filtered_samples, targetDict, alleleDict, dist_metric, prefix):
    '''
    Calculate MS distance between two samples and return distance matrix
    '''

    #Calculate "bulk" allelotype counts by merging all single cell calls
    alleleDict = calcBulk(alleleDict, targetDict)

    #We want to pre-compute shared targets between each sample
    sharedDict = {} #Contains all shared targets between each pairwise sample
    for sample1 in sorted(filtered_samples):
        for sample2 in sorted(filtered_samples):
            sample_pair = tuple(sorted([sample1, sample2]))
            if sample_pair not in sharedDict.keys():
                shared_targets = []
                for target_id in sorted(alleleDict.keys()):
                    if sample1 in alleleDict[target_id]["sample"].keys() and sample2 in alleleDict[target_id]["sample"].keys() and target_id in targetDict.keys(): #Only want to analyze target_id in targetDict
                        if "allelotype" in alleleDict[target_id]["sample"][sample1].keys() and "allelotype" in alleleDict[target_id]["sample"][sample2].keys():
                            shared_targets.append(target_id)
                sharedDict[sample_pair] = shared_targets

    #Calculate pairwise distance for all samples depending on shared targets
    distDict = {}
    for sample1 in tqdm(sorted(filtered_samples)):
        distDict[sample1] = []
        for sample2 in sorted(filtered_samples):
            MS_dist = calcDist_MS(sample1, sample2, alleleDict, sharedDict, dist_metric)
            distDict[sample1].append(MS_dist)

    MS_df = pd.DataFrame.from_dict(distDict, orient="index", columns=sorted(distDict.keys()))
    MS_df.to_csv(prefix + ".MSDist.csv")
