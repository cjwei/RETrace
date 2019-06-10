#!/usr/bin/env python3
from tqdm import tqdm
import re
import pandas as pd

def make_MethylDist(filtered_samples, methDict, regDict, prefix):
    '''
    Calculate Methyl distance between two samples and return distance matrix
    '''

    #We want to pre-filter only base loc that belong within regDict
    print("\tPre-filtering only bases found in Ensemble regulatory windows")
    for base_loc in tqdm(sorted(methDict["base"].keys())):
        (chr, pos, chain) = re.split(':|;', base_loc)
        if 'gl' in chr: #We need to make a special condition for chromosomes that contain 'gl'
            chr = chr.split('_')[1]
        if chr in regDict.keys():
            if int(pos) not in regDict[chr].keys():
                del methDict["base"][base_loc]
            else:
                for sample_name in sorted(methDict["base"][base_loc].keys()):
                    if sample_name not in filtered_samples: #We only want to analyze sample_name in sample_list
                        del methDict["base"][base_loc][sample_name]

    PD_dict = {} #Keep all pairwise comparisons
    #Next, we want to compare each individual CpGs shared between two samples
    print("\tCalculating pairwise dissimilarity between two samples")
    for base_loc in tqdm(sorted(methDict["base"].keys())):
        analyzed_pairs = []
        for sample1 in sorted(methDict["base"][base_loc].keys()):
            for sample2 in sorted(methDict["base"][base_loc].keys()):
                sample_pair = tuple(sorted([sample1, sample2]))
                if sample_pair not in analyzed_pairs:
                    analyzed_pairs.append(sample_pair)
                    if sample_pair not in PD_dict.keys():
                        PD_dict[sample_pair] = []
                    methRate1 = methDict["base"][base_loc][sample_pair[0]][0] / methDict["base"][base_loc][sample_pair[0]][1]
                    methRate2 = methDict["base"][base_loc][sample_pair[1]][0] / methDict["base"][base_loc][sample_pair[1]][1]
                    dis = abs(methRate1 - methRate2) * 100
                    PD_dict[sample_pair].append(dis)

    #Calculate pairwise distance for all sample pairs
    distDict = {}
    for sample1 in sorted(filtered_samples):
        distDict[sample1] = []
        for sample2 in sorted(filtered_samples):
            sample_pair = tuple(sorted([sample1, sample2]))
            Methyl_dist = float(sum(PD_dict[sample_pair]) / len(PD_dict[sample_pair]))
            distDict[sample1].append(Methyl_dist)

    Methyl_df = pd.DataFrame.from_dict(distDict, orient="index", columns=sorted(distDict.keys()))
    Methyl_df.to_csv(prefix + ".MethylDist.csv")
