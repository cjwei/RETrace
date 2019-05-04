#!/usr/bin/env python3
import peakutils #Similar to matlab findpeaks (as used in Biezuner, Gen Res 2006) <https://blog.ytotech.com/2015/11/01/findpeaks-in-python/>
from tqdm import tqdm
from collections import Counter
from sklearn.utils import resample
import numpy as np
import math

def make_stutter(sampleDict, targetDict, alleleDict, prefix):
    '''
    Go through each sample file and calculate aggregate stutter information (num_modal, num_error, num_increase, xbar [number of non-unit deviations from modal count]) for each sub_len
    '''
    stutterDict = {}
    #Calculate initial stutter statistics for each hemizygous chrX/Y target per male sample
    for target_id in sorted(alleleDict.keys()):
        if "chrX" in target_id or "chrY" in target_id: #We only want to calculate stutter from sex chromosomes
            for sample in sorted(alleleDict[target_id]["sample"].keys()):
                if sampleDict[sample]["sex"] == "M": #Only rely on hemizygous sex chromosomes in male samples
                    #We want to calculate the num_modal, num_error, xbar, and num_increase for each sample/target_id
                    sub_len = len(targetDict[target_id]["sub_seq"])
                    stutter_key = sub_len #Define stutter_key as sub_len (may use sub_seq in future for more detailed stutter model)
                    modal_msCount = Counter(alleleDict[target_id]["sample"][sample]["msCount"]).most_common(1)[0][0] #Modal msCount for given target_id/sample (<https://stackoverflow.com/questions/10797819/finding-the-mode-of-a-list>)
                    num_modal = alleleDict[target_id]["sample"][sample]["msCount"].count(modal_msCount)
                    num_error = len(alleleDict[target_id]["sample"][sample]["msCount"]) - num_modal
                    xbar = 0
                    num_increase = 0
                    for msCount in set(alleleDict[target_id]["sample"][sample]["msCount"]):
                        deviation = msCount - modal_msCount #This is the base difference from modal allele
                        xbar += (deviation % sub_len) * alleleDict[target_id]["sample"][sample]["msCount"].count(msCount) #This keeps track of the number of non-unit deviations from modal count
                        if deviation > 0:
                            num_increase += alleleDict[target_id]["sample"][sample]["msCount"].count(msCount)
                    #Save num_modal, num_error, xbar, and num_increase to stutterDict
                    if stutter_key not in stutterDict.keys():
                        stutterDict[stutter_key] = {}
                    if sample not in stutterDict[stutter_key].keys():
                        stutterDict[stutter_key][sample] = {}
                        stutterDict[stutter_key][sample]["num_modal"] = num_modal
                        stutterDict[stutter_key][sample]["num_error"] = num_error
                        stutterDict[stutter_key][sample]["xbar"] = xbar
                        stutterDict[stutter_key][sample]["num_increase"] = num_increase
                        stutterDict[stutter_key][sample]["sub_len"] = sub_len
                    else:
                        stutterDict[stutter_key][sample]["num_modal"] += num_modal
                        stutterDict[stutter_key][sample]["num_error"] += num_error
                        stutterDict[stutter_key][sample]["xbar"] += xbar
                        stutterDict[stutter_key][sample]["num_increase"] += num_increase
    #We want to do bootstrap resampling across all samples for each stutter_key in order to more accurately calculate stutter
    #We first initialize the necessary lists for storing the various stutter parameters from each bootstrap resampling
    for stutter_key in sorted(stutterDict.keys()):
        stutterDict[stutter_key]["bootstrap"] = {}
        stutterDict[stutter_key]["bootstrap"]["num_modal"] = []
        stutterDict[stutter_key]["bootstrap"]["num_error"] = []
        stutterDict[stutter_key]["bootstrap"]["pi"] = []
        stutterDict[stutter_key]["bootstrap"]["lambda"] = []
        stutterDict[stutter_key]["bootstrap"]["geom_p"] = []
        stutterDict[stutter_key]["bootstrap"]["q"] = []
    for i in range(1000): #We want to do bootstrap resampling of samples 1000 times
        bootstrap_samples = resample(sorted(sampleDict.keys()), replace=True, n_samples=len(sampleDict.keys()))
        stutterParam = {} #Make temperary dictionary containing stutter parameters for each bootstrap run
        for sample in bootstrap_samples:
            for stutter_key in sorted(stutterDict.keys()):
                if sample in stutterDict[stutter_key].keys():
                    if stutter_key not in stutterParam.keys():
                        stutterParam[stutter_key] = {}
                        stutterParam[stutter_key]["num_modal"] = stutterDict[stutter_key][sample]["num_modal"]
                        stutterParam[stutter_key]["num_error"] = stutterDict[stutter_key][sample]["num_error"]
                        stutterParam[stutter_key]["xbar"] = stutterDict[stutter_key][sample]["xbar"]
                        stutterParam[stutter_key]["num_increase"] = stutterDict[stutter_key][sample]["num_increase"]
                        stutterParam[stutter_key]["sub_len"] = stutterDict[stutter_key][sample]["sub_len"] #We want to define a sub_len variable, which should be constant for given stutter_key
                    else:
                        stutterParam[stutter_key]["num_modal"] += stutterDict[stutter_key][sample]["num_modal"]
                        stutterParam[stutter_key]["num_error"] += stutterDict[stutter_key][sample]["num_error"]
                        stutterParam[stutter_key]["xbar"] += stutterDict[stutter_key][sample]["xbar"]
                        stutterParam[stutter_key]["num_increase"] += stutterDict[stutter_key][sample]["num_increase"]
        #With the necessary statistics from a combination of the resampled files, we want to then calculate the stutter model parameters for the given bootstrap resampling
        for stutter_key in sorted(stutterDict.keys()):
            #0) We want to skip any sub_len that have too few reads to accurately call stutter
            if (stutterParam[stutter_key]["num_modal"] + stutterParam[stutter_key]["num_error"]) < 100:
                continue
            stutterDict[stutter_key]["bootstrap"]["num_modal"].append(int(stutterParam[stutter_key]["num_modal"]))
            stutterDict[stutter_key]["bootstrap"]["num_error"].append(int(stutterParam[stutter_key]["num_error"]))
            #1) We need the probability that a read is a product of stutter noise
            stutterDict[stutter_key]["bootstrap"]["pi"].append(float(stutterParam[stutter_key]["num_error"]/(stutterParam[stutter_key]["num_modal"] + stutterParam[stutter_key]["num_error"])))
            #2) We want the lambda variable for the poisson distribution used to determien the probability D(e,m) that noisy read deviates by e bases from teh original allele.  This is merely the length of the subunit
            stutterDict[stutter_key]["bootstrap"]["lambda"].append(int(stutterParam[stutter_key]["sub_len"]))
            #3) We want the probability p value (and hence the xbar value since p=1/(xbar+1)) that will be used for the geometric distribution for D(e,m)
            stutterDict[stutter_key]["bootstrap"]["geom_p"].append(float(1/(stutterParam[stutter_key]["xbar"]/stutterParam[stutter_key]["num_error"] + 1)))
            #4) We finally want to define q for D(e,m), which is the probability that stutter will increase the true allele length.  This is simply the division of num_increase/num_error
            stutterDict[stutter_key]["bootstrap"]["q"].append(float(stutterParam[stutter_key]["num_increase"]/stutterParam[stutter_key]["num_error"]))
    #We want to take the average of all these values in the bootstrap sub-dictionary and output them into the final stutterDict
    f_stutter = open(prefix + '.Custom.stutter.txt', 'w')
    for stutter_key in sorted(stutterDict.keys()):
        stutterDict[stutter_key]["final"] = {}
        stutterDict[stutter_key]["final"]["num_modal"] = np.mean(stutterDict[stutter_key]["bootstrap"]["num_modal"])
        stutterDict[stutter_key]["final"]["num_error"] = np.mean(stutterDict[stutter_key]["bootstrap"]["num_error"])
        stutterDict[stutter_key]["final"]["pi"] = np.mean(stutterDict[stutter_key]["bootstrap"]["pi"])
        stutterDict[stutter_key]["final"]["lambda"] = np.mean(stutterDict[stutter_key]["bootstrap"]["lambda"])
        stutterDict[stutter_key]["final"]["geom_p"] = np.mean(stutterDict[stutter_key]["bootstrap"]["geom_p"])
        stutterDict[stutter_key]["final"]["q"] = np.mean(stutterDict[stutter_key]["bootstrap"]["q"])
        f_stutter.write("----------\t" + str(stutter_key) + "\t----------\n")
        for param, value in stutterDict[stutter_key]["final"].items():
            f_stutter.write(param + "\t" + str(value) + "\n")
    f_stutter.close()
    return stutterDict

'''Allelotyping by using stutter model'''
def calcProb(e, stutterDict, sub_len): #This subroutine will calculate teh probability of observing the given UMI MS length (L) given the allele length (A) in question.  Consequently, it will take as input e = L - A, stutterDict, and sub_len
    if e == 0:
        prob = 1 - stutterDict[sub_len]["final"]["pi"]
    else:
        #We want to calculate the Poisson distribution.  Pois(stutterDict[sub_len]["final"]["lambda"])
        prob = (stutterDict[sub_len]["final"]["lambda"]**abs(e)) * math.exp(-stutterDict[sub_len]["final"]["lambda"]) / math.factorial(abs(e))
        #We next want to multiply by the geometric distribution of non-unit step sizes modulo sub_len
        prob *= stutterDict[sub_len]["final"]["geom_p"] * (1 - stutterDict[sub_len]["final"]["geom_p"])**(abs(e) % sub_len)
        #We want to account for increase vs decrease in length
        if e > 0:
            prob *= stutterDict[sub_len]["final"]["q"]
        else:
            prob *= (1 - stutterDict[sub_len]["final"]["q"])
        #We finally want to multiply by pi (the probaiblity a read results from stutter noise)
        prob *= stutterDict[sub_len]["final"]["pi"]
    return prob

def calcLogLike(alleleA, alleleB, msCount_list, stutterDict, sub_len): #We want to calculate the log likelihood of a given allelotype A/B
    logLike = 0
    for msCount in sorted(msCount_list):
        probA = calcProb(msCount - alleleA, stutterDict, sub_len)
        probB = calcProb(msCount - alleleB, stutterDict, sub_len)
        logLike += math.log10(max(probA, probB)) #We want the log10 of max(probA, probB)
    return logLike

def stutter_typing(sampleDict, targetDict, alleleDict, stutterDict, min_cov, min_ratio):
    for target_id in tqdm(sorted(alleleDict.keys())):
        for sample in sorted(alleleDict[target_id]["sample"]):
            msCount_list = alleleDict[target_id]["sample"][sample]["msCount"]
            if len(msCount_list) < min_cov:
                continue
            #Gather msCount frequency into msCount_dict where keys are microsatellite count and values are the number of times they occur
            msCount_dict = dict(Counter(msCount_list))
            logLikeDict = {} #We want to save all log likelihoods for allelotypes into dictionary
            likeDict = {} #We want to save all likelihoods into dictionary
            alleleLikeDict = {} #We want to save allele marginal likelihoods
            #We want to perform peak calling in order to determine most possible msCounts that could be alleles.  This is because we do not want to consider msCounts that are not local-maxima of count frequency, which most likely is not allele
            all_msCount = [] #List of all integer number subunits in msCount range for given sample
            msCount_freq = [] #List of number of reads for each of the above msCount in all_msCount
            for msCount in range(min(msCount_list), max(msCount_list) + 1):
                if msCount % len(targetDict[target_id]["sub_seq"]) == 0:
                    all_msCount.append(msCount)
                    if msCount in msCount_dict.keys():
                        msCount_freq.append(msCount_dict[msCount])
                    else:
                        msCount_freq.append(0)
                        msCount_dict[msCount] = 0
            if len(all_msCount) > 2:
                #We want to perform peak-calling on mscount and return the respective msCount call (possible allelotype) for each peak index
                allele_indx = peakutils.indexes(msCount_freq)
            else: #If number of msCounts is <= 2, we want to consider all given msCount as possible alleles
                allele_indx = range(len(all_msCount))
            possible_alleles = [all_msCount[indx] for indx in allele_indx]
            #We now want to calculate which of the possible_alleles is the most likely allelotype for given target/sample given all msCount
            for alleleA in sorted(possible_alleles):
                if alleleA not in alleleLikeDict.keys():
                    alleleLikeDict[alleleA] = 0
                for alleleB in sorted(possible_alleles):
                    if alleleB not in alleleLikeDict.keys():
                        alleleLikeDict[alleleB] = 0
                    if alleleB < alleleA:
                        continue #We do not want to repeat analysis on a poential allelelotype
                    elif alleleA == alleleB:
                        read_cov = msCount_dict[alleleA] #We want to keep track of the number of reads that support allelotype
                    else:
                        read_cov = msCount_dict[alleleA] + msCount_dict[alleleB]
                        if ("chrX" in target_id and sampleDict[sample]["sex"] == 'M') or ("chrY" in target_id and sampleDict[sample]["sex"] == 'M'):
                            continue #We only want to consider homozygous allelotypes if we are looking at male hemizygous sex chromosomes
                    if read_cov >= (min_ratio * len(msCount_list)): #We want to set a cutoff of the minimum ratio of total reads supporting the resulting allelotypes to min_ratio
                        like_key = '/'.join(str(allele) for allele in sorted([alleleA, alleleB]))
                        logLikeDict[like_key] = calcLogLike(alleleA, alleleB, msCount_list, stutterDict, len(targetDict[target_id]["sub_seq"]))
            if sum(logLikeDict.values()) != 0:
                allelotype = max(logLikeDict, key = logLikeDict.get)
                alleleDict[target_id]["sample"][sample]["allelotype"] = list(map(int, allelotype.split('/')))
    return alleleDict

def allelotype(sampleDict, prefix, targetDict, alleleDict, min_cov, min_ratio):
    '''
    This script will go through and calculate allelotype from the given alleleDict containing msCounts
    The allelotypes for all samples will be saved in the following alleleDict structure:
        alleleDict
            target_id (from targetDict)
                "sample"
                    sample (from sampleDict, which is already defined when labeling readGroups prior to HipSTR)
                        "msCount"
                            list of msCounts
                        "allelotype" - New
                            list of alleles (2 alleles) - New
    '''
    print("Making stutter model")
    stutterDict = make_stutter(sampleDict, targetDict, alleleDict, prefix)
    print("Allelotyping using stutter model")
    alleleDict = stutter_typing(sampleDict, targetDict, alleleDict, stutterDict, min_cov, min_ratio)

    return alleleDict
