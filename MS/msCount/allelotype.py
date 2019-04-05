#!/usr/bin/env python3
import argparse
import pickle
from importer import import_sampleDict
import peakutils #Similar to matlab findpeaks (as used in Biezuner, Gen Res 2006) <https://blog.ytotech.com/2015/11/01/findpeaks-in-python/>
from tqdm import tqdm
from count_visualizer import plotAlleles
from collections import Counter
from sklearn.utils import resample
import numpy as np
import math

'''
Usage: python script.py --input sample_list.txt --prefix output_prefix --allele_method [simple/stutter]
This script will go through and calculate allelotype from the given output_prefix.targetDict.pkl file, which has the following format:
    targetDict
        target_id
            "sample"
                sample
                    msCount_list = list containing all msCounts for given sample at target_id
We also want to input the same sample_list.txt file as used in msCount.py.  This will give us some useful information, particularly the sex of the cell, which will help with allelotyping.  Below are the information fields for each sample:
    1) File locaiton of sorted sample bam
    2) Sample name
    3) Sex (mostly used as information for allelotyping)
    4) [Optional] Sample type (i.e. known clone in ex vivo tree)
The allelotypes for all samples will be saved in the following output_prefix.alleleDict.pkl file:
    alleleDict
        target_id
            "sample"
                sample
                    allele_list = list containing allelotype for the given sample at target_id (will have either 1 or 2 alleles)
'''

'''Simple allelotyping by scipy peak-calling'''
def simple_typing(targetDict):
    '''Simple method for allelotyping by peak calling'''
    alleleDict = {}
    for target_id in tqdm(sorted(targetDict.keys())):
        alleleDict[target_id] = {}
        alleleDict[target_id]["sample"] = {}
        for sample in sorted(targetDict[target_id]["sample"].keys()):
            all_msCount = []  #List of all integer number subunits in msCount range
            for msCount in range(min(targetDict[target_id]["sample"][sample]), max(targetDict[target_id]["sample"][sample]) + 1):
                if msCount % len(targetDict[target_id]["sub_seq"]) == 0:
                    all_msCount.append(msCount)
            if len(all_msCount) > 2:
                msCount_freq = []
                for msCount in all_msCount:
                    msCount_freq.append(targetDict[target_id]["sample"][sample].count(msCount))
                #We want to perform peak-calling on msCount_freq and return the respective msCount call (allelotype) for each peak index
                allele_indx = peakutils.indexes(msCount_freq)
            elif len(all_msCount) == 2:
                if msCount_freq[0] > msCount_freq[1]:
                    allele_indx = [0]
                elif msCount_freq[0] < msCount_freq[1]:
                    allele_indx = [1]
                else:
                    allele_indx = [0, 1]
            elif len(all_msCount) == 1:
                allele_indx = [0]
            else:
                continue
            alleleDict[target_id]["sample"][sample] = [all_msCount[indx] for indx in allele_indx]
    return alleleDict

'''Stutter allelotyping by using method similar to LobSTR'''
def make_stutter(sampleDict, targetDict):
    '''Go through each sample file and calculate aggregate stutter information (num_modal, num_error, num_increase, xbar [number of non-unit deviations from modal count]) for each sub_len'''
    stutterDict = {}
    #Calculate initial stutter statistics for each hemizygous chrX/Y target per male sample
    for target_id in sorted(targetDict.keys()):
        if "chrX" in target_id or "chrY" in target_id: #We only want to calculate stutter from sex chromosomes
            for sample in sorted(targetDict[target_id]["sample"]):
                if sampleDict[sample]["sex"] == "M": #Only rely on hemizygous sex chromosomes in male samples
                    #We want to calculate the num_modal, num_error, xbar, and num_increase for each sample/target_id
                    sub_len = len(targetDict[target_id]["sub_seq"])
                    stutter_key = sub_len #Define stutter_key as sub_len (may use sub_seq in future for more detailed stutter model)
                    modal_msCount = Counter(targetDict[target_id]["sample"][sample]).most_common(1)[0][0] #Modal msCount for given target_id/sample (<https://stackoverflow.com/questions/10797819/finding-the-mode-of-a-list>)
                    num_modal = targetDict[target_id]["sample"][sample].count(modal_msCount)
                    num_error = len(targetDict[target_id]["sample"][sample]) - num_modal
                    xbar = 0
                    num_increase = 0
                    for msCount in set(targetDict[target_id]["sample"][sample]):
                        deviation = msCount - modal_msCount #This is the base difference from modal allele
                        xbar += (deviation % sub_len) * targetDict[target_id]["sample"][sample].count(msCount) #This keeps track of the number of non-unit deviations from modal count
                        if deviation > 0:
                            num_increase += targetDict[target_id]["sample"][sample].count(msCount)
                    #Save num_modal, num_error, xbar, and num_increase to stutterDict
                    if stutter_key not in stutterDict.keys():
                        stutterDict[stutter_key] = {}
                    else:
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
        bootstrap_samples = resample(sorted(sampleDict.keys()), replace=True, n_samples=len(sampleDict.keys()), random_state=1)
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
    for stutter_key in sorted(stutterDict.keys()):
        stutterDict[stutter_key]["final"] = {}
        stutterDict[stutter_key]["final"]["num_modal"] = np.mean(stutterDict[stutter_key]["bootstrap"]["num_modal"])
        stutterDict[stutter_key]["final"]["num_error"] = np.mean(stutterDict[stutter_key]["bootstrap"]["num_error"])
        stutterDict[stutter_key]["final"]["pi"] = np.mean(stutterDict[stutter_key]["bootstrap"]["pi"])
        stutterDict[stutter_key]["final"]["lambda"] = np.mean(stutterDict[stutter_key]["bootstrap"]["lambda"])
        stutterDict[stutter_key]["final"]["geom_p"] = np.mean(stutterDict[stutter_key]["bootstrap"]["geom_p"])
        stutterDict[stutter_key]["final"]["q"] = np.mean(stutterDict[stutter_key]["bootstrap"]["q"])
        print("----------\t" + str(stutter_key) + "\t----------")
        for param, value in stutterDict[stutter_key]["final"].items():
            print(param + "\t" + str(value))
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

def stutter_typing(sampleDict, targetDict, stutterDict, min_cov, min_ratio):
    alleleDict = {}
    for target_id in tqdm(sorted(targetDict.keys())):
        alleleDict[target_id] = {}
        alleleDict[target_id]["sample"] = {}
        for sample in sorted(targetDict[target_id]["sample"].keys()):
            #We want to check whether sample has sufficient number of reads/msCounts for given target
            msCount_list = targetDict[target_id]["sample"][sample]
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
            for msCount in range(min(targetDict[target_id]["sample"][sample]), max(targetDict[target_id]["sample"][sample]) + 1):
                if msCount % len(targetDict[target_id]["sub_seq"]) == 0:
                    all_msCount.append(msCount)
                    if msCount in msCount_dict.keys():
                        msCount_freq.append(msCount_dict[msCount])
                    else:
                        msCount_freq.append(0)
                        msCount_dict[msCount] = 0
            if len(all_msCount) > 2:
                #We want to perform peak-calling on mscount and return teh respective msCount call (possible allelotype) for each peak index
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
                alleleDict[target_id]["sample"][sample] = list(map(int, allelotype.split('/')))
    return alleleDict

def allelotype():
    parser = argparse.ArgumentParser(description="Allelotype all samples given msCounts per target_id")
    parser.add_argument('--input', action="store", dest="sample_info", help="Tab-delimited file containing sample information")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Specify output prefix (for targetDict and alleleDict, along with any stats or plot files)")
    parser.add_argument('--allele_method', action="store", dest="allele_method", default="simple", help="Specify allelotyping algorithm used [sample/stutter]")
    parser.add_argument('--min_cov', action="store", dest="min_cov", default=10, help="Specify minimum coverage for each target_id within each sample to call allelotype [default = 10]")
    parser.add_argument('--min_ratio', action="store", dest="min_ratio", default=0.0, help="Specify minimum percentae of reads supporting the resulting allelotype [default = 0.0; no filtering]")
    parser.add_argument('--plot_file', action="store", dest="plot_file", help="Flag for indicating whether we want to output a plot file (plot_file) visualzing msCounts per targetID")
    parser.add_argument('--num_plot', action="store", dest="num_plot", default=1784, help="If plotting allelotype, specify number of targets to plot")
    args = parser.parse_args()

    #Parse sample_info file
    sampleDict = import_sampleDict(args.sample_info)
    if not sampleDict:
        return

    #Import targetDict from pickle file
    targetDict = pickle.load(open(args.prefix + ".targetDict.pkl", "rb"))

    #Perform allelotyping
    if args.allele_method == "simple":
        print("Allelotyping using following method:\t" + args.allele_method)
        alleleDict = simple_typing(targetDict)
    elif args.allele_method == "stutter":
        print("Making stutter model")
        stutterDict = make_stutter(sampleDict, targetDict)
        print("Allelotyping using following method:\t" + args.allele_method)
        alleleDict = stutter_typing(sampleDict, targetDict, stutterDict, args.min_cov, float(args.min_ratio))
    else:
        print("Unknown allelotyping method:\t" + args.allele_method)
    pickle.dump(alleleDict, open(args.prefix + '.alleleDict.pkl', 'wb'))

    #Plot allelotype and msCounts
    if args.plot_file is not None:
        print("Plotting allelotypes/msCounts for each sample")
        plotAlleles(args.plot_file, int(args.num_plot), targetDict, alleleDict)

if __name__ == "__main__":
    allelotype()
