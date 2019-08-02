#!/usr/bin/env python3
import argparse
import pickle
from tqdm import tqdm
import re
import gzip

def compMethyl(infoDict, sampleDict, prefix):
    f_out = open(prefix + ".compMethyl.txt", 'w')
    for base_loc in tqdm(sorted(sampleDict["base"].keys())):
        #Only consider base_loc that have at least one test and one control sample covered
        control_samples = set(sampleDict["base"][base_loc].keys()).intersection(set(infoDict["group"]["control"]))
        test_samples = set(sampleDict["base"][base_loc].keys()).intersection(set(infoDict["group"]["test"]))
        if len(control_samples) > 0 and len(test_samples) > 0:
            #We want to calculate various methyl stats for each base_loc
            control_meth = 0
            control_total = 0
            hypo_list = []
            hypo_methRate = []
            hypo_PD = []
            hyper_list = []
            hyper_methRate = []
            hyper_PD = []
            for sample_name in control_samples: #We want to first gather statistics for control group in order to assign hypo/hyper samples
                control_meth += sampleDict["base"][base_loc][sample_name][0]
                control_total += sampleDict["base"][base_loc][sample_name][1]
            control_methRate = float(control_meth / control_total)
            for sample_name in test_samples:
                test_methRate = float(sampleDict["base"][base_loc][sample_name][0] / sampleDict["base"][base_loc][sample_name][1])
                if test_methRate > control_methRate: #Indicates hyper-methylated sample
                    hyper_list.append(sample_name)
                    hyper_methRate.append(test_methRate)
                    hyper_PD.append(abs(test_methRate - control_methRate) * 100)
                elif test_methRate < control_methRate: #Indicates hypo-methylated sample
                    hypo_list.append(sample_name)
                    hypo_methRate.append(test_methRate)
                    hypo_PD.append(abs(test_methRate - control_methRate) * 100)
            # f_out.write(base_loc + "\t" + str(round(control_methRate,2)) + "\t" + \
            #     str(len(hypo_list)) + "\t" + ','.join(hypo_list) + "\t" + ','.join(str(round(i,2)) for i in hypo_methRate) + "\t" + ','.join(str(round(i,2)) for i in hypo_PD) + "\t" + \
            #     str(len(hyper_list)) + "\t" + ','.join(hyper_list) + "\t" + ','.join(str(round(i,2)) for i in hyper_methRate) + "\t" + ','.join(str(round(i,2)) for i in hyper_PD)+ "\n")
            if len(hypo_list) == 0:
                avg_hypo_methRate = "NA"
                avg_hypo_PD = "NA"
            else:
                avg_hypo_methRate = round(sum(hypo_methRate) / len(hypo_methRate), 4)
                avg_hypo_PD = round(sum(hypo_PD) / len(hypo_PD), 4)
            if len(hyper_list) == 0:
                avg_hyper_methRate = "NA"
                avg_hyper_PD = "NA"
            else:
                avg_hyper_methRate = round(sum(hyper_methRate) / len(hyper_methRate), 4)
                avg_hyper_PD = round(sum(hyper_PD) / len(hyper_PD), 4)
            (chrom, chromStart, strand) = re.split(':|;', base_loc)
            chromEnd = str(int(chromStart) + 1)
            f_out.write("chr" + chrom + "\t" + chromStart + "\t" + chromEnd + "\t" + str(round(control_methRate, 4)) + "\t" + \
                str(len(hypo_list)) + "\t" + ','.join(hypo_list) + "\t" + str(avg_hypo_methRate) + "\t" + str(avg_hypo_PD) + "\t" + \
                str(len(hyper_list)) + "\t" + ','.join(hyper_list) + "\t" + str(avg_hyper_methRate) + "\t" + str(avg_hyper_PD) + "\n")
    f_out.close()

def main():
    '''
    This script will take as input a sample_info file that gives details on the samples we want to compare:
        sampleName1    test
        sampleName2    test
        sampleName3    control (this may be a reference cell line like HCT116)
        ...
    The script will then merge all control samples into a pseudo-bulk methylation signal by combining all methylation counts for each covered base position.
    We'll then compare all samples from test to this pseudo-bulk in order to determine whether there is any shared methylation sites that are different between samples and bulk.
    The output will be a file that contains:
        base_loc    control_methRate   num_hypo     hypo_samples    PD_hypo     num_hyper       hyper_samples   PD_hyper
    '''
    parser = argparse.ArgumentParser(description="Compare methyl signal from single cells against pseudo-bulk")
    parser.add_argument("--sample_info", action="store", dest="sample_info", help="Specify sample_info that contains the single cells and bulk to analyze")
    parser.add_argument("--methDict", action="store", dest="sample_methDict", help="methDict containing all meth calls for analysis")
    parser.add_argument("--prefix", action="store", dest="prefix", help="Output prefix")
    args = parser.parse_args()

    infoDict = {}
    infoDict["group"] = {}
    with open(args.sample_info) as f_info:
        for line in f_info:
            (sample_name, sample_group) = line.split()
            if sample_group not in infoDict["group"].keys():
                infoDict["group"][sample_group] = []
            infoDict["group"][sample_group].append(sample_name)

    #Import pre-computed sample methDict
    print("Import sampleDict")
    with open(args.sample_methDict, 'rb') as sample_methDict_file:
        sampleDict = pickle.load(sample_methDict_file)

    #Import methyl calls and asign to appropriate sample_group/sample_name
    compMethyl(infoDict, sampleDict, args.prefix)

if __name__ == "__main__":
    main()
