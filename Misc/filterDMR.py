#!/usr/bin/env python3
import argparse

def importDMR(DMR_file):
    DMRdict = {}
    with open(DMR_file) as f_DMR:
        next(f_DMR) #We want to skip the first line of the file, which contains column headers
        for line in f_DMR:
            info_list = line.split("\t")
            chr = info_list[0]
            chromStart = info_list[1]
            chromEnd = info_list[2]
            DMR_region = chr + ":" + chromStart + "-" + chromEnd
            DMRdict[DMR_region] = {}
            DMRdict[DMR_region]["bed_info"] = ["chr" + chr, chromStart, chromEnd]
            DMRdict[DMR_region]["hyper_samples"] = set(info_list[4].split(','))
            DMRdict[DMR_region]["hypo_samples"] = set(info_list[5].split(','))
    return DMRdict

def filterDMR(DMRdict, sample_list, prefix):
    f_hyper = open(prefix + ".hyperDMR.txt", 'w')
    f_hypo = open(prefix + ".hypoDMR.txt", 'w')
    for DMR_region in sorted(DMRdict.keys()):
        hyper_jaccard = len(sample_list.intersection(DMRdict[DMR_region]["hyper_samples"])) / len(sample_list.union(DMRdict[DMR_region]["hyper_samples"]))
        hypo_jaccard = len(sample_list.intersection(DMRdict[DMR_region]["hypo_samples"])) / len(sample_list.union(DMRdict[DMR_region]["hypo_samples"]))
        f_hyper.write("\t".join(DMRdict[DMR_region]["bed_info"]) + "\t" + str(hyper_jaccard) + "\n")
        f_hypo.write("\t".join(DMRdict[DMR_region]["bed_info"]) + "\t" + str(hypo_jaccard) + "\n")
    f_hyper.close()
    f_hypo.close()
    return

def main():
    parser = argparse.ArgumentParser(description="Filter DMR regions that are specific to specified samples")
    parser.add_argument('--DMR', action="store", dest="DMR_file", help="tsv file output from methylpy DMR find")
    parser.add_argument('--sample_list', action="store", dest="sample_list", help="File containing samples of interest per line")
    parser.add_argument('--prefix', action="store", dest="prefix", help="Prefix for output files")
    args = parser.parse_args()

    DMRdict = importDMR(args.DMR_file) #Import all DMR
    sample_list = set(open(args.sample_list).read().splitlines()) #Import samples of interest
    filterDMR(DMRdict, sample_list, args.prefix)

if __name__ == "__main__":
    main()
