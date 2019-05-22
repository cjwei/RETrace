#!/usr/bin/env python3
import argparse

'''
Usage: python script.py --DRM_tsv cell_line_DMR_rms_results_collapsed.tsv
This script will take as input the DMR results tsv file from methylpy and export all DMR bed files for each cell line in file.  This is used specifically for parsing the cell line methylpy DMR output file
The input format of the tsv file is as follows:
    #chr    start   end number_of_dms   hypermethylated_samples hypomethylated_samples  methylation_level_cellLine1 methylation_level_cellLine2 ...
It will then output cellLine.DMR.bed files for each cell line described in the tsv file.
NOTE: We want to include both hyper and hypomethylated DMR regions for each cell line, which is likewise used in the Luo, Science 2017 paper
'''

def makeDMRbed(tsv_file):
    DMRdict = {} #This dictionary contains cell lines as keys and list of DMRs found in that sample as values (hyper- and hypomethylated samples)
    with open(tsv_file) as f_tsv:
        for line in f_tsv:
            if line.startswith('#'):
                continue
            else:
                DMR_line = line.split("\t")
                DMR_region = 'chr' + DMR_line[0] + "\t" + DMR_line[1] + "\t" + DMR_line[2]
                hyper_samples = DMR_line[4].split(',')
                hypo_samples = DMR_line[5].split(',')
                all_samples = hyper_samples + hypo_samples
                for sample in sorted(all_samples):
                    if sample is "": #we want to filter out empty sample names, which occur when no samples were classified as either hyper or hypo methylated
                        continue
                    if sample in DMRdict.keys():
                        DMRdict[sample].append(DMR_region)
                    else:
                        DMRdict[sample] = [DMR_region]
    for sample in sorted(DMRdict.keys()):
        f_bed = open(sample + '.DMR.bed', 'w')
        f_bed.write("\n".join(sorted(DMRdict[sample])))
        f_bed.close()

def main():
    parser = argparse.ArgumentParser(description="Filter methlyation signal found in DMR regions")
    parser.add_argument('--DMR_tsv', action="store", dest="tsv_file", help="Methylpy DMR tsv file containing DMR regions for each cell line")
    args = parser.parse_args()

    makeDMRbed(args.tsv_file)

if __name__ == "__main__":
    main()
