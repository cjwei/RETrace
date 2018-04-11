#!/usr/bin/env python3
import argparse
import collections
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

'''
Usage: python script.py input1.vcf input2.vcf... --output plots.pdf
This script will take as input multiple vcf files: (delimited by spaces)
    input1.vcf input2.vcf ...
It will then parse these vcf files and extract the allelotypes/msCounts in order to plot shared targetIDs to visualize noise in counts.
'''

def parseVCF(file_list):
    alleleDict = {}
    for f in file_list:
        for line in f:
            if not line.startswith('#'):
                vcf_line = line.split()
                targetID = vcf_line[2]
                info_list = vcf_line[-1].split(':')
                if len(info_list)==1:
                    continue
                '''
                We expect the genotype info to have the following format:
                    GT:GB:Q:PQ:DP:DSNP:DSTUTTER:DFLANKINDEL:PDP:PSNP:GLDIFF:AB:DAB:FS:ALLREADS:MALLREADS
                    1|2:-2|10:1.00:0.50:40:0:0:0:20.50|19.50:0|0:8.87:-0.00:13:0.00:-2|6;0|1;10|2:-2|5;10|2
                '''
                allelotype = info_list[1].split('|')
                allele_list = []
                count_list = []
                for allele_info in info_list[-1].split(';'):
                    [allele,count] = allele_info.split('|')
                    allele_list.append(int(allele))
                    count_list.append(int(count))
                if targetID not in alleleDict.keys():
                    alleleDict[targetID] = {}
                alleleDict[targetID][f.name] = {}
                alleleDict[targetID][f.name]["allelotype"] = list(map(int,allelotype))
                alleleDict[targetID][f.name]["allele_list"] = allele_list
                alleleDict[targetID][f.name]["count_freq"] = [float(num)/float(sum(count_list)) for num in count_list]
    return alleleDict

def plotVCF(alleleDict, file_list, plot_file):
    plots = PdfPages(plot_file)
    for targetID in sorted(alleleDict.keys()):
        for file_name in sorted(alleleDict[targetID].keys()):
            #Plot raw alleles/msCounts (blue) along with final allelotypes (red)
            ax = plt.gca()
            ax.set_ylim([0,1])
            ax.set_xlim([min(alleleDict[targetID][file_name]["allele_list"])-5,max(alleleDict[targetID][file_name]["allele_list"])+5])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xlabel("Base pair differences of alleles from reference")
            plt.plot(alleleDict[targetID][file_name]["allelotype"],[0.5]*len(alleleDict[targetID][file_name]["allelotype"]),'ro')
            plt.vlines(alleleDict[targetID][file_name]["allele_list"],[0],alleleDict[targetID][file_name]["count_freq"],linestyle="dashed",color="b")
            plt.title(targetID + ", " + file_name)

            if len(alleleDict[targetID].keys()) == len(file_list):
                ax.set_facecolor('xkcd:light grey')
            plt.savefig(plots, format='pdf')
            plt.clf()

    plots.close()
    return

def main():
    parser = argparse.ArgumentParser(description="Parse HipSTR vcf output and visulaize msCounts/allelotype per targetID shared across files")
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    parser.add_argument('--output', action="store", dest="plot_file", help="Specify output plot_file")
    args = parser.parse_args()

    #Import targetDict from all vcf files
    alleleDict = parseVCF(args.file)

    #Plot msCounts for each targetID shared among all files
    plotVCF(alleleDict, args.file, args.plot_file)

#%%
if __name__ == "__main__":
    main()
