#!/usr/bin/env python3
import argparse

'''
Usage: python script.py --probes probes.info.txt --prefix probes
This script will take as input:
    1) probes.info.txt = contains all the probes designed per targetID
The script will then convert the probes.info.txt file into two bed formats (HipSTR and Picard) with the following example formats:
    prefix.HipSTR.bed
        (chromosomal position, period [subunit length], number of repeats, and name)
        chr1    6464771 6464818 4   12  chr1:6464771-6464818_12xNNNN
    prefix.Picard.bed
        (chromosoml position, name)
        chr1    6464771 6464818 chr1:6464771-6464818_12xNNNN
'''

#%%
def makeBed():
    parser = argparse.ArgumentParser(description="Convert probe file to bed format")
    parser.add_argument('--probes', action="store", dest="probe_file", help="File containing hybridization probe targets")
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    output_HipSTR = open(args.prefix + ".HipSTR.bed", 'w')
    output_Picard = open(args.prefix + ".Picard.bed", 'w')
    with open(args.probe_file) as f:
        for line in f:
            target_info = line.split(',')[0]
            info_list = target_info.split('_')
            #Obtain chromosomal position
            chrom = "_".join(info_list[0:len(info_list)-3])
            chromStart = info_list[-3]
            chromEnd = info_list[-2]

            #Obtain subunit information and reformat targetID
            (num_sub,sub_seq) = info_list[-1].split('x')
            targetID = chrom + ":" + chromStart + "-" + chromEnd + "_" + info_list[-1]


            #Print relevant information in bed format
            output_HipSTR.write(chrom + "\t" + chromStart + "\t" + chromEnd + "\t"
                         + str(len(sub_seq)) + "\t" + num_sub + "\t" + targetID + "\n")
            output_Picard.write(chrom + "\t" + chromStart + "\t" + chromEnd + "\t" + targetID + "\n")
    output_HipSTR.close()
    output_Picard.close()

#%%
if __name__ == "__main__":
    makeBed()
