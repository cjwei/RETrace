#!/usr/bin/env python
import argparse

'''
Usage: python script.py --probes probes.info.txt --output probes.bed
This script will take as input:
    1) probes.info.txt = contains all the probes designed per targetID
The script will then convert the probes.info.txt file into bed format with the following example format: (chromosomal position, period [subunit length], number of repeats, and name)
    chr1    6464771 6464818 4   12  MFD424-TTTA003
NOTE: For our probe files, we will replace the region name with chr1:6464771-6464818_12xNNNN
'''

#%%
def makeBed():
    parser = argparse.ArgumentParser(description="Convert probe file to bed format")
    parser.add_argument('--probes', action="store", dest="probe_file", help="File containing hybridization probe targets")
    parser.add_argument('--output', action="store", dest="output_file")
    args = parser.parse_args()
    
    output = open(args.output_file,"w")
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
            output.write(chrom + "\t" + chromStart + "\t" + chromEnd + "\t"
                         + str(len(sub_seq)) + "\t" + num_sub + "\t" + targetID + "\n")
    output.close()
    
#%%
if __name__ == "__main__":
    makeBed()