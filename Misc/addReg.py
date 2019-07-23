#!/usr/bin/env python3
import argparse
import gzip
from tqdm import tqdm
import pickle
import os

'''
Usage: python script.py --tsv ... --Reg /media/Home_Raid1/cjwei/software/RETrace/Data/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz --output ...
This script will take as input the reference tsv file created from methylpy (may also have DMR already annotated):
    10      94041   -       CGG     2       2       1
    10      94422   -       CGC     1       1       1
    10      94427   -       CGA     1       1       1
    10      94882   +       CGG     1       1       1       DMR:chr10:94882-95053
    10      94901   +       CGC     1       1       1       DMR:chr10:94882-95053
    10      94910   +       CGC     1       1       1       DMR:chr10:94882-95053
    ...
It'll then add regulator information to teh bases that occur within regulatory build from Ensembl
    GL000195.1      Regulatory_Build        TF_binding_site 51049   51622   .       .       .       ID=TF_binding_site:ENSR00000976637;bound_end=51622;bound_start=51029;description=Transcription factor binding site;feature_type=TF binding site
    GL000195.1      Regulatory_Build        TF_binding_site 54279   54721   .       .       .       ID=TF_binding_site:ENSR00000976638;bound_end=54865;bound_start=54279;description=Transcription factor binding site;feature_type=TF binding site
    ...
'''

def importReg(reg_file, reg_pkl):
    regDict = {}
    #We want to import the regulatory build window information into regDict
    with gzip.open(reg_file, 'rt') as f_reg:
        for i, line in enumerate(tqdm(f_reg)):
            (chr, source, feature, start, end, score, strand, frame, attribute) = line.split("\t")
            reg_id = attribute.split(';')[0].split(':')[-1]
            reg_name = "reg:" + reg_id + ':' + chr + ':' + start + '-' + end
            if chr.replace('chr','') not in regDict.keys():
                regDict[chr.replace('chr','')] = {}
            for pos in range(int(start), int(end) + 1): #We want to iterate through all base positions within reg window
                regDict[chr.replace('chr','')][pos] = reg_name
    pickle.dump(regDict, open(reg_pkl, 'wb'))
    return regDict

def filterReg(tsv_file, regDict, output_file):
    tsv_output = open(output_file, 'w')
    with open(tsv_file) as f_tsv:
        for i, line in enumerate(tqdm(f_tsv)):
            (chr, pos) = line.split()[0:2]
            if chr.replace('chr','') in regDict.keys():
                if int(pos) in regDict[chr.replace('chr','')].keys():
                    tsv_output.write(line.rstrip() + "\t" + regDict[chr.replace('chr','')][int(pos)] + "\n")
            else:
                tsv_output.write(line) #We want to keep the CpG information line in the file, but without label of Reg.  This allows us to have options for downstream analysis
    return

def main():
    parser = argparse.ArgumentParser(description="Annotate methlyation singal found in Ensembl regulatory build windows")
    parser.add_argument('--tsv', action="store", dest="tsv_file", help="Original methylpy tsv file for all CGN")
    parser.add_argument('--reg_file', action="store", dest="reg_file", help="Gff file containing Ensembl regulatory build information")
    parser.add_argument('--reg_pkl', action="store", dest="reg_pkl", help="Specify if regDict pickle file is already made")
    parser.add_argument('--output', action="store", dest="output_file", help="Filename for output filtered tsv file")
    args = parser.parse_args()

    if os.path.exists(args.reg_pkl):
        regDict = pickle.load(open(args.reg_pkl, 'rb'))
    else:
        regDict = importReg(args.reg_file, args.reg_pkl)
    filterReg(args.tsv_file, regDict, args.output_file)

if __name__ == "__main__":
    main()
