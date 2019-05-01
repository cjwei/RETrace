#!/usr/bin/env python3
import argparse
import pickle
from tqdm import tqdm

'''
Usage: python script.py --sample ... --gff ... --prefix ...
The script will calculate methRate across Ensembl Regulatory Build windows in order to perform NMF and tSNE plotting for single cells
To run this script, we must include the following input:
    1) sample_methDict = pre-computed methDict containing single cell sample information.  The structure of methDict is as follows:
    methDict
        "sample"
            sample_name
                "file_loc" = file_loc
        "base"
            base_loc (chr:pos;chain)
                sample_name = (num_meth, total_reads)
    2) gff = file containing location sof Ensembl Regulatory Build windows
    3) prefix = prefix for output files
The script will then go through and outpu the following files:
    1) prefix.methRate.csv = matrix containing methRate across Ensembl Regulatory Build
'''

def importReg(Ensembl_gff):
    regDict = {}
    #We want to import the Ensembl Regulatory Build windows from the gff file
    with open(Ensembl_gff) as f_gff:
        for line in f_gff:
            line_list = line.split()
            reg_name =

def main():
    parser = argparse.ArgumentParser(description="Calculate methRate across Ensembl Regulatory Build windows")
    parser.add_argument('--sample', action="store", dest="sample_pkl", help="Pre-computed methDict containing sample information")
    parser.add_argument('--gff', action="store", dest="Ensembl_gff", help="gff file containing Ensembl Regulatory Build windows")
    args = parser.parse_args()

    #Import pre-computed sample methDict
    print("Import sampleDict")
    with open(args.sample_pkl, 'rb') as sample_file:
        sampleDict = pickle.load(sample_file)

    #Import Ensemble Regulatory Build windows
    regDict = importReg(args.Ensembl_gff)

if __name__ == "__main__":
    main()
