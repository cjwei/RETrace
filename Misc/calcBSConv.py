#!/usr/bin/env python3
import argparse

'''
This script will take as input tsv files from methylpy against lambda genome in order to calculate the bisulfite conversion efficiency for various samples.
The tsv files are formatted as follows:
    chr pos strand  seq_context meth_count  read_cov ...
    Lambda_NEB	1904	+	CGG	0	1	1
    Lambda_NEB	1910	+	CGA	0	1	1
    Lambda_NEB	1913	+	CTC	0	1	1
Because we know that the lambda dna spiked into each reaction should be fully unmethylated, we want to calculate the number of C's that were counted as methylated (based off of a simple majority of reads)
'''

def main():
    parser = argparse.ArgumentParser(description="Determine bisulfite conversion efficiency from methylpy on lambda DNA")
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+', help="List of tsv files to be analyzed")
    parser.add_argument('--output', action="store", dest="output_file", help="Filename for output bisulfite conversion efficiency")
    args = parser.parse_args()

    f_out = open(args.output_file, 'w')

    for f in args.file:
        num_bases = 0
        num_meth = 0
        for line in f:
            line_list = line.split("\t")
            num_bases += 1
            if float(int(line_list[4]) / int(line_list[5])) > 0.5:
                num_meth += 1
        f_out.write(f.name + "\t" + str(num_bases) + "\t" + str(1 - float(num_meth / num_bases)) + "\n")

    f_out.close()

if __name__ == "__main__":
    main()
