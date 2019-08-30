#!/usr/bin/env python3
import argparse

def main():
    '''
    We will take as input the PD stats output for statistics of single cell coverage against a given cell line.  The format of input is as follows:
        Sample  Cell Type       Pairwise Dissimilarity  Num CpG Shared  Num Region Shared
    We will then output into the desired format: (classify cells as desired to series name and leave the rest as "All Cells")
        series  x-values    y-values
    '''
    parser = argparse.ArgumentParser(description="Format PD stats output to be compatible with plot_scatter.py")
    parser.add_argument('--input', action="store", dest="input", help="Specify PD.txt stats input used for reformatting")
    parser.add_argument('--series_name', action="store", dest="series_name", help="Series name for plotting legend")
    parser.add_argument('--series_samples', action="store", nargs='+', dest="series_samples", help="Samples/SC belonging to series")
    parser.add_argument('--x_col', action="store", dest="x_col", help="Specify x_column")
    parser.add_argument('--x_name', action="store", dest="x_name", help="Specify x axis name")
    parser.add_argument('--y_name', action="store", dest="y_name", help="Specify y axis name")
    parser.add_argument('--y_col', action="store", dest="y_col", help="Specify y_column")
    parser.add_argument('--output', action="store", dest="output", help="Output filename")
    args = parser.parse_args()

    f_output = open(args.output, 'w')
    f_output.write("Series\t" + args.x_name + "\t" + args.y_name + "\n")
    with open(args.input, 'r') as f_input:
        for line in f_input:
            x_value = line.split()[int(args.x_col)]
            y_value = line.split()[int(args.y_col)]
            if line.split()[0] in args.series_samples:
                f_output.write(args.series_name + "\t" + x_value + "\t" + y_value + "\n")
            else:
                f_output.write("All_Samples\t" + x_value + "\t" + y_value + "\n")
    f_output.close()

if __name__ == "__main__":
    main()
