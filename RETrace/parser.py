#!/usr/bin/env python3
import argparse
import RETrace.MS
import sys

def parse_args():
    # create the top-level parser
    parser = argparse.ArgumentParser(description="Running RETrace")

    subparsers = parser.add_subparsers(title="functions", dest="command", metavar="")

    add_HipSTR_allelotype_subparser(subparsers)
    add_Custom_allelotype_subparser(subparsers)
    add_buildPhylo_subparser(subparsers)
    add_evalPhylo_subparser(subparsers)

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-h"])
        exit()

    if args.command == "HipSTR_allelotype":
        from RETrace.MS.HipSTR import HipSTR_allelotype
        HipSTR_allelotype(args.sample_info, args.prefix,
            args.fasta_loc, args.picard_loc, args.HipSTR_loc,
            args.target_bed, args.target_info,
            args.min_qual, args.min_reads, args.max_stutter)

    elif args.command == "Custom_allelotype":
        from RETrace.MS.Custom import Custom_allelotype
        Custom_allelotype(args.sample_info, args.prefix,
            args.target_info, args.nproc, args.min_cov, args.min_ratio)

    elif args.command == "buildPhylo":
        from RETrace.MS.buildPhylo import buildPhylo
        buildPhylo(args.sample_info, args.prefix, args.target_info,
            args.alleleDict_file, args.dist_metric, args.outgroup, args.bootstrap)

    elif args.command == "evalPhylo":
        from RETrace.MS.evalPhylo import evalPhylo
        evalPhylo(args.sample_info, args.prefix, args.exVivo_dist, args.tree_file, args.nproc)

def add_HipSTR_allelotype_subparser(subparsers):
    # create the parser for "HipSTR_allelotype" command
    parser_HipSTR_allelotype = subparsers.add_parser(
        "HipSTR_allelotype",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Run HipSTR for microsatellite calling and allelotyping")

    # add options
    parser_HipSTR_allelotype_req = parser_HipSTR_allelotype.add_argument_group("required inputs")
    parser_HipSTR_allelotype_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)")
    parser_HipSTR_allelotype_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Output prefix for HipSTR vcf files along with any pickle dictionaries saving allelotype information")
    parser_HipSTR_allelotype_req.add_argument("--fasta_loc",
        action="store",
        dest="fasta_loc",
        help="Location of hg19 fasta file for running HipSTR")
    parser_HipSTR_allelotype_req.add_argument("--picard_loc",
        action="store",
        dest="picard_loc",
        default="~/software",
        help="Location of picard.jar for adding read groups to bam files")
    parser_HipSTR_allelotype_req.add_argument("--HipSTR_loc",
        action="store",
        dest="HipSTR_loc",
        default="~/software/HipSTR",
        help="Location of HipSTR for allelotype calling")
    parser_HipSTR_allelotype_req.add_argument("--target_bed",
        action="store",
        dest="target_bed",
        default="~/software/RETrace/Data/CA_order.20171211-20190301.info.bed",
        help="Location of probe bed file")
    parser_HipSTR_allelotype_req.add_argument("--target_info",
        action="store",
        dest="target_info",
        default="~/software/RETrace/Data/CA_order.20171211-20190301.info.txt",
        help="Location of probe info file")

    parser_HipSTR_allelotype_opt = parser_HipSTR_allelotype.add_argument_group("optional inputs")
    parser_HipSTR_allelotype_opt.add_argument("--min_qual",
        action="store",
        dest="min_qual",
        default=0.0,
        type=float,
        help="Specify the minimum posterior probability of genotype for filtering")
    parser_HipSTR_allelotype_opt.add_argument("--min_reads",
        action="store",
        dest="min_reads",
        default=10,
        type=int,
        help="Cutoff for minimum number of reads required for calling allelotype")
    parser_HipSTR_allelotype_opt.add_argument("--max_stutter",
        action="store",
        dest="max_stutter",
        default=1.0,
        type=float,
        help="Define maximum number of reads that can be classified as stutter")

def add_Custom_allelotype_subparser(subparsers):
    # create the parser for "Custom_allelotype" command
    parser_Custom_allelotype = subparsers.add_parser(
        "Custom_allelotype",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Run scripts for custom microsatellite calling and allelotyping (based off of LobSTR)")

    parser_Custom_allelotype_req = parser_Custom_allelotype.add_argument_group("required inputs")
    parser_Custom_allelotype_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)")
    parser_Custom_allelotype_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Output prefix for custom msCount/allelotype files along with any pickle dictionaries saving allelotype information")
    parser_Custom_allelotype_req.add_argument("--target_info",
        action="store",
        dest="target_info",
        default="~/software/RETrace/Data/CA_order.20171211-20190301.info.txt",
        help="Location of probe info file")

    parser_Custom_allelotype_opt = parser_Custom_allelotype.add_argument_group("optional inputs")
    parser_Custom_allelotype_opt.add_argument("--nproc",
        action="store",
        dest="nproc",
        default=10,
        type=int,
        help="Specify number of processors used for msCount (~2 days running on 10 processors for 10k targets)")
    parser_Custom_allelotype_opt.add_argument("--min_cov",
        action="store",
        dest="min_cov",
        default=10,
        type=int,
        help="Specify minimum coverage for each target_id within each sample to call allelotype")
    parser_Custom_allelotype_opt.add_argument("--min_ratio",
        action="store",
        dest="min_ratio",
        default=0.2,
        type=float,
        help="Specify minimum percentae of reads supporting the resulting allelotype")

def add_buildPhylo_subparser(subparsers):
    # create the parser for "buildPhylo" comand
    parser_buildPhylo = subparsers.add_parser(
        "buildPhylo",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Build phylogenetic tree given allelotype of single cells")

    parser_buildPhylo_req = parser_buildPhylo.add_argument_group("required inputs")
    parser_buildPhylo_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)")
    parser_buildPhylo_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Output prefix for MS phylogenetic tree")
    parser_buildPhylo_req.add_argument("--target_info",
        action="store",
        dest="target_info",
        default="~/software/RETrace/Data/CA_order.20171211-20190301.info.txt",
        help="Location of probe info file")
    parser_buildPhylo_req.add_argument("--alleleDict",
        action="store",
        dest="alleleDict_file",
        help="Pickle file containing alleleDict calculated from either HipSTR_allelotype or Custom_allelotype")

    parser_buildPhylo_opt = parser_buildPhylo.add_argument_group("optional inputs")
    parser_buildPhylo_opt.add_argument("--dist",
        action="store",
        dest="dist_metric",
        default="EqorNot",
        help="Specify distance metric for pairwise comparisons [Abs, EqorNot]")
    parser_buildPhylo_opt.add_argument("--outgroup",
        action="store",
        dest="outgroup",
        default="Midpoint",
        help="Specify outgroup for rooted NJ tree [NA, Midpoint]")
    parser_buildPhylo_opt.add_argument("-bootstrap",
        action="store_true",
        default=False,
        help="Flag for indicating whether we want to bootstrap the tree to determine node support (random sampling across all samples)")

def add_evalPhylo_subparser(subparsers):
    # create the parser for "buildPhylo" comand
    parser_evalPhylo = subparsers.add_parser(
        "evalPhylo",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Evaluate accuracy of phylogeny given exVivo tree data")

    parser_evalPhylo_req = parser_evalPhylo.add_argument_group("required inputs")
    parser_evalPhylo_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)")
    parser_evalPhylo_req.add_argument("--exVivo",
        action="store",
        dest="exVivo_dist",
        default="~/software/RETrace/Data/exVivo.rootDist.csv",
        help="csv file containing MRCA distance from root, as approximated in units of cell divisions")
    parser_evalPhylo_req.add_argument("--tree",
        action="store",
        dest="tree_file",
        help="Newick tree file")
    parser_evalPhylo_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Output prefix for plot file containing proportion of correct triplets in tree")
    parser_evalPhylo_req.add_argument("--nproc",
        action="store",
        dest="nproc",
        default=10,
        type=int,
        help="Specify number of processors for evaluating accuracy of tree")
