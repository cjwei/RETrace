#!/usr/bin/env python3
import argparse
import RETrace.MS
import RETrace.Methyl
import sys

def parse_args():
    # create the top-level parser
    parser = argparse.ArgumentParser(description="Running RETrace")

    subparsers = parser.add_subparsers(title="functions", dest="command", metavar="")

    '''
    MS subparsers
    '''
    add_HipSTR_allelotype_subparser(subparsers)
    add_Custom_allelotype_subparser(subparsers)
    add_merge_allelotype_subparser(subparsers)
    add_buildPhylo_subparser(subparsers)
    add_iterPhylo_subparser(subparsers)
    add_evalPhylo_subparser(subparsers)
    add_viewPhylo_subparser(subparsers)

    '''
    Methyl subparsers
    '''
    add_importMethyl_subparser(subparsers)
    add_refPD_subparser(subparsers)
    add_methRate_subparser(subparsers)

    '''
    MS+Methyl Combined subparsers
    '''
    add_combinedPhylo_subparser(subparsers)

    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-h"])
        exit()

    if args.command == "HipSTR_allelotype":
        from RETrace.MS.HipSTR import HipSTR_allelotype
        HipSTR_allelotype(args.sample_info, args.HipSTR_vcf, args.alleleDict_file,
            args.fasta_loc, args.picard_loc, args.HipSTR_loc,
            args.target_bed, args.target_info,
            args.min_qual, args.min_reads, args.max_stutter)

    elif args.command == "Custom_allelotype":
        from RETrace.MS.Custom import Custom_allelotype
        Custom_allelotype(args.sample_info, args.prefix,
            args.target_info, args.nproc, args.min_cov, args.min_ratio)

    elif args.command == "merge_allelotype":
        from RETrace.MS.merge import merge_allelotype
        merge_allelotype(args.file, args.output)

    elif args.command == "buildPhylo":
        from RETrace.MS.buildPhylo import buildPhylo
        buildPhylo(args.sample_info, args.sample_list, args.prefix, args.target_info,
            args.alleleDict_file, args.dist_metric, args.outgroup, args.bootstrap, args.merge)

    elif args.command == "iterPhylo":
        from RETrace.MS.iterPhylo import iterPhylo
        iterPhylo(args.sample_info, args.sample_list, args.prefix, args.target_info,
            args.alleleDict_file, args.dist_metric)

    elif args.command == "evalPhylo":
        from RETrace.MS.evalPhylo import evalPhylo
        evalPhylo(args.sample_info, args.prefix, args.exVivo_rootDist, args.tree_file, args.nproc, args.distDict_file)

    elif args.command == "viewPhylo":
        from RETrace.MS.viewPhylo import viewPhylo
        viewPhylo(args.sample_info, args.tree_file, args.prefix, args.bootstrap)

    elif args.command == "importMethyl":
        from RETrace.Methyl.importMethyl import importMethyl
        importMethyl(args.sample_info, args.prefix, args.ref_CGI)

    elif args.command == "refPD":
        from RETrace.Methyl.refPD import refPD
        refPD(args.sample_list, args.ref_info, args.sample_methDict, args.prefix, args.DMR, args.min_shared, args.min_rate)

    elif args.command == "pairwise_methRate":
        from RETrace.Methyl.pairwise_methRate import pairwise_methRate
        pairwise_methRate(args.sample_list, args.sample_methDict, args.Ensemble_gff, args.prefix)

    elif args.command == "combinedPhylo":
        from RETrace.Combined.combinedPhylo import combinedPhylo
        combinedPhylo(args.sample_list, args.target_info, args.alleleDict_file, args.methDict_file, args.Ensemble_gff, args.prefix,
            args.dist_metric, args.ratio, args.outgroup)

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
    parser_HipSTR_allelotype_req.add_argument("--HipSTR_vcf",
        action="store",
        dest="HipSTR_vcf",
        help="Name of HipSTR vcf output file")
    parser_HipSTR_allelotype_req.add_argument("--alleleDict",
        action="store",
        dest="alleleDict_file",
        help="Pickle file containing alleleDict calculated from HipSTR")
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
        help="Specify minimum percentage of reads supporting the resulting allelotype")

def add_merge_allelotype_subparser(subparsers):
    # create the parser for "merge_allelotype" command
    parser_merge_allelotype = subparsers.add_parser(
        "merge_allelotype",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Merge multiple alleleDict pickle files from either Custom or HipSTR allelotyping")

    parser_merge_allelotype_req = parser_merge_allelotype.add_argument_group("required inputs")
    parser_merge_allelotype_req.add_argument("file",
        nargs = "+",
        help = "Input pickle alleleDict files")
    parser_merge_allelotype_req.add_argument("--output",
        dest = "output",
        help = "Specify output pickle file name")

def add_buildPhylo_subparser(subparsers):
    # create the parser for "buildPhylo" comand
    parser_buildPhylo = subparsers.add_parser(
        "buildPhylo",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Build phylogenetic tree given allelotype of single cells")

    parser_buildPhylo_req = parser_buildPhylo.add_argument_group("required inputs")
    parser_buildPhylo_req.add_argument("--sample_list",
        action="store",
        dest="sample_list",
        help="File containing sample names to be included in calculations (one sample per line)")
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
    parser_buildPhylo_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone/cluster)")

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
    parser_buildPhylo_opt.add_argument("--merge",
        action="store",
        default=1,
        type=int,
        help="Indicate the number of single cells to merge together (based on clone/cluster as specified in sample_info file)")

def add_iterPhylo_subparser(subparsers):
    # create the parser for "iterPhylo" comand
    parser_iterPhylo = subparsers.add_parser(
        "iterPhylo",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Iteratively build phylogenetic tree given allelotype of single cells")

    parser_iterPhylo_req = parser_iterPhylo.add_argument_group("required inputs")
    parser_iterPhylo_req.add_argument("--sample_list",
        action="store",
        dest="sample_list",
        help="File containing sample names to be included in calculations (one sample per line)")
    parser_iterPhylo_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Output prefix for MS phylogenetic tree")
    parser_iterPhylo_req.add_argument("--target_info",
        action="store",
        dest="target_info",
        default="~/software/RETrace/Data/CA_order.20171211-20190301.info.txt",
        help="Location of probe info file")
    parser_iterPhylo_req.add_argument("--alleleDict",
        action="store",
        dest="alleleDict_file",
        help="Pickle file containing alleleDict calculated from either HipSTR_allelotype or Custom_allelotype")
    parser_iterPhylo_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone/cluster)")

    parser_iterPhylo_opt = parser_iterPhylo.add_argument_group("optional inputs")
    parser_iterPhylo_opt.add_argument("--dist",
        action="store",
        dest="dist_metric",
        default="EqorNot",
        help="Specify distance metric for pairwise comparisons [Abs, EqorNot, Chi]")

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
    parser_evalPhylo_req.add_argument("--exVivo_rootDist",
        action="store",
        dest="exVivo_rootDist",
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
    parser_evalPhylo_req.add_argument("--distDict",
        action="store",
        dest="distDict_file",
        help="Pickle file containing distDict calculated from buildPhylo command.  If specified, will calculate the expected vs actual number allele diff")

    parser_evalPhylo_opt = parser_evalPhylo.add_argument_group("optional inputs")
    parser_evalPhylo_opt.add_argument("--nproc",
        action="store",
        dest="nproc",
        default=10,
        type=int,
        help="Specify number of processors for evaluating accuracy of tree")

def add_viewPhylo_subparser(subparsers):
    # create the parser for "viewPhylo" command
    parser_viewPhylo = subparsers.add_parser(
        "viewPhylo",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Utilize ete3 to view phylogenetic tree derived from buildPhylo")

    parser_viewPhylo_req = parser_viewPhylo.add_argument_group("required inputs")
    parser_viewPhylo_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing sample information (bam, sample_name, sex, [optional] clone)")
    parser_viewPhylo_req.add_argument("--tree",
        action="store",
        dest="tree_file",
        help="Newick tree file")
    parser_viewPhylo_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Output prefix for tree file")

    parser_viewPhylo_opt = parser_viewPhylo.add_argument_group("required inputs")
    parser_viewPhylo_opt.add_argument("-bootstrap",
        action="store_true",
        default=False,
        help="Flag for indicating whether we want to bootstrap the tree to determine node support (random sampling across all samples)")

def add_importMethyl_subparser(subparsers):
    # create the parser for "importMethyl" command
    parser_importMethyl = subparsers.add_parser(
        "importMethyl",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Import methylation calls derived from methylpy TSV files for samples/cell type ref (make sure to run methylpy prior to import)")

    parser_importMethyl_req = parser_importMethyl.add_argument_group("required inputs")
    parser_importMethyl_req.add_argument("--sample_info",
        action="store",
        dest="sample_info",
        help="Tab-delimited file containing the sample/cell type tsv files/locations along with sample/cell type names")
    parser_importMethyl_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Prefix name for exporting methDict and coverage statistics information")
    parser_importMethyl_req.add_argument("--ref_CGI",
        action="store",
        dest="ref_CGI",
        help="Bed file containing reference genome CGI locations")

def add_refPD_subparser(subparsers):
    # create the parser for "refPD" command
    parser_refPD = subparsers.add_parser(
        "refPD",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Calculate and plot pairwise dissimilarity between single cells and reference cell types")

    parser_refPD_req = parser_refPD.add_argument_group("required inputs")
    parser_refPD_req.add_argument("--sample_list",
        action="store",
        dest="sample_list",
        help="File containing sample names to be included in PD calculations (one sample per line)")
    parser_refPD_req.add_argument("--ref_info",
        action="store",
        dest="ref_info",
        help="Tab-delmited file containing the cell type tsv files/locations along with cell type names")
    parser_refPD_req.add_argument("--sample_methDict",
        action="store",
        dest="sample_methDict",
        help="Pre-computed methDict pickle file containing sample information and methylation calls")
    parser_refPD_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Prefix name for exporting calPD output")

    parser_refPD_opt = parser_refPD.add_argument_group("optional inputs")
    parser_refPD_opt.add_argument("-DMR",
        action="store_true",
        dest="DMR",
        default=False,
        help="Flag to constrain bases of interest in refernce file to DMR regions")
    parser_refPD_opt.add_argument("--min_shared",
        action="store",
        dest="min_shared",
        default=100,
        type=int,
        help="Minimum number of CpG sites shared between each sample and cell type comparison")
    parser_refPD_opt.add_argument("--min_rate",
        action="store",
        dest="min_rate",
        default=0.5,
        type=float,
        help="Minimum methylated or unmethylated rate for cellType reference for filtering confident methyl calls (0.5 = no filtering; 1 = only perfect call)")

def add_methRate_subparser(subparsers):
    # create the parser for "pairwise_methRate" command
    parser_methRate = subparsers.add_parser(
        "pairwise_methRate",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Calculate methRate across Ensembl Regulatory Build windows")

    parser_methRate_req = parser_methRate.add_argument_group("required inputs")
    parser_methRate_req.add_argument("--sample_list",
        action="store",
        dest="sample_list",
        help="File containing sample names to be included in pairwise methRate calculations (one sample per line)")
    parser_methRate_req.add_argument("--sample_methDict",
        action="store",
        dest="sample_methDict",
        help="Pre-computed methDict pickle file containing sample information and methylation calls")
    parser_methRate_req.add_argument("--gff",
        action="store",
        dest="Ensemble_gff",
        help="Gzipped gff file containing Ensembl Regulatory Build windows")
    parser_methRate_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Prefix name for exporting pairwise methRate csv output")

def add_combinedPhylo_subparser(subparsers):
    # create the parser for "combinedPhylo" command
    parser_combinedPhylo = subparsers.add_parser(
        "combinedPhylo",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "Build phylogenty by combining MS and Methyl distances between single cell pairs")

    parser_combinedPhylo_req = parser_combinedPhylo.add_argument_group("required inputs")
    parser_combinedPhylo_req.add_argument("--sample_list",
        action="store",
        dest="sample_list",
        help="File containing sample names to be included in combinedPhylo calculations (one sample per line)")
    parser_combinedPhylo_req.add_argument("--target_info",
        action="store",
        dest="target_info",
        default="~/software/RETrace/Data/CA_order.20171211-20190301.info.txt",
        help="Location of probe info file")
    parser_combinedPhylo_req.add_argument("--alleleDict",
        action="store",
        dest="alleleDict_file",
        help="Pickle file containing alleleDict calculated from either HipSTR_allelotype or Custom_allelotype")
    parser_combinedPhylo_req.add_argument("--methDict",
        action="store",
        dest="methDict_file",
        help="Pre-computed methDict pickle file containing sample information and methylation calls")
    parser_combinedPhylo_req.add_argument("--gff",
        action="store",
        dest="Ensemble_gff",
        help="Gzipped gff file containing Ensembl Regulatory Build windows")
    parser_combinedPhylo_req.add_argument("--prefix",
        action="store",
        dest="prefix",
        help="Prefix name for exporting pairwise methRate csv output")

    parser_combinedPhylo_opt = parser_combinedPhylo.add_argument_group("optional inputs")
    parser_combinedPhylo_opt.add_argument("--dist",
        action="store",
        dest="dist_metric",
        default="EqorNot",
        help="Specify distance metric for pairwise comparisons [Abs, EqorNot]")
    parser_combinedPhylo_opt.add_argument("--ratio",
        action="store",
        dest="ratio",
        default=1.0,
        type=float,
        help="Ratio of MS to Methyl contribution for building phylogeny (0.0 = all methyl, 0.5 = equal MS/Methyl, 1.0 = all MS)")
    parser_combinedPhylo_opt.add_argument("--outgroup",
        action="store",
        dest="outgroup",
        default="Midpoint",
        help="Specify outgroup for rooted NJ tree [NA, Midpoint]")
