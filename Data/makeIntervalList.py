#!/usr/bin/env python3
import argparse

'''
Usage: python script.py --probes probes.info.txt --prefix
This script will take as input:
    1) probes.info.txt = contains all the probes designed per targetID
The script will then convert the probes.info.txt file into two bed files containing bait (probe) and target locations.  These will be used for Picard CollectHsMetrics.
Below is the format of the interval lists (based on <https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists>):

Picard-style interval files have a SAM-like header that includes a sequence dictionary. The intervals are given in the form <chr> <start> <stop> + <target_name>, with fields separated by tabs, and the coordinates are 1-based (first position in the genome is position 1, not position 0).

@HD     VN:1.0  SO:coordinate
@SQ     SN:1    LN:249250621    AS:GRCh37       UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   M5:1b22b98cdeb4a9304cb5d48026a85128     SP:Homo Sapiens
@SQ     SN:2    LN:243199373    AS:GRCh37       UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   M5:a0d9851da00400dec1098a9255ac712e     SP:Homo Sapiens
1       30366   30503   +       target_1
1       69089   70010   +       target_2
1       367657  368599  +       target_3
1       621094  622036  +       target_4
1       861320  861395  +       target_5
1       865533  865718  +       target_6
This is the preferred format because the explicit sequence dictionary safeguards against accidental misuse (e.g. apply hg18 intervals to an hg19 BAM file). Note that this file is 1-based, not 0-based (the first position in the genome is position 1).

'''

#%%
def makeIntervalList():
    parser = argparse.ArgumentParser(description="Convert probe file to bed format")
    parser.add_argument('--probes', action="store", dest="probe_file", help="Info txt file containing hybridization probe targets")
    parser.add_argument('--prefix', action="store", dest="prefix")
    args = parser.parse_args()

    output_target = open(args.prefix + ".targets.interval_list","w")
    output_probe = open(args.prefix + ".probes.interval_list","w")

    with open(args.probe_file) as f:
        for line in f:
            target_info = line.split(',')[0]
            info_list = target_info.split('_')
            #Obtain chromosomal position
            chrom = "_".join(info_list[0:len(info_list)-3])
            target_chromStart = int(info_list[-3])
            target_chromEnd = int(info_list[-2])

            #Obtain subunit information and reformat targetID
            (num_sub,sub_seq) = info_list[-1].split('x')
            targetID = chrom + ":" + str(target_chromStart) + "-" + str(target_chromEnd) + "_" + info_list[-1]
            probeID = "Probe_" + targetID

            #We next want to determine the chromosomal positions of the probe itself (not just the target).  This is done by looking at the length of flanking sequences around the probe
            probe_seq = line.split()[1]
            ms_seq = sub_seq * int(num_sub)
            #There are a few cases (i.e. chr6:55179847-55179875) in which microsatellite frag_seq is found more than once in reference fragment.  Thus, we need to choose the optimal up/down-seq based on which pair has the greater length
            frag_split = probe_seq.split(ms_seq)
            (up_flank, down_flank) = ('', '')
            for i in range(len(frag_split) - 1):
                if len(frag_split[i]) > len(up_flank) and len(frag_split[i + 1]) > len(down_flank):
                    (up_flank, down_flank) = frag_split[i:i+2]
            probe_chromStart = target_chromStart - len(up_flank)
            probe_chromEnd = target_chromEnd - len(down_flank)

            #Print relevant information in bed format
            output_target.write(chrom + "\t" + str(target_chromStart) + "\t" + str(target_chromEnd) + "\t+\t" + targetID + "\n")
            output_probe.write(chrom + "\t" + str(probe_chromStart) + "\t" + str(probe_chromEnd) + "\t+\t" + probeID + "\n")

    output_target.close()
    output_probe.close()

#%%
if __name__ == "__main__":
    makeIntervalList()
