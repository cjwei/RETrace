#Because Picard HsMetrics requires a bait interval list, we want to make a dummy interval list based off of the target.info.txt file.  Because we are using the target locationas as our bait locations, we are artificially changing our bait capture efficiency to 1, so we'll have to ignore any statistics that have to deal with bait information.  This is fine because we only care about our target capture efficiency.

#Within GenomeDB, we first have to create a sequence dictionary to be used for the BedToIntervalList command.  Below is the command we ran within folder containing hg19_reference.fa:
#java -jar ~/software/picard.jar CreateSequenceDictionary R=../raw_fasta/hg19_reference.fa O=hg19_reference.dict

for f in ./*.info.txt
do
	fbname=$(basename "$f" .txt)
	./makeBed.py --probes $f --prefix $fbname
done

for g in ./*.info.Picard.bed
do
	gbname=$(basename "$g" .bed)
	java -jar ~/software/picard.jar BedToIntervalList I=$g O=$gbname.interval_list SD=/media/Scratch_SSD/cjwei/GenomeDB/hg19/picard_seqDict/hg19_reference.dict
done
