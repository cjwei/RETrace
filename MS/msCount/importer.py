def import_targetDict(probe_file):
    '''Import probe file into targetDict'''
    targetDict = {}
    with open(probe_file) as f:
        for line in f:
            target_info = line.split(',')[0]
            info_list = target_info.split('_')
            #Reformat target_id and save to targetDict
            chrom = "_".join(info_list[0:len(info_list)-3])
            chromStart = info_list[-3]
            chromEnd = info_list[-2]
            target_id = chrom + ":" + chromStart + "-" + chromEnd

            #Save relevant information to targetDict
            targetDict[target_id] = {}
            (targetDict[target_id]["chrom"], targetDict[target_id]["chromStart"], targetDict[target_id]["chromEnd"]) = (chrom, chromStart, chromEnd)
            (targetDict[target_id]["num_sub"], targetDict[target_id]["sub_seq"]) = info_list[-1].split('x')
            #Also save the expected up/down_seq around the microsatellite
            MS_frag = line.split()[1]
            targetDict[target_id]["MS_frag"] = MS_frag
            frag_seq = targetDict[target_id]["sub_seq"]*int(targetDict[target_id]["num_sub"])
            (targetDict[target_id]["up_seq"], targetDict[target_id]["down_seq"]) = MS_frag.split(frag_seq)

            targetDict[target_id]["sample_msCount"] = {} #Place holder for sample msCounts in targetDict
    return targetDict

def import_sampleDict(info_file):
    sampleDict = {}
    with open(info_file) as f:
        for line in f:
            if len(line.split()) == 4: #If clone is specified
                (bam, sample, sex, clone) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
                sampleDict[sample]["sex"] = sex
                sampleDict[sample]["clone"] = clone
            elif len(line.split()) == 3: #If clone is not specified
                (bam, sample, sex) = line.split()
                sampleDict[sample] = {}
                sampleDict[sample]["bam"] = bam
                sampleDict[sample]["sex"] = sex
            else:
                print("Incorrect formatting for line (bam, sample, sex, [optional] clone):\n" + "\t".join(line.split()))
    return sampleDict
