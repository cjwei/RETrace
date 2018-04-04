def parseProbes(probe_file):
    '''Import probe file into targetDict'''
    targetDict = {}
    with open(probe_file) as f:
        for line in f:
            target_info = line.split(',')[0]
            info_list = target_info.split('_')
            #Reformat targetID and save to targetDict
            chrom = "_".join(info_list[0:len(info_list)-3])
            chromStart = info_list[-3]
            chromEnd = info_list[-2]
            targetID = chrom + ":" + chromStart + "-" + chromEnd
            
            #Save relevant information to targetDict
            targetDict[targetID] = {}
            (targetDict[targetID]["chrom"], targetDict[targetID]["chromStart"], targetDict[targetID]["chromEnd"]) = (chrom, chromStart, chromEnd)
            (targetDict[targetID]["num_sub"], targetDict[targetID]["sub_seq"]) = info_list[-1].split('x')
            #Also save the expected up/down_seq around the microsatellite
            MS_frag = line.split()[1]
            targetDict[targetID]["MS_frag"] = MS_frag
            frag_seq = targetDict[targetID]["sub_seq"]*int(targetDict[targetID]["num_sub"])
            (targetDict[targetID]["up_seq"], targetDict[targetID]["down_seq"]) = MS_frag.split(frag_seq)
    return targetDict