def parseProbes(probe_file):
    '''Import probe file into probeDict'''
    probeDict = {}
    with open(probe_file) as f:
        for line in f:
            target_info = line.split(',')[0]
            info_list = target_info.split('_')
            #Reformat targetID and save to probeDict
            chrom = "_".join(info_list[0:len(info_list)-3])
            chromStart = info_list[-3]
            chromEnd = info_list[-2]
            targetID = chrom + ":" + chromStart + "-" + chromEnd
            
            #Save relevant information to probeDict
            probeDict[targetID] = {}
            (probeDict[targetID]["num_sub"], probeDict[targetID]["sub_seq"]) = info_list[-1].split('x')
            #Also save the expected up/down_seq around the microsatellite
            MS_frag = line.split()[1]
            frag_seq = probeDict[targetID]["sub_seq"]*int(probeDict[targetID]["num_sub"])
            (probeDict[targetID]["up_seq"], probeDict[targetID]["down_seq"]) = MS_frag.split(frag_seq)
    return probeDict