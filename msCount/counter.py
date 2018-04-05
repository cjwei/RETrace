def simple_counter(read, targetDict, targetID):
    '''This is a simple counter to determine microsatellite subunits with exact match of up_seq and down_seq'''
    count = -1
    if targetDict[targetID]["up_seq"] in read and targetDict[targetID]["down_seq"] in read:    
        up_start = read.index(targetDict[targetID]["up_seq"])
        down_start = read.index(targetDict[targetID]["down_seq"])
        count = (down_start-(up_start+len(targetDict[targetID]["up_seq"])))/len(targetDict[targetID]["sub_seq"])
    return count

def counter(count_type, read_list, targetDict, targetID):
    count_list = []
    for read in read_list:
        if count_type=="simple":
            count_list.append(simple_counter(read, targetDict, targetID))
    count_list = [x for x in count_list if x>=0] #Filter out all counts that are <0, which indicates no up/down-seq 
    return count_list