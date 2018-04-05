import collections
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plotCounts(plot_file, targetDict):
    plots = PdfPages(plot_file)
    for targetID in sorted(targetDict.keys()):
        if "count_list" in targetDict[targetID].keys():
            count_list = targetDict[targetID]["count_list"]
            count_freq = collections.Counter(count_list)
            count_freq = {float(count): float(num)/float(len(count_list)) for count,num in count_freq.items()}
            x_values = count_freq.keys()
            y_values = count_freq.values()
        
            axes = plt.gca()
            axes.set_ylim([0,1])
            axes.set_xlim([min(count_list)-5, max(count_list)+5])
            plt.vlines(x_values,[0],y_values,linestyle="dashed",color='b')
    
            plt.title(targetID + " (Subunit: " + targetDict[targetID]["sub_seq"] + ", Num Reads: " + str(len(targetDict[targetID]["count_list"])) + ")")
            plt.savefig(plots, format='pdf')
            plt.clf()
    plots.close()
    return