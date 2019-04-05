from tqdm import tqdm
import collections
from itertools import cycle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

def plotCounts(plot_file, targetDict):
    #Plot MS counts
    plots = PdfPages(plot_file)
    color_options = ['xkcd:pale green','xkcd:pale blue','xkcd:light grey','xkcd:pale pink']
    color_cycle = cycle(color_options)
    next_color = next(color_cycle)
    for target_id in tqdm(sorted(targetDict.keys())):
        if len(targetDict[target_id]["sample"].keys()) >= 2:
            select_color, next_color = next_color, next(color_cycle)
        else:
            select_color = 'xkcd:white'
        for sample in sorted(targetDict[target_id]["sample"].keys()):
            #Create count_list and count_freq
            count_list = sorted(set(targetDict[target_id]["sample"][sample]))
            count_freq = [targetDict[target_id]["sample"][sample].count(i)/len(targetDict[target_id]["sample"][sample]) for i in count_list]
            #Plot raw msCounts (blue)
            ax = plt.gca()
            ax.set_ylim([0,1])
            ax.set_xlim([min(count_list) - 5, max(count_list) + 5])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.vlines(count_list, [0], count_freq, linestyle="dashed", color="b")
            plt.title(target_id + ", " + sample)
            ax.set_facecolor(select_color)
            plt.savefig(plots, format='pdf')
            plt.clf()
    plots.close()
    return

def plotAlleles(plot_file, num_plot, targetDict, alleleDict):
    #Plot MS counts and allelotype call
    plots = PdfPages(plot_file)
    color_options = ['xkcd:pale green','xkcd:pale blue','xkcd:light grey','xkcd:pale pink']
    color_cycle = cycle(color_options)
    next_color = next(color_cycle)
    for target_id in tqdm(sorted(targetDict.keys())[0:min(num_plot,len(targetDict.keys()))]): #Plot either a subset of target_ids (num_plot) or all target_id (can be specified by using inf)
        if len(targetDict[target_id]["sample"].keys()) >= 2:
            select_color, next_color = next_color, next(color_cycle)
        else:
            select_color = 'xkcd:white'
        for sample in sorted(targetDict[target_id]["sample"].keys()):
            # print(target_id + "\t" + sample)
            #Create count_list and count_freq
            count_list = sorted(set(targetDict[target_id]["sample"][sample]))
            count_freq = [targetDict[target_id]["sample"][sample].count(i)/len(targetDict[target_id]["sample"][sample]) for i in count_list]
            #Plot raw msCounts (blue)
            ax = plt.gca()
            ax.set_ylim([0,1])
            ax.set_xlim([min(count_list) - 5, max(count_list) + 5])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            if sample in alleleDict[target_id]["sample"].keys(): #We want to plot allelotypes if called for given target_id/sample
                # print(alleleDict[target_id]["sample"][sample])
                plt.plot(alleleDict[target_id]["sample"][sample], [0.5]*len(alleleDict[target_id]["sample"][sample]), 'ro')
            plt.vlines(count_list, [0], count_freq, linestyle="dashed", color="b")
            plt.title(target_id + ", " + sample)
            ax.set_facecolor(select_color)
            plt.savefig(plots, format='pdf')
            plt.clf()
    plots.close()
    return
