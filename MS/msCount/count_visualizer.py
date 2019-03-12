from tqdm import tqdm
import collections
from itertools import cycle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

def plotCounts(prefix, targetDict):
    #Plot MS counts
    plots = PdfPages(prefix + '.msCount.pdf')
    color_options = ['xkcd:pale green','xkcd:pale blue','xkcd:light grey','xkcd:pale pink']
    color_cycle = cycle(color_options)
    next_color = next(color_cycle)
    for target_id in tqdm(sorted(targetDict.keys())):
        if len(targetDict[target_id]["sample_msCount"].keys()) >= 2:
            select_color, next_color = next_color, next(color_cycle)
        else:
            select_color = 'xkcd:white'
        for sample in sorted(targetDict[target_id]["sample_msCount"].keys()):
            #Plot raw msCounts (blue)
            ax = plt.gca()
            ax.set_ylim([0,1])
            ax.set_xlim([min(targetDict[target_id]["sample_msCount"][sample]["count_list"]) - 5, max(targetDict[target_id]["sample_msCount"][sample]["count_list"]) + 5])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.vlines(targetDict[target_id]["sample_msCount"][sample]["count_list"], [0], targetDict[target_id]["sample_msCount"][sample]["count_freq"], linestyle="dashed", color="b")
            plt.title(target_id + ", " + sample)
            ax.set_facecolor(select_color)
            plt.savefig(plots, format='pdf')
            plt.clf()
    plots.close()
    return
