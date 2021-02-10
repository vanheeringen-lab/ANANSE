import os

import numpy as np
import pandas as pd
import seaborn as sns


# https://stackoverflow.com/questions/18534562/scipy-lognormal-fitting
# https://stackoverflow.com/questions/41940726/scipy-lognorm-fitting-to-histogram
# https://stackoverflow.com/questions/26406056/a-lognormal-distribution-in-python
# https://stackoverflow.com/questions/15630647/fitting-lognormal-distribution-using-scipy-vs-matlab

def distplot(infile, score_col="scaled_score", show=False):
    """
    generate simple distplot from bedfile
    """
    bed = pd.read_csv(infile, header=0, sep="\t")
    scores = pd.Series(bed[score_col])
    bins = min(30, len(scores))  # too many bins = bad
    fig = sns.histplot(scores, kde=True, stat="density", bins=bins, alpha=0.2)
    fig.set_yscale('log')  # most methods are log scaled

    # # exclude outliers from plot
    # y_min = np.percentile(scores, 1)
    # y_max = np.percentile(scores, 99)
    # fig.axes.set_ylim([y_min, y_max])

    title = os.path.splitext(os.path.basename(infile))[0].replace(".out", "")
    fig.set_title(f"{title} score distribution")
    fig.xaxis.set_label_text("Score")

    if show:
        fig.figure.show()
    else:
        outfile = infile.replace(".bed", ".png")
        fig.figure.savefig(outfile, orientation='landscape')


# distplot("../tests/output/ScorePeaks_scale.out.bed", show=True)
# distplot("../tests/output/ScorePeaks_logscale.out.bed", show=True)
# distplot("../tests/output/ScorePeaks_lognorm.out.bed", show=True)
# distplot("../tests/output/ScorePeaks_loglaplace.out.bed", show=True)
# distplot("../tests/output/ScorePeaks_peakrank.out.bed", show=True)
# distplot("../tests/output/ScorePeaks_peakrankfile.out.bed", show=True)


# from matplotlib import pyplot as plt
#
# fig, (ax1) = plt.subplots(1, 1)
# scores.plot(kind='hist', density=True, bins=bins, alpha=0.5, ax=ax1)
# ax1.set_title(f"{title} score distribution")
# ax1.set_xlabel('Score')
# fig.show()


# from matplotlib import pyplot as plt
# from matplotlib.ticker import FormatStrFormatter
#
# fig, (ax1, ax2) = plt.subplots(1, 2)
# fig.suptitle(f"{title} score distribution")
#
# sns.histplot(scores, ax=ax1, kde=True, stat="density")
# ax1.set_title("raw score")
# ax1.xaxis.set_label_text("Score")
# ax1.xaxis.set_major_locator(plt.MaxNLocator(6))
# ax1.xaxis.set_major_formatter(FormatStrFormatter('%i'))
#
# sns.histplot(np.log10(scores+1), ax=ax2, kde=True, stat="density")  # log_scale=10
# ax2.set_title(f"log10 score")
# ax2.xaxis.set_label_text("Score")
#
# fig.xaxis.set_major_locator(plt.MaxNLocator(10))
# fig.xaxis.set_tick_params(rotation=15)  # for long floats
# fig.xaxis.set_major_formatter(FormatStrFormatter('%i'))  # integer. for floats, use '%.3f'
# fig.set_size_inches(10, 5)
# fig.savefig(outfile, orientation='landscape')
