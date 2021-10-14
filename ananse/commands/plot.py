import os

from loguru import logger

from ananse.plot import plot_influence
from ananse.plot import plot_TF_GRN


@logger.catch
def plot(args):
    influence_plot = os.path.join(args.outdir, "influence.pdf")
    plot_influence(args.infile, influence_plot)
    if args.GRN_file:
        influence_grn = os.path.join(args.outdir, "topTF_GRN.png")
        plot_TF_GRN(
            args.infile,
            args.GRN_file,
            influence_grn,
            args.edge_info,
            args.network_algorithm,
            args.n_TFs,
            args.cmap,
            args.edge_min,
            args.full_output,
        )
