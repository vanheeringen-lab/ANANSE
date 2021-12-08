import os

from loguru import logger

from ananse.plot import plot_influence, plot_TF_GRN


@logger.catch
def plot(args):
    os.makedirs(args.outdir, exist_ok=True)
    influence_plot = os.path.join(args.outdir, f"influence.{args.ftype}")
    plot_influence(args.infile, influence_plot)
    if args.GRN_file:
        influence_grn = os.path.join(args.outdir, f"topTF_GRN.{args.ftype}")
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
