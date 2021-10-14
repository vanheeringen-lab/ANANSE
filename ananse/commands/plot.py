#!/usr/bin/env python
# Copyright (c) 2021 Simon van Heeringen
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from ananse.plot import plot_influence
from ananse.plot import plot_TF_GRN


def plot(args):
    print(f"args.GRNfile = {args.GRN_file}")
    influence_plot = args.outdir + "/influence.pdf"
    plot_influence(args.infile, influence_plot)
    if args.GRN_file:
        influence_grn = args.outdir + "/topTF_GRN.png"
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
