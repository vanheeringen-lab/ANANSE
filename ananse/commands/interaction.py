#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function
import sys
import os

import ananse.interaction
import ananse.config as cfg


def interaction(args):
    config = cfg.MotifConfig()
    params = config.get_default_params()

    if not os.path.exists(args.fin_rpkm):
        print("File %s does not exist!" % args.fin_rpkm)
        sys.exit(1)

    params = {
        # "pfmfile": args.pfmfile,
        # "fin_rpkm": args.fin_rpkm,
        # "outfile": args.outfile,
        # "genome": args.genome,
        # "gene_bed": args.annotation,
        # "fpomoter": args.fpomoter,
        # "detail": args.detail,
        "fin_rpkm": args.fin_rpkm,
        "pfmfile": args.pfmfile,
        "fin_expression": args.fin_expression,
        "outfile": args.outfile,
        "genome": args.genome,
        "gene_bed": args.annotation,
        "corrfiles": args.corrfiles,
        "binding": args.binding,
    }

    b = ananse.interaction.Interaction(
        genome=args.genome, gene_bed=args.annotation, pfmfile=args.pfmfile
    )
    b.run_interaction(
        args.fin_rpkm, args.binding, args.fin_expression, args.corrfiles, args.outfile
    )
