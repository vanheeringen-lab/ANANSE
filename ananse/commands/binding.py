#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function
import sys
import os

import ananse.binding
import ananse.config as cfg


def binding(args):
    config = cfg.MotifConfig()
    params = config.get_default_params()

    if not os.path.exists(args.fin_rpkm):
        print("File %s does not exist!" % args.fin_rpkm)
        sys.exit(1)

    params = {
        "pfmfile": args.pfmfile,
        "fin_rpkm": args.fin_rpkm,
        "outfile": args.outfile,
        "genome": args.genome,
        "gene_bed": args.annotation,
        "fpomoter": args.fpomoter,
        "detail": args.detail,
    }

    a = ananse.binding.Binding(
        genome=args.genome, gene_bed=args.annotation, pfmfile=args.pfmfile
    )
    a.run_binding(args.fin_rpkm, args.outfile)
