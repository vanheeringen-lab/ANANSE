#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function
import sys
import os

import ananse.influence


def influence(args):
    # if not os.path.exists(args.outfile):
    #     print("File %s does not exist!" % args.outfile)
    #     sys.exit(1)

    a = ananse.influence.Influence(
        ncore=args.ncore,  # --ncore (optional)
        Gbf=args.Gbf,  # --source (Gbf = GRN before)
        Gaf=args.Gaf,  # --target (Gaf = GRN after)
        outfile=args.outfile,  # --output
        degenes=args.expression,  # --degenes (HGNC gene names, padj and log2foldchanges)
        edges=args.edges,  # --edges (optional)
    )
    a.run_influence(
        args.plot
    )  # -p 
