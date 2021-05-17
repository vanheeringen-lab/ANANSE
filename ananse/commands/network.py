#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function


import ananse.network


def network(args):
    # if not os.path.exists(args.fin_rpkm):
    #     print("File %s does not exist!" % args.fin_rpkm)
    #     sys.exit(1)

    b = ananse.network.Network(
        ncore=args.ncore,
        genome=args.genome,
        gene_bed=args.annotation,
        include_promoter=args.include_promoter,
        include_enhancer=args.include_enhancer
        # pfmfile=args.pfmfile,
        # promoter=args.promoter
    )
    b.run_network(
        binding=args.binding,
        fin_expression=args.fin_expression,
        corrfiles=args.corrfiles,
        outfile=args.outfile,
    )
