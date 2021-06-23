#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function


import ananse.network
from dask.distributed import Client


def network(args):
    # if not os.path.exists(args.fin_rpkm):
    #     print("File %s does not exist!" % args.fin_rpkm)
    #     sys.exit(1)

    ncore = args.ncore
    if ncore is None:
        ncore = 4
    ncore = int(ncore)
    #
    #    memory_limit = args.memory_limit
    #    if args.memory_limit is None:
    #        memory_limit = '24GB'
    #
    memory_limit = "12GB"
    
    with Client(
        n_workers=ncore, threads_per_worker=1, memory_limit=memory_limit,
    ) as client:

        b = ananse.network.Network(
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
            outfile=args.outfile,
        )
