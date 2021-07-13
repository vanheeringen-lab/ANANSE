#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function

import os
import ananse.network
from ananse.utils import check_path
from dask.distributed import Client, LocalCluster


def network(args):
    ncore = args.ncore
    if ncore is None:
        ncore = min(os.cpu_count(), 4)
    ncore = int(ncore)

    memory_limit = "12GB"

    # With one core more memory is needed
    if ncore == 1:
        memory_limit = "20GB"

    b = ananse.network.Network(
        genome=args.genome,  # checked in CLI
        gene_bed=check_path(args.annotation),
        include_promoter=args.include_promoter,
        include_enhancer=args.include_enhancer
        # pfmfile=args.pfmfile,
        # promoter=args.promoter
    )

    cluster = LocalCluster(
        local_directory=os.environ.get("TMP", None),
        scheduler_port=0,
        dashboard_address=None,  # noqa
        n_workers=ncore,
        threads_per_worker=2,
        memory_limit=memory_limit,
    )
    client = Client(cluster)
    b.run_network(
        binding=check_path(args.binding),
        fin_expression=check_path(args.fin_expression),
        outfile=check_path(args.outfile, error_missing=False),
    )
    client.close()
