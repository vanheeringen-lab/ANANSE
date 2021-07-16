#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function

from loguru import logger

import ananse.influence
from ananse.utils import check_path


@logger.catch
def influence(args):
    a = ananse.influence.Influence(
        ncore=args.ncore,  # --ncore (optional)
        Gbf=check_path(args.Gbf),  # --source (Gbf = GRN before)
        Gaf=check_path(args.Gaf),  # --target (Gaf = GRN after)
        outfile=check_path(args.outfile, error_missing=False),  # --output
        degenes=check_path(
            args.expression
        ),  # --degenes (HGNC gene names, padj and log2foldchanges)
        edges=args.edges,  # --edges (optional)
    )
    a.run_influence(args.plot)  # -p
