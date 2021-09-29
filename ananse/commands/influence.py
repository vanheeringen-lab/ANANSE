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
        GRN_source_file=check_path(args.GRN_source_file),
        GRN_target_file=check_path(args.GRN_target_file),
        outfile=check_path(args.outfile, error_missing=False),
        full_output=args.full_output,
        GRNsort_column=args.GRNsort_column,
        edge_info=args.edge_info,
        padj_cutoff=args.padj_cutoff,
        degenes=check_path(
            args.expression
        ),  # --degenes (HGNC gene names, padj and log2foldchanges)
        edges=args.edges,  # --edges (optional)
    )
    a.run_influence(args.plot)  # -p
