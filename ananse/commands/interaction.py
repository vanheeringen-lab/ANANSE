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


def interaction(args):
    if not os.path.exists(args.fin_rpkm):
        print("File %s does not exist!" % args.fin_rpkm)
        sys.exit(1)

    b = ananse.interaction.Interaction(
        genome=args.genome, gene_bed=args.annotation, pfmfile=args.pfmfile
    )
    b.run_interaction(
        args.fin_rpkm, args.binding, args.fin_expression, args.corrfiles, args.outfile
    )
