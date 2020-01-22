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
        Gbf=args.Gbf,
        Gaf=args.Gaf,
        outfile=args.outfile,
        expression=args.expression,
        edges=args.edges,
    )
    a.run_influence(args.plot, args.fin_expression)
