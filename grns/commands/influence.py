#!/usr/bin/env python 
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from __future__ import print_function
import sys
import os

import grns.influence 
import grns.config as cfg

def influence(args):
    config = cfg.MotifConfig()
    params = config.get_default_params()

    if not os.path.exists(args.Gbf) and not os.path.exists(args.Gaf):
        print("At list need one network file")
        sys.exit(1)
    
    params = {
        "outfile": args.outfile,
        "expression": args.expression,
        "fin_expression": args.fin_expression,
        "Gbf": args.Gbf,
        "Gaf": args.Gaf,
        "plot": args.plot,
    }

    a = grns.influence.Influence(Gbf = args.Gbf, Gaf = args.Gaf, outfile= args.outfile, expression=args.expression, edges=10000)    
    a.run_influence(args.plot, args.fin_expression)
