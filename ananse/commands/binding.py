#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import sys
import os
import pandas as pd
from ananse.peakpredictor import predict_peaks


def binding(args):
    if args.pfmfile is not None:
        raise(NotImplementedError("pfmfile"))
    
    predict_peaks(
        args.outfile,
        atac_bams=args.atac_bams,
        histone_bams=args.histone_bams,
        regions=args.regions,
        factors=args.factors,
        genome=args.genome
    )
