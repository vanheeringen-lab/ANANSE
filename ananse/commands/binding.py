#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import sys
import os

from ananse.peakpredictor import predict_peaks


def binding(args):
    if args.pfmfile is not None:
        raise(NotImplementedError("pfmfile"))
    if args.regions is not None:
        raise(NotImplementedError("regions"))
    if args.genome is not "hg38":
        raise(NotImplementedError("genome"))
    if args.factors is not None:
        raise(NotImplementedError("factors"))
    
    predict_peaks(
        args.outfile,
        atac_bams=args.atac_bams,
        histone_bams=args.histone_bams,
        regions=args.regions,
        factors=args.factors,
    )
