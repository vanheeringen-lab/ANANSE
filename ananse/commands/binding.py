#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from ananse.peakpredictor import predict_peaks


def binding(args):
    predict_peaks(
        args.outdir,
        atac_bams=args.atac_bams,
        histone_bams=args.histone_bams,
        regionfiles=args.regionfiles,
        reference=args.reference,
        factors=args.factors,
        genome=args.genome,
        pfmfile=args.pfmfile,
        pfmscorefile=args.pfmscorefile,
        ncpus=args.ncpus,
    )
