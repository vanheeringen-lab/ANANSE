#!/usr/bin/env python
# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from loguru import logger

from ananse.peakpredictor import predict_peaks
from ananse.utils import check_path, check_input_factors


@logger.catch
def binding(args):
    predict_peaks(
        check_path(args.outdir, error_missing=False),
        atac_bams=check_path(args.atac_bams),
        histone_bams=check_path(args.histone_bams),
        regionfiles=check_path(args.regionfiles),
        reference=check_path(args.reference),
        factors=check_input_factors(args.factors),
        genome=args.genome,  # checked in CLI
        pfmfile=check_path(args.pfmfile),
        pfmscorefile=check_path(args.pfmscorefile),
        ncpus=args.ncpus,
    )
