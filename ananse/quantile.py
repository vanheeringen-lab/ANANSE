#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

"""Quantile normalize enhancer reads"""

import os

import numpy as np
import pandas as pd

from ananse import mytmpdir
import ananse

class Quantile(object):
    def __init__(self, bed_input, bed_output):
        
        package_dir = os.path.dirname(ananse.__file__)

        self.bed_input = bed_input
        self.bed_output = bed_output
        self.peak_rank =  os.path.join(package_dir, "db", "peak_rank.txt")

    def quantileNormalize(self, bed_input, bed_output):
        rank=[]
        with open(self.peak_rank) as p:
            for i in p:
                rank.append(float(i[:-1]))

        bed = pd.read_csv(bed_input, header=None, sep="\t")
        t = np.searchsorted(np.sort(bed[3]), bed[3])
        bed[3] = [rank[i] for i in t]
        bed[1] = [int(i)+900 for i in bed[1].tolist()]
        bed[2] = [int(i)-900 for i in bed[2].tolist()]
        bed.to_csv(bed_output,sep="\t",header=False,index=False)

    def run_quantile(self, bed_input, bed_output):
        self.quantileNormalize(bed_input, bed_output)

