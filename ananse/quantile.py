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
from tempfile import NamedTemporaryFile
import subprocess

from ananse import mytmpdir
import ananse

class Quantile(object):
    def __init__(self, bam_input, broadPeak, bed_output):
        
        package_dir = os.path.dirname(ananse.__file__)

        self.bam_input = bam_input
        self.broadPeak = broadPeak
        self.bed_output = bed_output
        self.peak_2k =  os.path.join(package_dir, "db", "enhancer_peak_2000.bed")
        self.peak_rank =  os.path.join(package_dir, "db", "peak_rank.txt")

    def runCov(self, bam_input):
        covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        covcmd = f"multiBamCov -bams {bam_input} -bed {self.peak_2k} > {covfile.name}"
        process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        return covfile.name

    def quantileNormalize(self, bed_input):
        rank=[]
        with open(self.peak_rank) as p:
            for i in p:
                rank.append(float(i[:-1]))

        bed = pd.read_csv(bed_input, header=None, sep="\t")
        t = np.searchsorted(np.sort(bed[3]), bed[3])
        bed[3] = [rank[i] for i in t]
        bed[1] = [int(i)+900 for i in bed[1].tolist()]
        bed[2] = [int(i)-900 for i in bed[2].tolist()]
        
        quantile_bed = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        bed.to_csv(quantile_bed,sep="\t",header=False,index=False)
        return quantile_bed.name

    def runIntersect(self, broadPeak, bed_input, bed_output):
        intercmd = f"bedtools intersect -a {bed_input} -b {broadPeak} -wa > {bed_output}"
        process = subprocess.Popen(intercmd, shell=True, stdout=subprocess.PIPE)
        process.wait()

    def run_quantile(self, bed_input, broadPeak, bed_output):
        bed_input = self.runCov(self.bam_input)
        quantile_bed = self.quantileNormalize(bed_input)
        self.runIntersect(self.broadPeak, quantile_bed, bed_output)

