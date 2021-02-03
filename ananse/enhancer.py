#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import os
import logging

from pybedtools import BedTool
from genomepy import Genome, utils
from ananse.utils import set_width

# use logger settings defined in __init__
logger = logging.getLogger(__name__)


class Enhancer:
    """
     Create a BED3+1 file, the 4th column containing a scoring method.

     The regions in the peak file(s) are normalized to a specified width.
     The output width of the regions is set to a specified width, centered around a summit if possible.
     Regions that don't fit after normalization are shifted/shrunk.

     Reads in the bam files are summed for the peak regions, then normalized to an experimentally determined distribution.
     """

    def __init__(self, genome, bams, peaks, output, width=200, summit_column=-1):
        self.genome = genome
        self.bams = bams
        self.peaks = peaks
        self.output = output
        self.outdir = os.path.dirname(self.output)
        self.width = width
        self.summit_column = summit_column

    @staticmethod
    def combine_peaks(peaks_list, peaks_file):
        """
        sum bedfiles
        """
        raise NotImplementedError()

    @staticmethod
    def count_reads(bam_list, peaks_file, bed_output):
        """
        Count bam reads in putative enhancer regions
        """
        bed = BedTool(peaks_file)
        bed.multi_bam_coverage(bams=bam_list, output=bed_output)

    @staticmethod
    def normalize(self, bed_input, bed_output):
        """
        Fit the scores to the model-compatible distribution
        """
        raise NotImplementedError()

    def run_enhancer(self):
        utils.mkdir_p(os.path.join(self.outdir, "intermediate"))

        combined_peaks = os.path.join(self.outdir, "intermediate", "combined_peaks.bed")
        self.combine_peaks(self.peaks, combined_peaks)

        normalized_peaks = os.path.join(self.outdir, "intermediate", "normalized_peaks.bed")
        set_width(genome=self.genome, bed_in=combined_peaks, bed_out=normalized_peaks, witdh=self.width, summit_col=self.summit_column)

        bam_reads = os.path.join(self.outdir, "intermediate", "bam_reads")
        self.count_reads(self.bams, normalized_peaks, bam_reads)

        self.normalize(bam_reads, self.output)

        if os.path.exists(self.output):
            utils.rm_rf(os.path.join(self.outdir, "intermediate"))

# class Enhancer(object):
#     def __init__(self, bam_input, epeak, bed_output, genome="hg38"):
#
#         package_dir = os.path.dirname(ananse.__file__)
#
#         self.genome = genome
#         self.bam_input = bam_input
#         self.epeak = epeak
#         self.bed_output = bed_output
#         self.peak_2k =  os.path.join(package_dir, "db", "enhancer_peak_2000.bed")
#         self.peak_rank =  os.path.join(package_dir, "db", "peak_rank_hg38_h3k27ac.txt")
#
#     def runCov(self, bam_input):
#         covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         covcmd = f"multiBamCov -bams {bam_input} -bed {self.peak_2k} > {covfile.name}"
#         process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#         return covfile.name
#
#     def quantileNormalize(self, bed_input):
#         rank=[]
#         with open(self.peak_rank) as p:
#             for i in p:
#                 rank.append(float(i[:-1]))
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         bed[1] = [int(i)+900 for i in bed[1].tolist()]
#         bed[2] = [int(i)-900 for i in bed[2].tolist()]
#
#         quantile_bed = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         bed.to_csv(quantile_bed, sep="\t", header=False, index=False)
#         return quantile_bed.name
#
#     def runIntersect(self, epeak, bed_input, bed_output):
#         intercmd = f"bedtools intersect -a {bed_input} -b {epeak} -wa > {bed_output}"
#         process = subprocess.Popen(intercmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#
#     def run_enhancer(self, bed_input, epeak, bed_output):
#         bed_input = self.runCov(self.bam_input)
#         quantile_bed = self.quantileNormalize(bed_input)
#         self.runIntersect(self.epeak, quantile_bed, bed_output)
#
# class P300Enhancer(object):
#     def __init__(self, bam_input, epeak, bed_output, genome="hg38"):
#
#         package_dir = os.path.dirname(ananse.__file__)
#
#         self.genome = genome
#         self.bam_input = bam_input
#         self.epeak = epeak
#         self.bed_output = bed_output
#         self.peak_rank =  os.path.join(package_dir, "db", "peak_rank.txt")
#
#     def set_peak_size(self, peak_bed, seqlen=200):
#
#         """set all input peaks to 200bp
#         Arguments:
#             peak_bed {[bed]} -- [input peak bed file]
#
#         Keyword Arguments:
#             seqlen {int} -- [peak length] (default: {200})
#
#         Returns:
#             [type] -- [200bp peak file]
#         """
#         gsizedic = Genome(self.genome).sizes
#
#         peaks = BedTool(peak_bed)
#         fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#
#         for peak in peaks:
#             if peak.length < seqlen or peak.length > seqlen:
#                 # get the summit and the flanking low and high sequences
#                 summit = (peak.start + peak.end) // 2
#                 start, end = summit - seqlen // 2, summit + seqlen // 2
#             else:
#                 start, end = peak.start, peak.end
#
#             # remove seq which langer than chromosome length or smaller than 0
#             if start > 0 and end < int(gsizedic[peak.chrom]):
#                 fl2.write(f"{peak.chrom}\t{start}\t{end}\n")
#
#         return fl2.name
#
#     def mk_peak(self, epeak):
#         epeak200 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         with open(epeak) as peakfile, open(epeak200.name,"w") as npeakfile:
#             for line in peakfile:
#                 a=line.split()
#                 chrm=a[0]
#                 start=int(a[1])
#                 summit=int(a[9])
#                 nsummit=start+summit
#                 if nsummit<100:
#                     nsummit=100
#                 npeakfile.write(f"{chrm}\t{nsummit-100}\t{nsummit+100}\n")
#         return epeak200.name
#
#     def quantileNormalize(self, epeak, bed_input, bed_output):
#         enahcer_number = pd.read_csv(epeak, sep="\t", header=None).shape[0]
#         rank = pd.read_csv(self.peak_rank, header=None).sample(n = enahcer_number, random_state = 1).sort_values(0)[0].tolist()
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         bed[1] = [int(i) for i in bed[1].tolist()]
#         bed[2] = [int(i) for i in bed[2].tolist()]
#
#         bed.to_csv(bed_output, sep="\t", header=False, index=False)
#
#     def runCov(self, bam_input, clear_epeak200):
#         covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         covcmd = f"multiBamCov -bams {bam_input} -bed {clear_epeak200} > {covfile.name}"
#         process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#         return covfile.name
#
#     def run_enhancer(self, bed_input, epeak, bed_output):
#         epeak200 = self.mk_peak(epeak)
#         clear_epeak200 = self.set_peak_size(epeak200)
#         # os.system(f"cp {clear_epeak200} ./")
#         bed_cov = self.runCov(self.bam_input, clear_epeak200)
#         quantile_bed = self.quantileNormalize(epeak, bed_cov, bed_output)
#
#
# class AtacEnhancer(object):
#     def __init__(self, bam_input, epeak, bed_output, genome="hg38"):
#
#         package_dir = os.path.dirname(ananse.__file__)
#
#         self.genome = genome
#         self.bam_input = bam_input
#         self.epeak = epeak
#         self.bed_output = bed_output
#         self.peak_rank =  os.path.join(package_dir, "db", "peak_rank.txt")
#
#     def set_peak_size(self, peak_bed, seqlen=200):
#         """set all input peaks to 200bp
#
#         Arguments:
#             peak_bed {[bed]} -- [input peak bed file]
#
#         Keyword Arguments:
#             seqlen {int} -- [peak length] (default: {200})
#
#         Returns:
#             [type] -- [200bp peak file]
#         """
#         gsizedic = Genome(self.genome).sizes
#
#         peaks = BedTool(peak_bed)
#         fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#
#         for peak in peaks:
#
#             if peak.length < seqlen or peak.length > seqlen:
#                 # get the summit and the flanking low and high sequences
#                 summit = (peak.start + peak.end) // 2
#                 start, end = summit - seqlen // 2, summit + seqlen // 2
#             else:
#                 start, end = peak.start, peak.end
#             # remove seq which langer than chromosome length or smaller than 0
#             if start > 0 and end < int(gsizedic[peak.chrom]):
#                 fl2.write(f"{peak.chrom}\t{start}\t{end}\n")
#
#         # return npeaks
#         return fl2.name
#
#     def mk_peak(self, epeak):
#         epeak200 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         with open(epeak) as peakfile, open(epeak200.name,"w") as npeakfile:
#             for line in peakfile:
#                 a=line.split()
#                 chrm=a[0]
#                 start=int(a[1])
#                 summit=int(a[9])
#                 nsummit=start+summit
#                 if nsummit<100:
#                     nsummit=100
#                 npeakfile.write(f"{chrm}\t{nsummit-100}\t{nsummit+100}\n")
#         return epeak200.name
#
#     def quantileNormalize(self, epeak, bed_input, bed_output):
#         enahcer_number = pd.read_csv(epeak, sep="\t", header=None).shape[0]
#         rank = pd.read_csv(self.peak_rank, header=None).sample(n = enahcer_number, random_state = 1).sort_values(0)[0].tolist()
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         bed[1] = [int(i)+900 for i in bed[1].tolist()]
#         bed[2] = [int(i)-900 for i in bed[2].tolist()]
#
#         bed.to_csv(bed_output, sep="\t", header=False, index=False)
#
#     def runCov(self, bam_input, clear_epeak2k):
#         covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         covcmd = f"multiBamCov -bams {bam_input} -bed {clear_epeak2k} > {covfile.name}"
#         process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#         return covfile.name
#
#     def run_enhancer(self, bed_input, epeak, bed_output):
#         epeak200 = self.mk_peak(epeak)
#         clear_epeak2k = self.set_peak_size(epeak200, 2000)
#         bed_cov = self.runCov(self.bam_input, clear_epeak2k)
#         quantile_bed = self.quantileNormalize(epeak, bed_cov, bed_output)
#
#
