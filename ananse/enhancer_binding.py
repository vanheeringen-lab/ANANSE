import os
import tempfile
import shutil
import warnings

import genomepy
import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy import stats
from sklearn.preprocessing import minmax_scale


def shhh_bedtool(func):
    """
    Decorator that silences this pybedtools warning:
    RuntimeWarning: line buffering (buffering=1) isn't supported in binary mode
    """
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            func(*args, **kwargs)
    return wrapper


def is_narrowpeak(bed, check_values=True):
    """
    Check BED type by column count.
    Check if peak values are not all zeroes unless check_values is False.

    Accepts a BED file (including narrowPeak, broadPeak, etc.)
    Returns bool
    """
    with open(bed) as b:
        for n, line in enumerate(b):
            if line.startswith("#"):
                continue
            line = line.split("\t")
            cols = len(line)
            break

    # narrowPeak has 10 columns
    # and the peak column is >= 0
    if cols != 10 or int(line[9]) < 0:
        return False

    if not check_values:
        return True

    # check if the peak values aren't all zeroes
    summit_values = 0
    sample_size = 20  # check an arbitrary number of lines
    with open(bed) as b:
        for n, line in enumerate(b):
            if line.startswith("#"):
                continue
            line = line.split("\t")
            peak_val = int(line[9])

            # value must be >=0
            if peak_val < 0:
                return False
            summit_values += peak_val
            if n >= sample_size:
                break

    if summit_values > 0:
        return True
    return False


def set_bed_width(genome, bed_in, bed_out, width=200, narrowpeak=False, fix_outliers=False, output_bed3=True):
    """
    Set bed region width.

    If the input bed is a narrowPeak file (narrowpeak=True),
    center region on the summit (start+peak).
    Otherwise center on the middle of the region.

    If fix_outliers is set to True, shift regions to fit their chromosomes.
    Otherwise drop these regions.

    If output_bed3 is set to False, output the whole bed file.
    """
    half_seqlen = width // 2
    chrom_sizes = genomepy.Genome(genome).sizes
    missing_chrm = []

    if narrowpeak:
        def get_summit(_start, _, summit_offset):
            return _start + int(summit_offset)
        summit_col = 9
    else:
        def get_summit(_start, _end, _):
            return (_start + _end) // 2
        summit_col = 0  # unused

    with open(bed_in) as old, open(bed_out, "w") as new:
        for line in old:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            chrm = str(line[0])
            if chrm not in chrom_sizes.keys():
                missing_chrm.append(chrm)
                continue
            start = int(line[1])
            end = int(line[2])
            rest = line[3:] if not output_bed3 else []

            chrm_len = chrom_sizes[chrm]
            if width == end - start:
                nstart = str(start)
                nend = str(end)
            elif chrm_len <= width:
                if not fix_outliers:
                    continue
                nstart = str(0)
                nend = str(chrm_len)
            else:
                summit = get_summit(start, end, line[summit_col])
                if not fix_outliers:
                    nstart = str(summit - half_seqlen)
                    nend = str(summit + half_seqlen)
                    if int(nstart) < 0 or int(nend) > chrm_len:
                        continue
                else:
                    # adjust the summit for the chromosome boundaries
                    summit = max(summit, 0 + half_seqlen)
                    summit = min(summit, chrm_len - half_seqlen)

                    nstart = str(summit - half_seqlen)
                    nend = str(summit + half_seqlen)

            new.write("\t".join([chrm, nstart, nend] + rest) + "\n")

    if missing_chrm:
        print(
            f"The following contigs were present in "
            f"'{os.path.basename(bed_in)}',\n"
            f"but were missing in the genome file:",
            ", ".join(list(set(missing_chrm))), "\n"
        )
    return bed_out


@shhh_bedtool
def sort_bed(bed):
    """
    Sort a bed file.
    """
    tmpdir = tempfile.mkdtemp()
    tmpfile = os.path.join(tmpdir, "tmpfile")
    try:
        BedTool(bed).sort(output=tmpfile)
        shutil.copy2(tmpfile, bed)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@shhh_bedtool
def merge_bed_files(list_of_beds, merged_bed):
    """
    merge any number of bed files
    """
    bed = BedTool(list_of_beds[0])
    if list_of_beds[1:]:
        bed = bed.cat(*list_of_beds[1:])
    # bed.sort(output=merged_bed)
    bed.moveto(merged_bed)


class CombinePeakFiles:
    def __init__(self, genome, peakfiles, outfile):
        self.genome = genome
        self.list_of_peakfiles = peakfiles if isinstance(peakfiles, list) else [peakfiles]
        self.outfile = outfile

    def resize_n_combine_peakfiles(self, width=200):
        tmpdir = tempfile.mkdtemp()
        try:
            list_of_beds = []
            for peakfile in self.list_of_peakfiles:
                # use narrowPeak Peak location for region centering if possible
                narrowpeak = is_narrowpeak(peakfile)
                resized_peakfile = os.path.join(tmpdir, os.path.basename(peakfile))

                # resize each BED region to 200 BP
                set_bed_width(genome=self.genome, bed_in=peakfile, bed_out=resized_peakfile, width=width, narrowpeak=narrowpeak)
                sort_bed(resized_peakfile)
                list_of_beds.append(resized_peakfile)

            # merge resized beds into one
            merged_bed = os.path.join(tmpdir, "merged")
            merge_bed_files(list_of_beds=list_of_beds, merged_bed=merged_bed)

            sort_bed(merged_bed)
            shutil.copy2(merged_bed, self.outfile)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def run(self, width=200):
        self.resize_n_combine_peakfiles(width=width)


@shhh_bedtool
def count_reads(bam_list, peakfile, bed_output):
    """
    Count bam reads in putative enhancer regions
    """
    bed = BedTool(peakfile)
    bed.multi_bam_coverage(bams=bam_list, output=bed_output)


def sum_coverages(multi_bam_coverage, sum_bam_coverage):
    """
    MultiBamCov returns a BED3+n with one column per bam.
    This function sums up all bam counts and returns a BED3+1.
    """
    bed = pd.read_csv(multi_bam_coverage, header=None, sep="\t")
    columns = bed.shape[1]
    if columns != 4:
        bed3 = bed.iloc[:, :3]
        scores = bed.iloc[:, 3:].sum(axis=1)
        bed = pd.concat([bed3, scores], axis=1)
    bed.to_csv(sum_bam_coverage, sep="\t", header=False, index=False)


def lognorm_dist(scores, **kwargs):
    """
    fit scores to a log normal distribution
    """
    scores = scores + 1  # add pseudocount
    s, loc, scale = stats.lognorm.fit(scores, floc=0)
    x = range(len(scores))
    dist = stats.lognorm.pdf(x=x, s=s, loc=loc, scale=scale)
    return dist


def normalize(bam_coverage, bed_output, dist_func=lognorm_dist, **kwargs):
    """
    Fit the bam coverage scores the a distribution
    """
    bed = pd.read_csv(bam_coverage, header=None, sep="\t")
    bed3 = bed.iloc[:, :3]
    scores = bed[3]

    # obtain a distribution
    dist = dist_func(scores, **kwargs)

    # replace scores with distribution values
    ascending_dist = np.sort(dist)
    ascending_scores_index = np.searchsorted(np.sort(scores), scores)
    norm_scores = [ascending_dist[i] for i in ascending_scores_index]

    scaled = minmax_scale(norm_scores)
    ranked = minmax_scale(stats.rankdata(norm_scores))

    bed = bed3
    for col in [norm_scores, scaled, ranked]:
        bed = pd.concat([bed, pd.Series(col)], axis=1)
    bed.columns = ["chrom", "start", "end", "logscore", "scaledscore", "rankedscore"]
    bed.to_csv(bed_output, sep="\t", index=False)


class ScorePeaks:
    def __init__(self, bams, peaks, outfile):
        self.list_of_bams = bams if isinstance(bams, list) else [bams]
        self.peaks = peaks
        self.outfile = outfile

    def count_reads(self, coverage_file):
        tmpdir = tempfile.mkdtemp()
        try:
            # bam read counts per bam file
            multi_bam_coverage = os.path.join(tmpdir, "multi_bam_coverage")
            count_reads(self.list_of_bams, self.peaks, multi_bam_coverage)

            # sum bam read counts
            sum_bam_coverage = os.path.join(tmpdir, "sum_bam_coverage")
            sum_coverages(multi_bam_coverage, sum_bam_coverage)

            shutil.copy2(sum_bam_coverage, coverage_file)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def run(self, force=False, dist_func=lognorm_dist):
        # save the results as it takes ages to run
        coverage_file = os.path.join(os.path.dirname(self.outfile), "raw_enhancer_coverage.bed")
        if force or not os.path.exists(coverage_file):
            print("running multi bam cov (slow)")
            self.count_reads(coverage_file)

        # fit bam read counts to specified distribution
        normalize(coverage_file, self.outfile, dist_func=dist_func)
