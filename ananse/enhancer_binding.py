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

    peak_score = np.expm1(norm_scores)
    log10_peak_score = norm_scores
    Scaled_peak_score = minmax_scale(norm_scores)
    Ranked_peak_score = minmax_scale(stats.rankdata(norm_scores))

    # scaled = minmax_scale(norm_scores)
    # ranked = minmax_scale(stats.rankdata(norm_scores))

    bed = bed3
    # for col in [norm_scores, scaled, ranked]:
    for col in [peak_score, log10_peak_score, Scaled_peak_score, Ranked_peak_score]:
        bed = pd.concat([bed, pd.Series(col)], axis=1)
    bed.columns = ["chrom", "start", "end", "peak_score", "log10_peak_score", "Scaled_peak_score", "Ranked_peak_score"]
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

# TODO: still messy from here

from tqdm import tqdm
from gimmemotifs.scanner import Scanner
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import as_fasta, pfmfile_location
import pickle
import dask.dataframe as dd


class Binding(object):
    def __init__(self, genome, scored_peaks, outfile, model=None, pfmfile=None, ncore=1):
        self.genome = genome
        self.scored_peaks = scored_peaks
        self.outfile = outfile
        self.model = model
        self.ncore = ncore

        if self.model is None:
            # dream_model.txt is a 2D logistic regression model.
            package_dir = os.path.dirname(__file__)
            self.model = os.path.join(package_dir, "db", "dream_model_p300.pickle")

        # Motif information file
        self.pfmfile = pfmfile_location(pfmfile)
        self.motifs2factors = self.pfmfile.replace(".pfm", ".motif2factors.txt")

    def setup_gimme_scanner(self):
        s = Scanner(ncpus=self.ncore)
        s.set_motifs(self.pfmfile)
        s.set_genome(self.genome)
        s.set_threshold(threshold=0.0)

        # generate GC background index
        _ = s.best_score([], zscore=True, gc=True)
        return s

    def get_motif_scores(self, enhancer_regions_bed, pfmscorefile):
        """
        Scan for TF binding motifs in potential enhancer regions.
        """
        s = self.setup_gimme_scanner()
        with open(self.pfmfile) as f:
            motifs = read_motifs(f)

        # open an empty file (append results in chunks)
        with open(pfmscorefile, 'w') as f:
            # Quan: When we built model, rank and minmax normalization was used.
            cols = ["enhancer", "zscore", "zscoreRank"]
            f.write("\t".join(cols)+"\n")

        seqs = [s.split(" ")[0] for s in as_fasta(enhancer_regions_bed, genome=self.genome).ids]
        with tqdm(total=len(seqs)) as pbar:
            # Run 10k regions per scan.
            chunksize = 10000
            for chunk in range(0, len(seqs), chunksize):
                pfm_score = []
                chunk_seqs = seqs[chunk:chunk+chunksize]
                # We are using GC-normalization for motif scanning as many enhancer binding regions are GC-enriched.
                chunk_scores = s.best_score(chunk_seqs, zscore=True, gc=True)
                for seq, scores in zip(chunk_seqs, chunk_scores):
                    for motif, score in zip(motifs, scores):
                        pfm_score.append([motif.id, seq, score])
                    pbar.update(1)
                pfm_score = pd.DataFrame(pfm_score, columns=["motif", "enhancer", "zscore"])
                pfm_score = pfm_score.set_index("motif")
                pfm_score["zscoreRank"] = minmax_scale(stats.rankdata(pfm_score["zscore"]))

                pfm_score[cols].to_csv(pfmscorefile, sep="\t", header=False, mode='a')

    def get_binding_score(self, pfm, peak, outfile):
        """Infer TF binding score from motif z-score and peak intensity.
        Arguments:
            pfm {[type]} -- [motif scan result]
            peak {[type]} -- [peak intensity]
        Returns:
            [type] -- [the predicted tf binding table]
        """
        peak = dd.read_csv(peak, sep="\t", blocksize=200e6)
        pfm = dd.read_csv(pfm, sep="\t")

        # Load model
        with open(self.model, "rb") as f:
            clf = pickle.load(f)

        # ft = self.filtermotifs2factors

        r = pfm.merge(peak, left_on="enhancer", right_on="peak")[
            ["motif", "enhancer", "zscore", "log10_peakRPKM"]
        ]
        # r = r.merge(ft, left_on="motif", right_on="Motif")
        r = r.groupby(["factor", "enhancer"])[["zscore", "log10_peakRPKM"]].mean()
        r = r.dropna().reset_index()

        table = r.compute(num_workers=self.ncore)
        table["binding"] = clf.predict_proba(table[["zscore", "log10_peakRPKM"]])[:, 1]

        table.to_csv(outfile, sep="\t", index=False)

    def run(self):
        tmpdir = tempfile.mkdtemp()
        try:
            pfm_weight = os.path.join(tmpdir, "peak_weight")
            self.get_motif_scores(self.scored_peaks, pfm_weight)

            peak_weight = self.scored_peaks
            table = os.path.join(tmpdir, "table")
            self.get_binding_score(pfm_weight, peak_weight, table)

            shutil.copy2(table, self.outfile)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
