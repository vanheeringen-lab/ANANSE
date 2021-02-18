import multiprocessing as mp
import os
import tempfile
import shutil

import dask.dataframe as dd
import dask.diagnostics
import genomepy
from gimmemotifs.scanner import scan_regionfile_to_table
from gimmemotifs.utils import pfmfile_location
from loguru import logger
import numpy as np
import pandas as pd
import pickle
import qnorm
from scipy import stats
from sklearn.preprocessing import minmax_scale

from ananse.utils import (
    bed_sort,
    bed_merge,
    bam_index,
    bam_sort,
    mosdepth,
)
from ananse.distributions import Distributions


class CombineBedFiles:
    def __init__(self, genome, peakfiles, verbose=True):
        self.genome = genome
        self.list_of_peakfiles = (
            peakfiles if isinstance(peakfiles, list) else [peakfiles]
        )
        self.verbose = verbose

    @staticmethod
    def is_narrowpeak(bed, check_values=True):
        """
        Check BED type by column count.
        Check if peak values are not all zeroes unless check_values is False.

        Accepts a BED file (including narrowPeak, broadPeak, etc.)
        Returns bool
        """
        with open(bed) as b:
            for line in b:
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

    @staticmethod
    def bed_resize(
        genome,
        bed_in,
        bed_out,
        width=200,
        narrowpeak=False,
        fix_outliers=False,
        output_bed3=True,
        verbose=True,
    ):
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

        if missing_chrm and verbose:
            logger.warning(
                "The following contigs were present in "
                + f"'{os.path.basename(bed_in)}', "
                + "but were missing in the genome file: "
                + f"{', '.join(list(set(missing_chrm)))}\n"
            )
        return bed_out

    def run(self, outfile, width=200, force=False):
        if force or not os.path.exists(outfile):
            if self.verbose:
                logger.info("Combining bed files")
            tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
            try:
                list_of_beds = []
                for peakfile in self.list_of_peakfiles:
                    # use narrowPeak Peak location for region centering if possible
                    is_np = self.is_narrowpeak(peakfile)
                    resized_peakfile = os.path.join(tmpdir, os.path.basename(peakfile))

                    # resize each BED region to 200 BP
                    self.bed_resize(
                        genome=self.genome,
                        bed_in=peakfile,
                        bed_out=resized_peakfile,
                        width=width,
                        narrowpeak=is_np,
                        verbose=self.verbose,
                    )
                    bed_sort(resized_peakfile)
                    list_of_beds.append(resized_peakfile)

                # merge resized beds into one
                merged_bed = os.path.join(tmpdir, "merged")
                bed_merge(list_of_beds=list_of_beds, merged_bed=merged_bed)

                shutil.copy2(merged_bed, outfile)
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)


class ScorePeaks:
    def __init__(self, bams, bed, ncore=1, verbose=True):
        self.list_of_bams = bams if isinstance(bams, list) else [bams]
        self.bed = bed  # one bed file with all putative enhancer binding regions
        self.verbose = verbose
        self.ncore = ncore

    def peaks_count(self, outdir):
        """
        count bam reads in the bed regions
        returns one bed file for each bam in outdir
        """
        # linear script:
        # for bam in self.list_of_bams:
        #     bed_output = os.path.join(outdir, os.path.basename(bam).replace(".bam", ".regions.bed"))
        #     mosdepth(self.bed, bam, bed_output, self.ncore)

        # parallel script:
        nbams = len(self.list_of_bams)
        npool = min(self.ncore, nbams)
        ncore = min(4, self.ncore // npool)  # 1-4 cores/bam

        # list with tuples. each tuple = one run
        mosdepth_params = []
        for bam in self.list_of_bams:
            bed_output = os.path.join(
                outdir, os.path.basename(bam).replace(".bam", ".regions.bed")
            )
            mosdepth_params.append((self.bed, bam, bed_output, ncore))

        pool = mp.Pool(npool)
        try:
            pool.starmap_async(mosdepth, mosdepth_params)
        finally:  # To make sure processes are closed in the end, even if errors happen
            pool.close()
            pool.join()

    @staticmethod
    def peaks_merge(indir, bed_output, ncore=1):
        """
        averages all peaks_count outputs
        uses quantile normalization to normalize for read depth
        returns one BED 3+1 file
        """
        ncore = min(4, ncore)
        files = [
            os.path.join(indir, f)
            for f in os.listdir(indir)
            if f.endswith(".regions.bed")
        ]
        bed = pd.read_csv(files[0], header=None, sep="\t")
        if len(files) > 1:
            for file in files[1:]:
                scores = pd.read_csv(file, header=None, sep="\t")[3]
                bed = pd.concat([bed, scores], axis=1)

            scores = bed.iloc[:, 3:]
            scores = qnorm.quantile_normalize(scores, axis=1, ncpus=ncore)
            scores = scores.mean(axis=1)

            bed = pd.concat([bed.iloc[:, :3], scores], axis=1)
        bed.to_csv(bed_output, sep="\t", header=False, index=False)

    @staticmethod
    def peaks_fit(bam_coverage, bed_output, dist_func="lognorm_dist", **kwargs):
        """
        fit the peak scores to a distribution
        """
        bed = pd.read_csv(bam_coverage, header=None, sep="\t")
        region = (
            bed[0].astype(str) + ":" + bed[1].astype(str) + "-" + bed[2].astype(str)
        )
        score = bed[3]

        # obtain a distribution
        dist_func = Distributions().set(dist_func)
        # with np.errstate(divide="ignore", invalid="ignore"):
        #     dist = dist_func(score, **kwargs)
        dist = dist_func(score, **kwargs)

        # replace scores with distribution values
        ascending_dist = np.sort(dist)
        ascending_scores_index = np.searchsorted(np.sort(score), score)
        norm_score = np.array([ascending_dist[i] for i in ascending_scores_index])

        logn_score = np.log(norm_score + 1)
        scaled_score = minmax_scale(logn_score)
        log10_score = np.log10(norm_score + 1)

        data = {
            "region": region,  # ex: "chr1:0-200"
            "score": score,
            "norm_score": norm_score,
            "logn_score": logn_score,
            "scaled_score": scaled_score,
            "log10_score": log10_score,  # used by the original function
        }
        bed = pd.DataFrame(data=data)
        bed.to_csv(bed_output, sep="\t", index=False)

    def run(self, outfile, dist_func="peak_rank_file_dist", force=False, **kwargs):
        # save the results as it takes ages to run
        raw_peak_scores = os.path.join(os.path.dirname(outfile), "raw_scoredpeaks.bed")
        if force or not os.path.exists(raw_peak_scores):
            tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
            try:
                if self.verbose:
                    logger.info("Scoring peaks (slow)")
                try:  # assumes sorted
                    for bam in self.list_of_bams:
                        bam_index(bam, force=False, ncore=self.ncore)
                    self.peaks_count(tmpdir)
                except Exception:  # sort, index & try again
                    for bam in self.list_of_bams:
                        bam_sort(bam, self.ncore)
                    self.peaks_count(tmpdir)

                tmp_peak_scores = os.path.join(tmpdir, "raw_scoredpeaks.bed")
                self.peaks_merge(tmpdir, tmp_peak_scores, self.ncore)

                shutil.copy2(tmp_peak_scores, raw_peak_scores)
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        # fit bam read counts to specified distribution
        if force or not os.path.exists(outfile):
            self.peaks_fit(raw_peak_scores, outfile, dist_func=dist_func, **kwargs)


class ScoreMotifs:
    def __init__(self, genome, bed, pfmfile=None, ncore=1, verbose=True):
        self.genome = genome
        self.bed = bed  # putative enhancer regions in format chr:start-end (in column 0 with header)
        self.pfm_file = pfmfile_location(pfmfile)
        self.ncore = ncore
        self.verbose = verbose

    def motifs_get_scores(self, pfmscorefile, debug=False):
        """
        Scan for TF binding motifs in potential enhancer regions.
        """
        if not debug:
            df = scan_regionfile_to_table(
                input_table=self.bed,
                genome=self.genome,
                scoring="score",
                pfmfile=self.pfm_file,
                ncpus=self.ncore,
                zscore=True,
                gc=True,
            )
        else:  # test output
            df = pd.DataFrame(
                {
                    "region": ["chr1:400-600", "chr1:2400-2600", "chr1:10003-10203"],
                    "GM.5.0.Sox.0001": [-0.544, -2.496, -0.544],
                    "GM.5.0.Homeodomain.0001": [-0.750, -0.377, -7.544],
                }
            ).set_index("region")

        df["motif"] = df.idxmax(axis=1)
        df["zscore"] = df.max(axis=1)
        df.reset_index(inplace=True)

        df.to_csv(
            pfmscorefile,
            sep="\t",
            header=True,
            index=False,
            columns=["motif", "region", "zscore"],  # filter + order columns
        )

    @staticmethod
    def motifs_normalize(bed_input, bed_output):
        """
        Add normalized scores to the scored motifs
        """
        bed = pd.read_csv(bed_input, sep="\t")
        bed["rank_zscore"] = minmax_scale(stats.rankdata(bed["zscore"]))
        bed.to_csv(bed_output, sep="\t", index=False)

    def run(self, outfile, force=False):
        raw_motif_scores = os.path.join(
            os.path.dirname(outfile), "raw_scoredmotifs.bed"
        )
        if force or not os.path.exists(raw_motif_scores):
            if self.verbose:
                logger.info("Scoring motifs (really slow)")
            self.motifs_get_scores(raw_motif_scores)
            # tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
            # try:
            #     if self.verbose:
            #         logger.info("Scoring motifs (really slow)")
            #     tmp_motif_scores = os.path.join(tmpdir, "motif_scores")
            #     self.motifs_get_scores(tmp_motif_scores)
            #
            #     shutil.copy2(tmp_motif_scores, raw_motif_scores)
            # finally:
            #     shutil.rmtree(tmpdir, ignore_errors=True)

        if force or not os.path.exists(outfile):
            self.motifs_normalize(raw_motif_scores, outfile)


class Binding:
    def __init__(
        self,
        peak_weights,
        motif_weights,
        pfmfile=None,
        model=None,
        curation_filter=None,
        tf_list=None,
        whitelist=True,
        ncore=1,
        verbose=True,
    ):
        self.peak_weights = peak_weights  # output from ScorePeaks
        self.motif_weights = motif_weights  # output from ScoreMotifs

        self.motifs2factors_file = pfmfile_location(pfmfile).replace(
            ".pfm", ".motif2factors.txt"
        )
        self.motifs2factors = self.filter_transcription_factors(
            curation_filter, tf_list, whitelist
        )

        self.model = model
        if self.model is None:
            # dream_model.txt is a 2D logistic regression model.
            package_dir = os.path.dirname(__file__)
            self.model = os.path.join(package_dir, "db", "dream_model_p300.pickle")

        self.ncore = ncore
        self.verbose = verbose

    def filter_transcription_factors(
        self, curation_filter=None, tf_list=None, whitelist=True
    ):
        """
        filter transcription factors from the motif database

        curation_filter: If None (default), keep all factors.
        If True, keep only curated factors. If False, keep only non-curated factors.
        Note: "Curated" TFs have direct evidence for binding or are manually selected for likely binding.

        tf_list: an optional, single-column file with (case-insensitive) transcription factor names.
        whitelist: if True (default), tf_list is used as a whitelist. If False, as a blacklist.
        """
        m2f = pd.read_csv(self.motifs2factors_file, sep="\t")

        # rename stuff
        m2f.rename(
            columns={"Motif": "motif", "Factor": "factor", "Curated": "curated"},
            inplace=True,
        )
        m2f["factor"] = m2f.factor.str.upper()  # make case-insensitive
        m2f.replace("T", "TBXT", inplace=True)  # rename T to TBXT

        # filter by curation
        if curation_filter is True:
            m2f = m2f.loc[m2f.curated == "Y"]
        elif curation_filter is False:
            m2f = m2f.loc[m2f.curated == "N"]

        # shrink table
        m2f = m2f[["motif", "factor"]]  # subset
        m2f.drop_duplicates(
            inplace=True
        )  # case-insensitivity adds loads of duplicates (ex: Sox9 and SOX9)

        # filter by white/blacklist
        if tf_list:
            tfs = (
                pd.read_csv(tf_list, header=None)[0].str.upper().tolist()
            )  # make case-insensitive
            m2f = (
                m2f.loc[m2f.factor.isin(tfs)]
                if whitelist
                else m2f.loc[~m2f.factor.isin(tfs)]
            )

        return m2f

    def get_binding_score(self, motif_weights, peak_weights, outfile):
        """
        Infer TF binding score from motif z-score and peak intensity.
        """
        # merge the scoring tables
        m = dd.read_csv(motif_weights, sep="\t")
        m = m.merge(dd.read_csv(peak_weights, sep="\t", blocksize=200e6), on="region")[
            ["motif", "region", "zscore", "log10_score"]
        ]

        # filter scoring tables for motifs found in motifs2factors
        m = m.merge(self.motifs2factors, on="motif")  # also adds "factor" column

        # combine scores
        m = m.groupby(["factor", "region"])[["zscore", "log10_score"]].mean()
        m = m.dropna().reset_index()

        with dask.diagnostics.ProgressBar():
            m = m.compute(num_workers=self.ncore)

        # Load model
        with open(self.model, "rb") as f:
            clf = pickle.load(f)
        m["binding"] = clf.predict_proba(m[["zscore", "log10_score"]])[:, 1]

        # "region" renames to "enhancer" for consistency with ANANSE network
        m.rename(columns={"region": "enhancer"}, inplace=True)
        m.to_csv(
            outfile, sep="\t", index=False, columns=["factor", "enhancer", "binding"]
        )

    def run(self, outfile, force=False):
        if force or not os.path.exists(outfile):
            if self.verbose:
                logger.info("Predict TF binding")
            self.get_binding_score(self.peak_weights, self.motif_weights, outfile)
