#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

"""Predict TF binding network"""

# Python imports
import os
import pickle
from tqdm import tqdm
import warnings
from tempfile import NamedTemporaryFile
from loguru import logger

import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
import dask.dataframe as dd
from scipy.stats import rankdata

from pybedtools import BedTool
from genomepy import Genome
from gimmemotifs.scanner import Scanner
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import as_fasta, pfmfile_location

from ananse import mytmpdir
import ananse

warnings.filterwarnings("ignore")


def clear_tfs(motifs2factors, tffile, include_notfs=False, rm_curated=True):
    """filter unreal TFs from motif database

    Arguments:
        motifs2factors {[type]} -- [motifs2factors]
        tffile {[type]} -- [real tfs]

    Returns:
        [type] -- [motifs2factors]
    """
    ft = pd.read_csv(motifs2factors, sep="\t")

    ft['Factor'] = ft.Factor.str.upper()
    if rm_curated:
        # "Curated" is manually curated or direct evidence for binding. For instance a ChIP-seq predicted motif is an N in this column
        ft = ft.loc[ft.Curated == "Y"]
    if not include_notfs:
        tfs = pd.read_csv(tffile, header = None)[0].tolist()
        ft = ft.loc[ft.Factor.isin (tfs)]
    # replace T to TBXT
    ft = ft.replace("T" , "TBXT")
    ft = ft.replace("t" , "tbxt")

    ft.rename(columns = {"Factor":"factor"}, inplace = True)
    return ft

class Binding(object):
    def __init__(self, ncore=1, genome="hg38", gene_bed=None, pfmfile=None, include_notfs=False, rm_curated=True, etype="hg38H3K27ac", tffile=None):

        self.ncore = ncore
        self.genome = genome

        # dream_model.txt is the logistic regression model.
        package_dir = os.path.dirname(ananse.__file__)
        self.etype = etype

        if self.genome == "hg38" and self.etype == "hg38H3K27ac":
            self.model = os.path.join(package_dir, "db", "dream_model_h3k27ac.txt")
        elif self.etype == "p300" or self.etype == "ATAC":
            self.model = os.path.join(package_dir, "db", "dream_model_p300.txt")
        else:
            raise TypeError("""The input enhancer data type should hg38H3K27ac, p300 or ATAC. 
            It is not possible set -e to hg38H3K27ac if the genome is not hg38. 
            Please provide a enhancer type with -e argument. By default is hg38H3K27ac.""")

        # filter tfs?
        self.include_notfs = include_notfs
        # remove curated?
        self.rm_curated = rm_curated

        # load real tfs
        self.tffile = tffile
        if self.tffile is None:
            self.tffile = os.path.join(package_dir, "db", "tfs.txt")
        # self.tffile = "db/tfs.txt"

        # Motif information file
        self.pfmfile = pfmfile_location(pfmfile)
        self.motifs2factors = self.pfmfile.replace(".pfm", ".motif2factors.txt")
        self.filtermotifs2factors = clear_tfs(self.motifs2factors, self.tffile, self.include_notfs, self.rm_curated)

    def set_peak_size(self, peak_bed, seqlen=200):
        """set all input peaks to 200bp

        Arguments:
            peak_bed {[bed]} -- [input peak bed file]

        Keyword Arguments:
            seqlen {int} -- [peak length] (default: {200})

        Returns:
            [type] -- [200bp peak file]
        """
        gsizedic = Genome(self.genome).sizes

        peaks = BedTool(peak_bed)
        fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)

        for peak in peaks:

            if peak.length < seqlen or peak.length > seqlen:
                # get the summit and the flanking low and high sequences
                summit = (peak.start + peak.end) // 2
                start, end = summit - seqlen // 2, summit + seqlen // 2
            else:
                start, end = peak.start, peak.end
            # remove seq which langer than chromosome length or smaller than 0
            if start > 0 and end < int(gsizedic[peak.chrom]):
                fl2.write(
                    str(peak.chrom)
                    + "\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t"
                    + str(peak.fields[-1])
                    + "\n"
                )
        # return npeaks
        return fl2.name

    def get_peakRPKM(self, fin_rpkm):
        # When we built model, the peak intensity was ranked and scaled.
        peaks = pd.read_table(fin_rpkm, names=["chrom", "start", "end", "peakRPKM"])
        peaks["peak"] = (
            peaks["chrom"]
            + ":"
            + peaks["start"].astype(str)
            + "-"
            + peaks["end"].astype(str)
        )
        add = peaks["peakRPKM"][peaks["peakRPKM"] > 0].min()
        peaks["log10_peakRPKM"] = np.log10(peaks["peakRPKM"] + add)
        peaks["peakRPKMScale"] = minmax_scale(peaks["log10_peakRPKM"])
        peaks["peakRPKMRank"] = minmax_scale(rankdata(peaks["log10_peakRPKM"]))

        peakrpkmfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        cols = ["peak", "peakRPKM", "log10_peakRPKM", "peakRPKMScale", "peakRPKMRank"]
        peaks[cols].to_csv(peakrpkmfile, sep="\t", index=False)

        return peakrpkmfile.name

    def get_PWMScore(self, fin_regions_fa):
        """ Scan motif in every peak.

        Arguments:
            fin_regions_fa {[type]} -- [input fasta file]

        Returns:
            [type] -- [pfmscorefile]
        """
        pfmscorefile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=self.genome).ids]

        s = Scanner(ncpus=self.ncore)
        s.set_motifs(self.pfmfile)
        s.set_threshold(threshold=0.0)
        s.set_genome(self.genome)

        with open(self.pfmfile) as f:
            motifs = read_motifs(f)

        chunksize = 10000
        # Run 10k peaks one time.

        with tqdm(total=len(seqs)) as pbar:
            for chunk in range(0, len(seqs), chunksize):
                chunk_seqs = seqs[chunk : chunk + chunksize]
                # print(chunk, "-", chunk + chunksize, "enhancers")
                pfm_score = []
                it = s.best_score(chunk_seqs, zscore=True, gc=True)
                # We are using GC-normalization for motif scan because many sequence is GC-enriched.
                # GimmeMotif develop branch already include GC-normalization option now.
                for seq, scores in zip(chunk_seqs, it):
                    for motif, score in zip(motifs, scores):
                        pfm_score.append([motif.id, seq, score])
                    pbar.update(1)
                pfm_score = pd.DataFrame(pfm_score, columns=["motif", "enhancer", "zscore"])
                pfm_score = pfm_score.set_index("motif")

                # print("\tCombine")
                pfm_score["zscoreRank"] = minmax_scale(rankdata(pfm_score["zscore"]))
                # When we built model, rank and minmax normalization was used.
                cols = ["enhancer", "zscore", "zscoreRank"]
                write_header = False
                if chunk == 0:
                    write_header = True
                pfm_score[cols].to_csv(pfmscorefile, sep="\t", header=write_header)
                # pbar.update(chunk + chunksize)

        return pfmscorefile.name

    def get_binding_score(self, pfm, peak):
        """Infer TF binding score from motif z-score and peak intensity.

        Arguments:
            pfm {[type]} -- [motif scan result]
            peak {[type]} -- [peak intensity]

        Returns:
            [type] -- [the predicted tf binding table]
        """

        # Load model
        with open(self.model, "rb") as f:
            clf = pickle.load(f)

        ft = self.filtermotifs2factors

        if self.etype == "hg38H3K27ac":
            r = pfm.merge(peak, left_on="enhancer", right_on="peak")[
                ["motif", "enhancer", "zscore", "log10_peakRPKM"]
            ]
            r = r.merge(ft, left_on="motif", right_on="Motif")
            r = r.groupby(["factor", "enhancer"])[["zscore", "log10_peakRPKM"]].mean()
            r = r.dropna().reset_index()
            table = r.compute(num_workers=self.ncore)
            table["binding"] = clf.predict_proba(table[["zscore", "log10_peakRPKM"]])[:, 1]

        elif self.etype == "p300" or self.etype == "ATAC":
            r = pfm.merge(peak, left_on="enhancer", right_on="peak")[
                ["motif", "enhancer", "zscore", "log10_peakRPKM"]
            ]
            r = r.merge(ft, left_on="motif", right_on="Motif")
            r = r.groupby(["factor", "enhancer"])[["zscore", "log10_peakRPKM"]].mean()
            r = r.dropna().reset_index()

            table = r.compute(num_workers=self.ncore)
            table["binding"] = clf.predict_proba(table[["zscore", "log10_peakRPKM"]])[:, 1]
        else:
            raise TypeError("The input enhancer data type should hg38H3K27ac, p300 or ATAC. Please provide a enhancer type with -e argument. By default is hg38H3K27ac.")

        return table

    def run_binding(self, peak_bed, outfile):

        logger.info("Peak initialization")

        filter_bed = self.set_peak_size(peak_bed)

        logger.info("Motif scan")
        pfm_weight = self.get_PWMScore(filter_bed)
        pfm = dd.read_csv(pfm_weight, sep="\t")

        logger.info("Predicting TF binding sites")
        peak_weight = self.get_peakRPKM(filter_bed)
        peak = dd.read_csv(peak_weight, sep="\t", blocksize=200e6)
        table = self.get_binding_score(pfm, peak)

        logger.info("Save results")
        table.to_csv(outfile, sep="\t", index=False)
