#!/usr/bin/env python

# Copyright (c) 2009-2016 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Predict TF binding site"""

# Python imports
import os
import pickle
import warnings
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import dask.dataframe as dd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale

from pybedtools import BedTool
from genomepy import Genome
from gimmemotifs.scanner import Scanner
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import as_fasta, pfmfile_location

from ananse import mytmpdir
import ananse

warnings.filterwarnings("ignore")


class Binding(object):
    def __init__(self, genome="hg19", gene_bed=None, pfmfile=None):
        self.genome = genome
        g = Genome(self.genome)
        self.gsize = g.props["sizes"]["sizes"]

        pfmfile = pfmfile_location(pfmfile)

        # Motif information file
        if pfmfile is None:
            self.pfmfile = "../data/gimme.vertebrate.v5.1.pfm"
            self.motifs2factors = pfmfile.replace(".pfm", ".motif2factors.txt")
            self.factortable = pfmfile.replace(".pfm", ".factortable.txt")
        else:
            self.pfmfile = pfmfile
            self.motifs2factors = pfmfile.replace(".pfm", ".motif2factors.txt")
            self.factortable = pfmfile.replace(".pfm", ".factortable.txt")

        self.gene_bed = gene_bed

        package_dir = os.path.dirname(ananse.__file__)
        self.model = os.path.join(package_dir, "db", "dream_model.txt")

        # dream_model.txt is the logistic regression model.

    def set_peak_size(self, peaks, seqlen=200):

        gsizedic = {}
        with open(self.gsize) as gsizefile:
            for chrom in gsizefile:
                gsizedic[chrom.split()[0]] = int(chrom.split()[1])

        s = ""
        for peak in peaks:

            if peak.length < seqlen or peak.length > seqlen:
                # get the summit and the flanking low and high sequences
                summit = (peak.start + peak.end) // 2
                start, end = summit - seqlen // 2, summit + seqlen // 2
            else:
                start, end = peak.start, peak.start

            # remove seq which langer than chromosome length or smaller than 0
            if start > 0 and end < gsizedic[peak.chrom]:
                s += (
                    str(peak.chrom)
                    + "\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t"
                    + str(peak.fields[-1])
                    + "\n"
                )

        npeaks = BedTool(s, from_string=True)

        return npeaks

    def clear_peak(self, peak_bed, filter_promoter=True, up=2000, down=2000):
        """
        Filter the enhancer peaks in promoter range.
        """
        # set all seq to 200bp
        peaks = BedTool(peak_bed)
        peaks = self.set_peak_size(peaks, 200)

        # remove all peaks that overlap with TSS(up2000 to down2000).
        b = BedTool(self.gene_bed)
        b = b.flank(l=1, r=0, s=True, g=self.gsize).slop(  # noqa: E741
            l=up, r=down, g=self.gsize, s=True  # noqa: E741
        )
        vals = []
        # for f in b.intersect(peaks, wo=True, nonamecheck=True):
        # Bedtools don't have nonamecheck option now?
        for f in b.intersect(peaks, wo=True):
            chrom = f[0]
            peak_start, peak_end = int(f[13]), int(f[14])
            vals.append(chrom + ":" + str(peak_start) + "-" + str(peak_end))
        fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        with open(peak_bed) as pbed:
            for line in pbed:
                if filter_promoter:
                    if (
                        line.split()[0] + ":" + line.split()[1] + "-" + line.split()[2]
                        not in vals
                    ):
                        fl2.write(line)
                else:
                    fl2.write(line)
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
        """
        Scan motif in every peak.
        """
        pfmscorefile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        seqs = [
            s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=self.genome).ids
        ]
        s = Scanner()
        s.set_motifs(self.pfmfile)
        s.set_threshold(threshold=0.0)
        s.set_genome(self.genome)
        with open(self.pfmfile) as f:
            motifs = read_motifs(f)

        chunksize = 10000
        # Run 10k peaks one time.
        for chunk in range(0, len(seqs), chunksize):
            chunk_seqs = seqs[chunk : chunk + chunksize]
            print("\t", chunk, "-", chunk + chunksize, "enhancers")
            pfm_score = []
            it = s.best_score(chunk_seqs, zscore=True, gc=True)
            # We are using GC-normalization for motif scan because many sequence is GC-enriched.
            # GimmeMotif develop branch already include GC-normalization option now.
            for seq, scores in zip(chunk_seqs, it):
                for motif, score in zip(motifs, scores):
                    pfm_score.append([motif.id, seq, score])
            pfm_score = pd.DataFrame(pfm_score, columns=["motif", "enhancer", "zscore"])
            pfm_score = pfm_score.set_index("motif")

            print("\tCombine")
            pfm_score["zscoreRank"] = minmax_scale(rankdata(pfm_score["zscore"]))
            # When we built model, rank and minmax normalization was used.
            cols = ["enhancer", "zscore", "zscoreRank"]
            write_header = False
            if chunk == 0:
                write_header = True
            pfm_score[cols].to_csv(pfmscorefile, sep="\t", header=write_header)
        return pfmscorefile.name

    def get_binding_score(self, pfm, peak):
        """
        Infer TF binding score from motif z-score and peak intensity.
        """
        # Load model
        with open(self.model, "rb") as f:
            clf = pickle.load(f)

        ft = dd.read_csv(self.factortable, sep="\t")
        r = pfm.merge(peak, left_on="enhancer", right_on="peak")[
            ["motif", "enhancer", "zscore", "peakRPKMScale"]
        ]
        r = r.merge(ft, left_on="motif", right_on="motif")
        r = r.groupby(["factor", "enhancer"])[["zscore", "peakRPKMScale"]].max()
        r = r.dropna().reset_index()

        table = r.compute()
        print("Predicting TF binding sites")
        table["binding"] = clf.predict_proba(table[["zscore", "peakRPKMScale"]])[:, 1]
        print("Save results")
        return table

    def run_binding(self, peak_bed, outfile):
        filter_bed = self.clear_peak(peak_bed)
        print("Motif scanning")
        pfm_weight = self.get_PWMScore(filter_bed)
        pfm = dd.read_csv(pfm_weight, sep="\t")
        peak_weight = self.get_peakRPKM(filter_bed)
        peak = dd.read_csv(peak_weight, sep="\t")
        table = self.get_binding_score(pfm, peak)
        table.to_csv(outfile, sep="\t", index=False)
