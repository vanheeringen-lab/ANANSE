#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

"""Predict TF influence score"""

# Python imports
from __future__ import print_function
import sys
import warnings
from collections import namedtuple
from loguru import logger
from tqdm import tqdm

import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from sklearn.preprocessing import minmax_scale
from scipy.stats import rankdata, mannwhitneyu
from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns


warnings.filterwarnings("ignore")

# Here because of multiprocessing and pickling
Expression = namedtuple("Expression", ["score", "absfc", "realfc"])


def read_network(fname, edges=100000):
    """Read network file and return networkx DiGraph."""

    G = nx.DiGraph()

    rnet = pd.read_csv(fname, sep="\t")
    nrnet = rnet.sort_values("prob", ascending=False)
    if len(nrnet) < edges:
        usenet = nrnet
    else:
        usenet = nrnet[:edges]

    for vals in usenet.iterrows():
        source, target = vals[1][0].split("_", 1)
        try:
            if len(vals[1]) > 1:
                # weight = 1 - float(vals[1])
                weight = float(vals[1][1])
                # if weight < 0 or weight > 1:
                #    sys.stderr.write("expect weight between 0 and 1")
                #    sys.exit(1)
            else:
                weight = 0
            G.add_edge(source, target, weight=weight, n=1)
        except Exception:
            sys.stderr.write("could not parse edge weight\n")
            raise
    return G


def difference(S, R):
    """Calculate the network different between two cell types."""
    DIF = nx.create_empty_copy(R)
    for (u, v, d) in S.edges(data=True):
        if (u, v) not in R.edges:
            if S.edges[u, v]["weight"] > 0.5:
                DIF.add_edge(u, v, weight=d["weight"], n=1)
        else:
            diff_weight = S.edges[u, v]["weight"] - R.edges[u, v]["weight"]
            if diff_weight >= 0.3:
                DIF.add_edge(
                    u, v, weight=diff_weight, n=1, neglogweight=-np.log(diff_weight)
                )
    return DIF


def read_expression(fname):
    """Read differential gene expression analysis output, return dictionary with namedtuples of scores, absolute fold
    change and "real" (directional) fold change.

    input:
    a tab-separated file containing 3 columns (HGNC gene symbols, (adjusted) p-values and log2foldchange)
    header is omitted if starting with "resid"
    """
    expression_change = dict()

    df = pd.read_table(
        fname,
        index_col=0,
        header=0,
        dtype={"resid": str, "log2FoldChange": float, "padj": float},
    )

    # convert to upper case (todo: this is not strictly necessary)
    df.index = [index.upper() for index in df.index]

    # absolute fold change
    df["fc"] = df["log2FoldChange"].abs()

    # get the gscore (absolute fold change if significanlty differential)
    df["score"] = df["fc"] * (df["padj"] < 0.05)

    for k, row in df.iterrows():
        expression_change[row.name] = Expression(
            score=row.score, absfc=row.fc, realfc=row.log2FoldChange
        )

    return expression_change


def targetScore(node, G, expression_change, max_degree=3):
    """Calculate the influence score."""

    # debug only.
    # todo
    # if expression_change is None:
    #     expression_change = {"score": {}, "fc": {}}

    total_score = 0

    # Get the targets that are within a certain number of steps from TF
    lengths, paths = nx.single_source_dijkstra(G, node, cutoff=max_degree - 1)
    targets = [t for t in lengths if 0 < lengths[t] <= max_degree]

    for target in paths:
        all_paths = {}
        # Calculate all paths from TF to target to select to path with the lowest total weight
        for path in nx.all_simple_paths(G, node, target, cutoff=max_degree - 1):
            if len(path) <= max_degree:
                weight = np.cumprod(
                    [G[s][t]["weight"] for s, t in zip(path, path[1:])]
                )[-1]
                # Add weight, corrected for the length of the path
                all_paths[tuple(path)] = weight / (len(path) - 1)
        if len(all_paths) > 0:
            path, weight = sorted(all_paths.items(), key=lambda p: p[1])[-1]

            # print(target, path, weight)

            # outdegree of parent node of the target
            # d = np.log(G.out_degree(path[-2]) + 1)
            # d = G.out_degree(path[-2])

            # the level (or the number of steps) that gene is away from transcription factor
            pathlen = len(path)

            # expression score of the target
            g = expression_change[target].score if target in expression_change else 0

            # weight is cumulative product of probabilities
            # weight = [G[s][t]["weight"] for s, t in zip(path[:-1], path[1:])]

            # cumulative sum of weight
            # weight = np.cumprod(weight)[-1]

            # score = g / len(path) / d * weight
            score = g / pathlen * weight
            total_score += score

    # Get Mann-Whitney U p-value of direct targets vs. non-direct targets
    direct_targets = [n for n in G[node] if n in expression_change]
    non_direct_targets = [
        n for n in list(G.nodes) if n in expression_change and n not in direct_targets
    ]

    target_fc = [expression_change[t].absfc for t in direct_targets]
    non_target_fc = [expression_change[t].absfc for t in non_direct_targets]

    pval = mannwhitneyu(target_fc, non_target_fc)[1]
    target_fc_diff = np.mean(target_fc) - np.mean(non_target_fc)

    # factor, targetScore, directTargets, totalTargets, Gscore, pval, target_fc
    return (
        node,
        total_score,
        G.out_degree(node),
        len(targets),
        expression_change[node].absfc if node in expression_change else 0,
        pval,
        target_fc_diff,
    )


def filter_TF(scores_df, network=None, tpmfile=None, tpm=20, overlap=0.98):
    """Filter TFs:
    1) it have high expression in origin cell type;
    2) 98% of its target genes are also regulated by previous TFs.
    """

    tpmscore = {}
    with open(tpmfile) as tpf:
        next(tpf)
        for line in tpf:
            tpmscore[line.split()[0]] = float(line.split()[1])

    tftarget = {}
    for tf in scores_df.index:
        tftarget[tf] = set(network[tf]) if tf in network else set()

    ltf = list(scores_df.index)

    keeptf = []
    for i in ltf:
        passtf = []
        if len(tftarget[i]) > 0:
            for j in ltf[: ltf.index(i)]:
                if len(tftarget[i] & tftarget[j]) / len(tftarget[i]) > overlap:
                    break
                else:
                    passtf.append(j)
            if passtf == ltf[: ltf.index(i)] and i in tpmscore and tpmscore[i] < tpm:
                keeptf.append(i)
    scores_df = scores_df.loc[keeptf]
    scores_df.sort_values("sumScaled", inplace=True, ascending=False)
    return scores_df


def plot_influscore(infile, outfile):
    """Plot TF influence score to expression."""

    mogrify = pd.read_table(infile, index_col="factor")
    mogrify = mogrify.dropna()
    factors = list(mogrify.sort_values("sumScaled").tail(20).index)
    # factors = list(mogrify.sort_values("sumScaled").tail(20).index)
    xcol = "factor_fc"
    plt.figure(figsize=(8, 6))
    sns.regplot(
        data=mogrify,
        x=xcol,
        y="sumScaled",
        fit_reg=False,
        scatter_kws={"s": mogrify["directTargets"] / 10, "alpha": 0.5},
    )
    x = mogrify.loc[factors, xcol]
    y = mogrify.loc[factors, "sumScaled"]
    texts = []
    for s, xt, yt in zip(factors, x, y):
        texts.append(plt.text(xt, yt, s))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black"))
    plt.xlabel("Log2 fold change of TF")
    plt.ylabel("Influence score")
    plt.savefig(outfile, dpi=300)


class Influence(object):
    def __init__(
        self, outfile, degenes, Gbf=None, Gaf=None, filter=False, edges=100000, ncore=1
    ):

        self.ncore = ncore
        logger.info("Reading network(s)")
        # Load GRNs
        if Gbf is None and Gaf is not None:
            self.G = read_network(Gaf, edges=edges)
            logger.warning("You only provide the target network!")
        elif Gaf is None and Gbf is not None:
            self.G = read_network(Gbf, edges=edges)
            logger.warning("You only provided the source network!")
        elif Gaf is None and Gbf is None:
            logger.warning("You should provide at least one ANANSE network file!")
        else:
            G1 = read_network(Gbf, edges=edges)
            G2 = read_network(Gaf, edges=edges)
            self.G = difference(G2, G1)

        # Load expression file
        self.expression_change = read_expression(degenes)

        self.outfile = outfile

        # Filter TFs
        self.filter = filter

    def save_reg_network(self, filename):
        """Save the network difference between two cell types to a file."""

        with open(filename, "w") as nw:
            for (u, v, d) in self.G.edges(data=True):
                nw.write(u + "\t" + v + "\t" + str(d["weight"]) + "\n")

    def run_target_score(self, max_degree=3):
        """Run target score for all TFs."""

        pool = mp.Pool(self.ncore)
        jobs = []

        tfs = [node for node in self.G.nodes() if self.G.out_degree(node) > 0]

        # differentially expressed TFs
        detfs = [tf for tf in tfs if tf in self.expression_change]
        if len(detfs) == 0:
            sys.stderr.write(
                "no overlapping transcription factors found between the network file(s) "
                "(-s/--source, -t/--target) and the differential expression data (-d/--degenes)\n"
            )
            sys.exit(1)

        detfs = [tf for tf in detfs if self.expression_change[tf].realfc > 0]
        if len(detfs) == 0:
            sys.stderr.write(
                "no differentially expressed TFs found with a log2 fold change above 0\n"
            )
            sys.exit(1)

        for tf in detfs:
            jobs.append(
                pool.apply_async(
                    targetScore, (tf, self.G, self.expression_change, max_degree)
                )
            )

        # Get results and write to file
        influence_file = open(self.outfile, "w")
        influence_file.write(
            "factor\tdirectTargets\ttotalTargets\ttargetsore\tGscore\tfactor_fc\tpval\ttarget_fc\n"
        )

        with tqdm(total=len(jobs)) as pbar:
            for j in jobs:
                (
                    factor,
                    score,
                    direct_targets,
                    total_targets,
                    factor_fc,
                    pval,
                    target_fc,
                ) = j.get()
                print(
                    factor,
                    direct_targets,
                    total_targets,
                    score,
                    self.expression_change[factor].score,
                    factor_fc,
                    pval,
                    target_fc,
                    file=influence_file,
                    sep="\t",
                )
                pbar.update(1)
        print("\n", file=influence_file)

        pool.close()
        influence_file.close()

        scores_df = pd.read_table(self.outfile, index_col=0)
        scores_df["targetScaled"] = minmax_scale(
            rankdata(scores_df["targetsore"], method="dense")
        )
        scores_df.sort_values("targetScaled", inplace=True, ascending=False)

        return self.outfile

    def run_influence_score(self, influence_file, fin_expression=None):
        """Calculate influence score from target score and gscore"""

        scores_df = pd.read_table(influence_file, index_col=0)

        scores_df["targetScaled"] = minmax_scale(
            rankdata(scores_df["targetsore"], method="dense")
        )
        scores_df["GscoreScaled"] = minmax_scale(
            rankdata(scores_df["Gscore"], method="dense")
        )
        scores_df["sumScaled"] = minmax_scale(
            rankdata(scores_df.targetScaled + scores_df.GscoreScaled, method="dense")
        )

        scores_df.sort_values("sumScaled", inplace=True, ascending=False)
        scores_df = scores_df[
            [
                "targetScaled",
                "GscoreScaled",
                "sumScaled",
                "directTargets",
                "targetsore",
                "factor_fc",
            ]
        ]

        scores_df.to_csv(self.outfile, sep="\t")

        if self.filter:
            scores_df2 = filter_TF(
                network=self.G, scores_df=scores_df, tpmfile=fin_expression
            )
            scores_df2.to_csv(
                ".".join(self.outfile.split(".")[:-1]) + "_filtered.txt", sep="\t"
            )

    def run_influence(self, plot=True, fin_expression=None):

        logger.info("Run target score")
        influence_file = self.run_target_score()

        logger.info("Run influence score")
        self.run_influence_score(influence_file, fin_expression=fin_expression)

        logger.info("Save results")
        self.save_reg_network(
            ".".join(self.outfile.split(".")[:-1]) + "_diffnetwork.txt"
        )

        if plot is True:
            logger.info("Plot results")
            plot_influscore(
                self.outfile, ".".join(self.outfile.split(".")[:-1]) + ".pdf"
            )
