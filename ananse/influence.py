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
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from scipy.stats import rankdata, mannwhitneyu
from sklearn.preprocessing import minmax_scale
from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns

from ananse import mytmpdir

warnings.filterwarnings("ignore")


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
        source, target = vals[1][0].split("_")
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
    DIF = nx.create_empty_copy(R)
    for (u, v, d) in S.edges(data=True):

        if (u, v) not in R.edges and S.edges[u, v]["weight"] > 0.5:
            DIF.add_edge(u, v, weight=d["weight"], n=1)
        elif S.edges[u, v]["weight"] - R.edges[u, v]["weight"] >= 0.3:
            DIF.add_edge(
                u, v, weight=R.edges[u, v]["weight"] - R.edges[u, v]["weight"], n=1
            )

    return DIF


def read_expression(fname):
    """Read kallisto output, return dictionary with abs fold change."""
    # Expression change
    expression_change = {"score": {}, "fc": {}, "realfc": {}}
    for line in open(fname):
        if not line.startswith("resid"):
            gene = line.split("\t")[0].strip().upper()
            foldchange = abs(float(line.split("\t")[1]))
            realFC = float(line.split("\t")[1])
            padj = float(line.split("\t")[2])
            # if padj==0:
            #     padj=1e-300
            # gscore =foldchange * (-np.log10(padj))
            if padj < 0.05:
                # gscore = np.log2(foldchange)
                gscore = foldchange
            else:
                gscore = 0
            expression_change["score"][gene] = gscore
            expression_change["fc"][gene] = foldchange
            expression_change["realfc"][gene] = realFC
    return expression_change


def influenceScore(node, G, max_degree=3, expression=None):
    """Calculate the influence score."""
    if expression is None:
        expression = {"score": {}, "fc": {}}
    total_score = 0
    # get all targets up to max_degree degree
    lengths, paths = nx.single_source_dijkstra(G, node, weight="n")
    targets = [t for t in lengths if 0 < lengths[t] <= max_degree]
    # get shortest paths based on edge weight
    lengths, paths = nx.single_source_dijkstra(G, node, weight="weight")
    # calculate influence score
    for target in targets:
        path = paths[target]
        # outdegree of parent node of the target
        # d = np.log(G.out_degree(path[-2]) + 1)
        # d = G.out_degree(path[-2])
        # expression score of the target
        # g = expression["score"].get(target, 1)
        g = expression["score"].get(target, 0)
        # Weight is cumulative product of probabilities
        # weight is a list of all path node one by one weight
        # weight = [1 - G[s][t]["weight"] for s, t in zip(path[:-1], path[1:])]
        weight = [G[s][t]["weight"] for s, t in zip(path[:-1], path[1:])]
        # cumulative sum of weight
        weight = np.cumprod(weight)[-1]
        # score = g / len(path) / d * weight
        score = g / len(path) * weight
        total_score += score
    # Get Mann-Whitney U p-value of direct targets vs. non-direct targets
    direct_targets = [n for n in G[node] if n in expression["fc"]]
    non_direct_targets = [
        n for n in list(G.nodes) if n in expression["fc"] and n not in direct_targets
    ]
    target_fc = [expression["fc"][t] for t in direct_targets]
    non_target_fc = [expression["fc"][t] for t in non_direct_targets]
    pval = mannwhitneyu(target_fc, non_target_fc)[1]
    target_fc_diff = np.mean(target_fc) - np.mean(non_target_fc)
    return (
        node,
        total_score,
        G.out_degree(node),
        len(targets),
        expression["fc"].get(node, 0),
        pval,
        target_fc_diff,
    )


def filter_TF(scores_df, network=None, tpmfile=None, tpm=20, overlap=0.98):
    tpmscore = {}
    with open(tpmfile) as tpf:
        next(tpf)
        for line in tpf:
            tpmscore[line.split()[0]] = float(line.split()[1])

    meg = lambda tf: set(network[tf]) if tf in network else set()

    tftarget = {}
    for tf in scores_df.index:
        tftarget[tf] = meg(tf)

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
    def __init__(self, Gbf=None, Gaf=None, outfile=None, expression=None, edges=10000):

        # Load GRNs
        if Gbf is None and Gaf is not None:
            self.G = read_network(Gaf, edges=edges)
            print("You only previde one network file in second cell!")
        elif Gaf is None and Gbf is not None:
            self.G = read_network(Gbf, edges=edges)
            print("You only previde one network file in first cell!")
        elif Gaf is None and Gbf is None:
            print("You should previde at list one network file!")
        else:
            G1 = read_network(Gbf, edges=edges)
            G2 = read_network(Gaf, edges=edges)
            self.G = difference(G2, G1)

        # Load expression file
        self.expression_change = read_expression(expression)

        self.outfile = outfile

    def save_reg_network(self, filename):
        with open(filename, "w") as nw:
            for (u, v, d) in self.G.edges(data=True):
                nw.write(u + "\t" + v + "\t" + str(d["weight"]) + "\n")

    def run_influence_score(self, max_degree=3):
        pool = mp.Pool()

        jobs = []
        tfs = [node for node in self.G.nodes() if self.G.out_degree(node) > 0]
        detfs = [
            g
            for g in tfs
            if g in self.expression_change["score"]
            and self.expression_change["realfc"][g] < 0
        ]
        for tf in detfs:
            jobs.append(
                pool.apply_async(
                    influenceScore, (tf, self.G, max_degree, self.expression_change)
                )
            )

        # Get results and write to file
        influence_file = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        influence_file.write(
            "factor\tdirectTargets\ttotalTargets\tinflscore\tGscore\tfactor_fc\tpval\ttarget_fc\n"
        )
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
                self.expression_change["score"][factor],
                factor_fc,
                pval,
                target_fc,
                file=influence_file,
                sep="\t",
            )

        pool.close()

        scores_df = pd.read_table(influence_file.name, index_col=0)
        scores_df["influenceScaled"] = minmax_scale(
            rankdata(scores_df["inflscore"], method="dense")
        )
        scores_df.sort_values("influenceScaled", inplace=True, ascending=False)
        # scores_df.to_csv(influence_file.name, sep='\t')

        return influence_file.name

    def rank_TF(self, influence_file, filter=None, fin_expression=None):

        scores_df1 = pd.read_table(influence_file, index_col=0)
        scores_df1["influenceScaled"] = minmax_scale(
            rankdata(scores_df1["inflscore"], method="dense")
        )
        scores_df1["GscoreScaled"] = minmax_scale(
            rankdata(scores_df1["Gscore"], method="dense")
        )
        scores_df1["sumScaled"] = minmax_scale(
            rankdata(
                scores_df1.influenceScaled + scores_df1.GscoreScaled, method="dense"
            )
        )

        scores_df1.sort_values("sumScaled", inplace=True, ascending=False)
        scores_df1 = scores_df1[
            [
                "influenceScaled",
                "GscoreScaled",
                "sumScaled",
                "directTargets",
                "inflscore",
                "factor_fc",
            ]
        ]
        scores_df1.to_csv(self.outfile, sep="\t")

        if filter:
            scores_df2 = filter_TF(
                network=self.G, scores_df=scores_df1, tpmfile=fin_expression
            )
            scores_df2.to_csv(
                ".".join(self.outfile.split(".")[:-1]) + "_filtered.txt", sep="\t"
            )

    def run_influence(self, plot=True, filter=None, fin_expression=None):
        influence_file = self.run_influence_score()
        self.rank_TF(influence_file, filter=None, fin_expression=None)

        if plot is True:
            plot_influscore(
                self.outfile, ".".join(self.outfile.split(".")[:-1]) + ".pdf"
            )
