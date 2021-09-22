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


def read_top_interactions(fname, edges=100000):
    """Read network file and return the top interactions"""
    rnet = pd.read_csv(fname, sep="\t")  # read the GRN file
    rnet.sort_values(
        "prob", ascending=False, inplace=True
    )  # sort based on probability score
    rnet = rnet.head(edges)
    top_int = set(rnet["tf_target"])
    return top_int


def read_network(fname, edges=100000, interactions=None):
    """Read network file and return networkx DiGraph"""
    G = nx.DiGraph()  # initiate empty network
    rnet = pd.read_csv(fname, sep="\t")  # read the GRN file
    rnet.sort_values(
        "prob", ascending=False, inplace=True
    )  # sort based on probability score

    # check method of selecting edges, if a interaction set is provided, export those from the GRN
    # if not export the edges of the network.
    if interactions == None:
        rnet = rnet.head(edges)
    else:
        rnet = rnet[rnet.tf_target.isin(interactions)]

    edge_atribute_length = len(
        rnet.columns
    )  # check if network is full output or simple output
    for _, row in rnet.iterrows():
        source, target = row[0].split("_", 1)
        weight = 0 if len(row) < edge_atribute_length else float(row["prob"])
        # when the GRN with the full output is generated, more data can be exported to the diffnetwork
        # all the extra data is stored in edge attributes
        if edge_atribute_length == 6:
            tf_expression = (
                0 if len(row) < edge_atribute_length else float(row["tf_expression"])
            )
            target_expression = (
                0
                if len(row) < edge_atribute_length
                else float(row["target_expression"])
            )
            weighted_binding = (
                0 if len(row) < edge_atribute_length else float(row["weighted_binding"])
            )
            activity = 0 if len(row) < edge_atribute_length else float(row["activity"])
            try:
                G.add_edge(
                    source,
                    target,
                    weight=weight,
                    weighted_binding=weighted_binding,
                    TF_expression=tf_expression,
                    TF_activity=activity,
                    TG_expression=target_expression,
                )
            except Exception:
                print("Could not parse edge weight.")
                raise

        # if there is only prob score, only that is inc in the GRN
        else:
            try:
                G.add_edge(source, target, weight=weight, n=1)
            except Exception:
                print("Could not parse edge weight.")
                raise
    return G


def difference(GRN_source, GRN_target):
    """Calculate the network different between two GRNs.
    It first takes the nodes from both networks, and subsequently
    First adds the edges from the target network that are missing in
    the source network.
    Secondly add the edges present in both but with a higher interaction
    score in the target network
    """
    DIF = nx.create_empty_copy(GRN_source)  # copying source GRN nodes
    nodes_target = nx.create_empty_copy(GRN_target)  # copying target GRN nodes
    DIF.add_nodes_from(
        list(nodes_target.nodes)
    )  # merge the nodes if there are differences
    # (which can happen) when taking the header instead of the --unnion-interactions flagg

    # lets check if the full GRN output is loaded and if so output all atributes to the diffnetwork:
    full_output_GRN = bool(nx.get_edge_attributes(GRN_source, "TF_expression"))
    if full_output_GRN:
        logger.info("calculating diff GRN with full output")
        # lets load all the differential edges into the diffnetwork
        for (u, v, a) in GRN_target.edges(data=True):
            if (u, v) not in GRN_source.edges:
                # if the edge is not in the target network (when using head instead of union of interactions)
                # , add all source atributes to the diff network and put the target atributes to 0
                DIF.add_edge(
                    u,
                    v,
                    weight=a["weight"],
                    n=1,
                    source_weight="NA",
                    target_weight=GRN_target.edges[u, v]["weight"],
                    tf_expr_target=a["TF_expression"],
                    tf_expr_source="NA",
                    tg_expr_target=a["TG_expression"],
                    tg_expr_source="NA",
                    target_wb=a["weighted_binding"],
                    source_wb="NA",
                    TF_act_target=a["TF_activity"],
                    TF_act_source="NA",
                )
            else:  # if the edge is present in both networks and higher in the target than in the source,
                # calculate the weight difference and output all atributes
                source_weight = GRN_source.edges[u, v]["weight"]
                target_weight = GRN_target.edges[u, v]["weight"]
                wb_source = GRN_source.edges[u, v]["weighted_binding"]
                wb_target = GRN_target.edges[u, v]["weighted_binding"]

                diff_weight = target_weight - source_weight

                if (
                    diff_weight > 0
                ):  # if the interaction probability is higher in the target than in the
                    # source, add the interaction to the diffnetwork:
                    DIF.add_edge(
                        u,
                        v,
                        weight=diff_weight,
                        source_weight=source_weight,
                        target_weight=target_weight,
                        neglogweight=-np.log(diff_weight),
                        n=1,
                        tf_expr_target=a["TF_expression"],
                        tf_expr_source=GRN_source[u][v]["TF_expression"],
                        tg_expr_target=a["TG_expression"],
                        tg_expr_source=GRN_source[u][v]["TG_expression"],
                        target_wb=a["weighted_binding"],
                        source_wb=GRN_source[u][v]["weighted_binding"],
                        TF_act_target=a["TF_activity"],
                        TF_act_source=GRN_source[u][v]["TF_activity"],
                    )
    else:
        logger.info("calculating diff GRN with simple output")
        # if only the weight is loaded, lets load only that in the diffnetwork:
        for (u, v, d) in GRN_target.edges(data=True):
            if (u, v) not in GRN_source.edges:
                DIF.add_edge(u, v, weight=d["weight"], n=1)
            else:
                diff_weight = (
                    GRN_target.edges[u, v]["weight"] - GRN_source.edges[u, v]["weight"]
                )
                if diff_weight > 0:
                    DIF.add_edge(
                        u, v, weight=diff_weight, n=1, neglogweight=-np.log(diff_weight)
                    )
    return DIF


def read_expression(fname, padj_cutoff=0.05):
    """
    Read differential gene expression analysis output,
    return dictionary with namedtuples of scores, absolute fold
    change and "real" (directional) fold change.

    Parameters
    ----------
    fname: str
        DESeq2 output file.
        Tab-separated, containing (at least) 3 columns:
        1. a column with names/IDs (column name is ignored),
        2. named column "padj" (adjusted p-values)
        3. named column "log2FoldChange"
    pval = float of the cutoff of which genes will be flagged as differential

    Returns
    -------
    dict
        namedtuples of scores, absolute fold change and "real" (directional) fold change.
    """
    expression_change = dict()

    df = pd.read_table(
        fname,
        index_col=0,
        header=0,
        dtype={
            "resid": str,
            "Unnamed: 0": str,
            "log2FoldChange": float,
            "padj": float,
        },
    )[["log2FoldChange", "padj"]]
    # removes NaNs
    df.dropna(inplace=True)
    # check for duplicated index rows and return an error (but continue running)
    dup_df = df[df.index.duplicated()]
    if len(dup_df) > 0:
        dupped_gene = str(dup_df.index[0])
        logger.warning(
            f"duplicated gene names detected in differential expression file e.g.{dupped_gene}"
        )
        logger.warning(f"taking mean of values of duplicated genes")
    # average values for duplicate gene names (none hopefully)
    df = df.groupby(by=df.index, dropna=True).mean(0)

    # absolute fold change
    df["fc"] = df["log2FoldChange"].abs()

    # get the gscore (absolute fold change if significantly differential)
    df["score"] = df["fc"] * (df["padj"] < padj_cutoff)

    for k, row in df.iterrows():
        expression_change[k] = Expression(
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

    try:
        pval = mannwhitneyu(target_fc, non_target_fc)[1]
    except RecursionError:
        pval = np.NAN
        logger.warning(
            f"Could not calculate p-val (target vs non-target fold-change) for {node}."
        )
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


def plot_influence(infile, outfile):
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
        self,
        outfile,
        degenes,
        GRN_source_file=None,
        GRN_target_file=None,
        union_grns=False,
        filter=False,
        edges=100000,
        ncore=1,
        padj_cutoff=0.05,
    ):
        self.ncore = ncore
        if union_grns == False:
            logger.info(f"Reading network(s), using top {edges} edges.")
            # Load GRNs
            if GRN_source_file is None and GRN_target_file is not None:
                self.G = read_network(GRN_target_file, edges=edges)
                logger.warning("You only provide the target network!")
            elif GRN_target_file is None and GRN_source_file is not None:
                self.G = read_network(GRN_source_file, edges=edges)
                logger.warning("You only provided the source network!")
            elif GRN_target_file is None and GRN_source_file is None:
                logger.warning("You should provide at least one ANANSE network file!")
            else:
                G1_source = read_network(GRN_source_file, edges=edges)
                G2_target = read_network(GRN_target_file, edges=edges)
        if union_grns == True:
            if (GRN_source_file is None) or (GRN_target_file is None):
                logger.warning(
                    "You should provide at least two ANANSE network files to take the interaction union!"
                )
            else:
                logger.info(
                    f"Reading network(s), using the union of the top {edges} edges of each GRN."
                )
                top_int_source = read_top_interactions(GRN_source_file, edges=edges)
                top_int_target = read_top_interactions(GRN_target_file, edges=edges)
                top_int = set.union(top_int_source, top_int_target)
                G1_source = read_network(GRN_source_file, interactions=top_int)
                G2_target = read_network(GRN_target_file, interactions=top_int)
        self.G = difference(G1_source, G2_target)
        logger.info(f"Differential network has {len(self.G.edges)} edges.")
        # Load expression file
        self.expression_change = read_expression(degenes)
        self.outfile = outfile
        # Filter TFs
        self.filter = filter

    def save_reg_network(self, filename):
        """Save the network difference between two cell types to a file."""
        full_output_GRN = bool(nx.get_edge_attributes(self.G, "tf_expr_source"))
        with open(filename, "w") as nw:
            if full_output_GRN:
                logger.info("output full diff network")
                nw.write(
                    f"TF \t target \t weight \t source_weight \t target_weight \t tf_expr_source \t tf_expr_target \t tg_expr_source \t tg_expr_target \t source_wb \t target_wb \t source_Tf_act \t target_Tf_act \n"
                )
                for (u, v, d) in self.G.edges(data=True):
                    nw.write(
                        u
                        + "\t"
                        + v
                        + "\t"
                        + str(d["weight"])
                        + "\t"
                        + str(d["source_weight"])
                        + "\t"
                        + str(d["target_weight"])
                        + "\t"
                        + str(d["tf_expr_source"])
                        + "\t"
                        + str(d["tf_expr_target"])
                        + "\t"
                        + str(d["tg_expr_source"])
                        + "\t"
                        + str(d["tg_expr_target"])
                        + "\t"
                        + str(d["source_wb"])
                        + "\t"
                        + str(d["target_wb"])
                        + "\t"
                        + str(d["TF_act_source"])
                        + "\t"
                        + str(d["TF_act_target"])
                        + "\n"
                    )
            else:
                nw.write(f"TF \t target \t weight  \n")
                for (u, v, d) in self.G.edges(data=True):
                    nw.write(u + "\t" + v + "\t" + str(d["weight"]) + "\n")

    def run_target_score(self, max_degree=3):
        """Run target score for all TFs."""

        tfs = [node for node in self.G.nodes() if self.G.out_degree(node) > 0]
        logger.info(f"Differential network contains {len(tfs)} transcription factors.")

        # differentially expressed TFs
        detfs = [tf for tf in tfs if tf in self.expression_change]

        if len(detfs) == 0:
            logger.error(
                "No overlapping transcription factors found between the network file(s) "
                "(-s/--source, -t/--target) and the differential expression data (-d/--degenes)!"
            )
            sys.exit(1)

        detfs = [tf for tf in detfs if self.expression_change[tf].realfc > 0]
        if len(detfs) == 0:
            logger.error(
                "No differentially expressed TFs found with a log2 fold change above 0!"
            )
            sys.exit(1)
        else:
            logger.info(f"Out of these, {len(detfs)} are differentially expressed.")

        result = []
        if self.ncore > 1:
            try:
                pool = mp.Pool(self.ncore)
                jobs = []
                for tf in detfs:
                    jobs.append(
                        pool.apply_async(
                            targetScore,
                            (tf, self.G, self.expression_change, max_degree),
                        )
                    )
                with tqdm(total=len(jobs)) as pbar:
                    for j in jobs:
                        result.append(j.get())
                        pbar.update(1)
                pool.close()
            except Exception as e:
                if "multiprocessing" in e.__repr__():
                    logger.error(str(e))
                    logger.error("Error seems to be related to multiprocessing.")
                    logger.error(
                        "In some cases running `ananse influence` with `-n 1` will solve this issue."
                    )
                    logger.error(
                        "If it doesn't, then please file a bug report (with the output of the command run with `-n 1`) at:"
                    )
                    logger.error("https://github.com/vanheeringen-lab/ANANSE/issues")
                    sys.exit(1)
                raise
        else:
            for tf in tqdm(detfs):
                result.append(
                    targetScore(tf, self.G, self.expression_change, max_degree)
                )

        # Get results and write to file
        influence_file = open(self.outfile, "w")
        influence_file.write(
            "factor\tdirectTargets\ttotalTargets\ttargetsore\tGscore\tfactor_fc\tpval\ttarget_fc\n"
        )

        for (
            factor,
            score,
            direct_targets,
            total_targets,
            factor_fc,
            pval,
            target_fc,
        ) in result:
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

        print("\n", file=influence_file)

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
        logger.info("Saving differential network.")
        self.save_reg_network(
            ".".join(self.outfile.split(".")[:-1]) + "_diffnetwork.txt"
        )

        logger.info("Calculating target scores.")
        influence_file = self.run_target_score()

        logger.info("Calculating influence scores.")
        self.run_influence_score(influence_file, fin_expression=fin_expression)

        if plot is True:
            logger.info("Plotting results.")
            plot_influence(
                self.outfile, ".".join(self.outfile.split(".")[:-1]) + ".pdf"
            )
