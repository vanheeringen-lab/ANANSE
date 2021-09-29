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
from matplotlib.lines import Line2D


warnings.filterwarnings("ignore")

# Here because of multiprocessing and pickling
Expression = namedtuple("Expression", ["score", "absfc", "realfc"])


def read_top_interactions(fname, edges=100000, GRNsort_column="prob"):
    """Read network file and return the top interactions
    by default it incorporates all the data (prob score) to select the most
    certain interactions, however other sorting options such as the weighted binding score"""
    rnet = pd.read_csv(
        fname, usecols=["tf_target", GRNsort_column], sep="\t"
    )  # read the GRN file
    rnet.sort_values(
        GRNsort_column, ascending=False, inplace=True
    )  # sort based on the selected value (default probability) score
    rnet = rnet.head(edges)
    top_int = set(rnet["tf_target"])
    return top_int


def read_network(
    fname, edges=100000, interactions=None, GRNsort_column="prob", full_output=False
):
    """Read network file and return networkx DiGraph"""

    # read GRN files
    if full_output == False:
        data_columns=["tf_target", "prob"]
    if full_output:
        data_columns=[
            "tf_target",
            "prob",
            "tf_expression",
            "target_expression",
            "weighted_binding",
            "activity",
        ]
    rnet = pd.read_csv(
        fname,
        sep="\t",
        usecols=data_columns,
        dtype="float64",
        converters={"tf_target": str},
    )  # read the GRN file
    # sort on selection variable
    rnet.sort_values(GRNsort_column, ascending=False, inplace=True)
    if interactions == None:
        rnet = rnet.head(edges)  # no interactions so take the top head of interactions
    else:
        rnet = rnet[rnet.tf_target.isin(interactions)]

    G = nx.DiGraph()  # initiate empty network
    for _, row in rnet.iterrows():
        source, target = row[0].split("_", 1)
        if full_output == False:
            try:
                G.add_edge(source, target, weight=row["prob"], n=1)
            except Exception:
                logger.error(f"Could not parse edge weight of edge {source}:{target}")
                raise
        if full_output:
            try:
                G.add_edge(
                    source,
                    target,
                    weight=row["prob"],
                    weighted_binding=row["weighted_binding"],
                    tf_expression=row["tf_expression"],
                    tg_expression=row["target_expression"],
                    tf_activity=row["activity"],
                )
            except Exception:
                logger.error(f"Could not parse edge weight {source}:{target}")
                raise
    return G


def difference(GRN_source, GRN_target, full_output=False):
    """Calculate the network different between two GRNs.
    It first takes the nodes from both networks, and subsequently
    first adds the edges from the target network that are missing in
    the source network.
    Secondly it add the edges present in both but with a higher interaction
    score in the target network
    """
    DIF = nx.create_empty_copy(GRN_source)  # copying source GRN nodes
    nodes_target = nx.create_empty_copy(GRN_target)  # copying target GRN nodes
    DIF.add_nodes_from(
        list(nodes_target.nodes)
    )  # merge the nodes if there are differences
    # (which can happen) when taking the header instead of the --unnion-interactions flagg

    # lets check if the full GRN output is loaded and if so output all atributes to the diffnetwork:
    if full_output:
        logger.info("calculating diff GRN with full output")
        # lets load all the  edges into the diffnetwork
        for (u, v, a) in GRN_target.edges(
            data=True
        ):  # u = source node, v = target node, a = dictionary ofedge atributes
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
                    tf_expr_target=a["tf_expression"],
                    tf_expr_source=GRN_source[u][v]["tf_expression"],
                    tg_expr_target=a["tg_expression"],
                    tg_expr_source=GRN_source[u][v]["tg_expression"],
                    target_wb=a["weighted_binding"],
                    source_wb=GRN_source[u][v]["weighted_binding"],
                    TF_act_target=a["tf_activity"],
                    TF_act_source=GRN_source[u][v]["tf_activity"],
                )
    else:
        logger.info("calculating diff GRN")
        # if only the weight is loaded, lets load only that in the diffnetwork:
        for (u, v, d) in GRN_target.edges(data=True):
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
            f"Could not calculate p-val (target vs non-target fold-change) for {node}, targets = {len(target_fc)}, non-target = {len(non_target_fc)}."
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
    plt.clf()


class Influence(object):
    def __init__(
        self,
        outfile,
        degenes,
        GRN_source_file=None,
        GRN_target_file=None,
        filter=False,
        edges=100000,
        ncore=1,
        GRNsort_column="prob",
        padj_cutoff=0.05,
        full_output=False,
        GRN_wb = False,
    ):
        self.ncore = ncore
        self.full_output = full_output
        self.GRNsort_column = GRNsort_column
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
            logger.info(
                f"Reading network(s), using the union of the top {edges} edges of each GRN."
            )
            top_int_source = read_top_interactions(
                GRN_source_file, edges=edges, GRNsort_column=GRNsort_column
            )
            top_int_target = read_top_interactions(
                GRN_target_file, edges=edges, GRNsort_column=GRNsort_column
            )
            top_int = set.union(top_int_source, top_int_target)
            G1_source = read_network(
                GRN_source_file,
                interactions=top_int,
                GRNsort_column=GRNsort_column,
                full_output=full_output,
            )
            G2_target = read_network(
                GRN_target_file,
                interactions=top_int,
                GRNsort_column=GRNsort_column,
                full_output=full_output,
            )
            self.G = difference(G1_source, G2_target, full_output=full_output)
            logger.info(f"Differential network has {len(self.G.edges)} edges.")
        # Load expression file
        self.expression_change = read_expression(degenes)
        self.outfile = outfile
        # Filter TFs
        self.filter = filter

    def save_reg_network(self, filename, full_output=False):
        """Save the network difference between two cell types to a file."""
        with open(filename, "w") as nw:
            if full_output == False:
                nw.write("TF\ttarget\tweight\n")
                for (u, v, d) in self.G.edges(data=True):
                    nw.write(u + "\t" + v + "\t" + str(d["weight"]) + "\n")
            if full_output:
                logger.info("output full diff network")
                nw.write(
                    "tf\ttarget\tweight\tweight_source\tweight_target\ttf_expr_source\ttf_expr_target\ttg_expr_source\ttg_expr_target\twb_source\twb_target\tsource_tf_act\ttarget_tf_act\n"
                )
                for (u, v, d) in self.G.edges(data=True):
                    cols = [
                        u,
                        v,
                        d["weight"],
                        d["source_weight"],
                        d["target_weight"],
                        d["tf_expr_source"],
                        d["tf_expr_target"],
                        d["tg_expr_source"],
                        d["tg_expr_target"],
                        d["source_wb"],
                        d["target_wb"],
                        d["TF_act_source"],
                        d["TF_act_target"],
                    ]
                    nw.write("\t".join(str(v) for v in cols) + "\n")
#def plot_TF_GRN(self, infile, outfile, GRN_wb = False, network_algorithm = 'neato', n_TFs = 20, cmap = 'viridis', full_output = False):
    def plot_TF_GRN(self, infile, outfile, GRN_wb = False, network_algorithm = 'neato', n_TFs = 20, cmap = 'viridis', full_output = False, cutoff_val = 0.1):
        """
        Plot the top20 of differential TFs and their interactions with the highest interaction score into a GRN network image

        Parameters
        ----------
        self: influence_object to extract the diff network self.G, #future optimization = read diffnetwork instead
        infile: influence output file
        outfile: output file location
        GRN_wb: when true and a -full--output influence GRN is suplied, use weighted binding interaction score to vizualize
                edges instead of the probability score
        network_algorithm: pyviz cluster algorithm used for node placement, options include: neato, dot, fdp, twopi, sfdp, circo 
        n_TFs: number of (significantly differential expressed) Tfs to vizualize
        cmap: matplotlib colour pallet used
        full_output: if the GRN is a full output GRN (usefull for ploting weighted binding instead of prob score)
        cutoff_val: minimum value the prob/wb score needs to have to be included in the GRN, 0.1 seems like a good number.
    """
        
        cmap=plt.get_cmap(cmap)
        mogrify = pd.read_table(infile, index_col="factor")
        mogrify = mogrify[mogrify.GscoreScaled > 0] #plot only TFs that are differential 
        factors = list(mogrify.sort_values("sumScaled").tail(n_TFs).index)#get the top TFs to include in the GRN    
        #get the mean interaction score of the selected top TFs
        def filter_node(n1):
            return n1 in  factors
        TF_G = nx.subgraph_view(self.G, filter_node=filter_node)#filter the diffnetwork to only contain topTF-topTF
        #filter the network to only contain interactions above the cutoff value

        if not full_output: 
            if GRN_wb:
                logger.error(f"weighted binding edgescores require a --full-output GRN generated in both the network and influence step, using interaction score instead ")
            edge_info = 'weight'
            edge_info_title = 'interaction score'
            TF_G_large_int = nx.DiGraph(((u, v, e) for u,v,e in TF_G.edges(data=True) if e[edge_info] > cutoff_val))
            TF_G2 = TF_G_large_int
            edge_atribute = list(nx.get_edge_attributes(TF_G2,'weight').values())
        if full_output:
            if GRN_wb:
                edge_info = 'wb'
                edge_info_title = 'weighted binding score'
                logger.info(f"using weighted binding edge scores for GRN edges")
            else:
                edge_info = 'weight'
                edge_info_title = 'interaction score'
                logger.info(f"using interaction score edge scores for GRN edges")
            #sort all interactions on interaction or weighted binding score and filter away interactions below the cutoff:
            TF_G_large_int = nx.DiGraph(((u, v, e) for u,v,e in TF_G.edges(data=True) if e[edge_info] > cutoff_val))#does the interaction exceed the cutoff value
            TF_G2 = nx.DiGraph(((u, v, e) for u,v,e in TF_G_large_int.edges(data=True) if e[f'target_{edge_info}'] > e[f'source_{edge_info}'])) #is the interaction higher in target then source
            target_wb=nx.get_edge_attributes(TF_G2,f'target_{edge_info}')
            source_wb=nx.get_edge_attributes(TF_G2,f'source_{edge_info}')
            edge_atribute = {key: target_wb[key] - source_wb.get(key, 0) for key in target_wb}#calculate the source-target difference
            edge_atribute = list(edge_atribute.values())
        TF_G2 = TF_G2.to_directed()#make the network directed agian
        TF_G2.remove_nodes_from(list(nx.isolates(TF_G2)))#remove TFs with no interactions
        edge_atribute_scaled = minmax_scale(edge_atribute, feature_range=(0, 1), axis=0, copy=True) #scale all values
        #calculate quantile regions to dipslay 4 relevant numbers within the edge legend
        #normzalize the edge atributes (interaction scores) to be between 0 and 1 (where 0-1 corresponds to edge width)
        edges_norm_weight = [0,
                         round((np.quantile(sorted(edge_atribute_scaled), 1)/4),3),
                         round((np.quantile(sorted(edge_atribute_scaled), 1)/2),3),
                         round(np.quantile(sorted(edge_atribute_scaled), 1),3)]
        #Also get the numbers of the unnormalized edge numbers to display within the legend next to the lines
        edges_weight = [0,
                         round((np.quantile(sorted(edge_atribute), 1)/4),3),
                         round((np.quantile(sorted(edge_atribute), 1)/2),3),
                         round(np.quantile(sorted(edge_atribute), 1),3)]
        #lets calculate the nodes their outdegree (edges regulating other TFs):
        if edge_info == 'weight':
            outdegree = pd.DataFrame(TF_G2.out_degree(weight = edge_info))
            outdegree = outdegree[1]
        else: #calculate the source target difference, add this as node atribute and then calculate the outdegree
            source_edgeweight=nx.get_edge_attributes(TF_G2,f'source_{edge_info}')
            target_edgeweight=nx.get_edge_attributes(TF_G2,f'target_{edge_info}')
            diff_edgeweight = {key: target_edgeweight[key] - source_edgeweight.get(key, 0) for key in target_edgeweight}
            nx.set_edge_attributes(TF_G2, diff_edgeweight, name = edge_info)    
            outdegree = pd.DataFrame(TF_G2.out_degree(weight = edge_info))
            outdegree = outdegree[1]
        node_outdegree_size = 600 + outdegree*100
        colors=outdegree
        vmin = min(colors)
        vmax = max(colors)
        plt.figure(figsize=(10, 10))
        #calculate node position of the graph
        pos = nx.drawing.nx_agraph.graphviz_layout(TF_G2, prog=network_algorithm)
        #plot the TF nodes
        nx.draw_networkx_nodes(TF_G2,
                               pos,
                               node_size = node_outdegree_size,
                               node_color = outdegree, 
                               cmap=cmap)
        #plot the TF name labels:
        nx.draw_networkx_labels(TF_G2, pos,font_color ='black', bbox = dict(facecolor='white', alpha=0.5, pad=0))
        #plot the node colour legend:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
        cbar = plt.colorbar(sm)
        cbar.ax.set_ylabel('outdegree (regulation other TFs)', rotation=270, labelpad = 25)
        #plot the TF-TF edges:
        nx.draw_networkx_edges(TF_G2, 
                               pos, 
                               arrows=True, 
                               arrowstyle = '->',
                               arrowsize = 20,
                               width = edge_atribute_scaled,
                               node_size = node_outdegree_size,
                               connectionstyle='arc3, rad = 0.1')
        #add edge width legend:
        lines = []
        for i, width in enumerate(edges_norm_weight):
            lines.append(Line2D([],[], linewidth=width, color='black'))
        legend2 = plt.legend(lines, edges_weight, bbox_to_anchor=(0, 0.5), frameon=False, title = f'{edge_info_title}') 
        #save the plot
        plt.savefig(outfile,dpi=300)
        plt.clf()

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
            (".".join(self.outfile.split(".")[:-1]) + "_diffnetwork.txt"),
            full_output=self.full_output,
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
            logger.info("Plotting top TF GRN.")
            self.plot_TF_GRN(infile = self.outfile, outfile = (".".join(self.outfile.split(".")[:-1]) + "_topTF_GRN.png"))
