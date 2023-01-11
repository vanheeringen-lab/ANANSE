"""Predict TF influence score"""
from heapq import heappush, heappop
from itertools import count
import os
import shutil
import sys
import warnings
import genomepy
from typing import Union
from collections import namedtuple
from loguru import logger
from tqdm.auto import tqdm
import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from sklearn.preprocessing import minmax_scale
from scipy.stats import rankdata, mannwhitneyu

from . import SEPARATOR
from .utils import mytmpdir


warnings.filterwarnings("ignore")

# Here because of multiprocessing and pickling
Expression = namedtuple("Expression", ["score", "absfc", "realfc"])

# dictionary used to convert column names from ANANSE network
# to column names for ANANSE influence
GRN_COLUMNS = {
    "prob": "weight",
    # "weighted_binding": "weighted_binding",
    # "tf_expression": "tf_expression",
    # "target_expression": "tg_expression",
    "activity": "tf_activity",
}

# This piece of code is adapted from the networkx code licensed under a 3-clause license:
# https://networkx.org/documentation/networkx-2.7/#license
def dijkstra_prob_length(G, source, weight, cutoff=None, target=None):  # noqa
    """Uses Dijkstra's algorithm to find shortest weighted paths

    Parameters
    ----------
    G : NetworkX graph
    source : node
        Starting node for paths.
    weight: string
        Name of attribute that contains weight (as a probability between 0 and 1,
        where 0 is the minimum and 1 the maximum).
    target : node label, optional
        Ending node for path. Search is halted when target is found.
    cutoff : integer or float, optional
        Minimum combined, weighted probability.
        If cutoff is provided, only return paths with summed weight >= cutoff.

    Returns
    -------
    paths, distance : dictionaries
        Dictionary of shortest paths keyed by target, and
        a mapping from node to shortest distance to that node from one
        of the source nodes.

    Raises
    ------
    NodeNotFound
        If the `source` is not in `G`.
    """
    # cumulative weight function for values < 1
    weight_function = lambda u, v, data: -np.log10(data[weight])  # noqa
    paths = {source: [source]}

    if cutoff == 0:
        cutoff = None
    if cutoff is not None:
        if cutoff < 0 or cutoff > 1:
            raise ValueError(
                "cutoff (combined, weighted probability) should be between 0 and 1"
            )
        cutoff = -np.log10(cutoff)

    G_succ = G._succ if G.is_directed() else G._adj  # noqa

    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {}
    # fringe is heapq with 3-tuples (distance,c,node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    fringe = []

    seen[source] = 0
    push(fringe, (0, next(c), source))
    while fringe:
        (d, _, v) = pop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        if v == target:
            break
        for u, e in G_succ[v].items():
            cost = weight_function(v, u, e)
            if cost is None:
                continue

            # This is the major difference, both probability and length are taken
            # into account
            length = len(paths[v])
            vu_dist = dist[v] + cost - np.log10(1 / length)  # noqa
            if length > 1:
                vu_dist = vu_dist + np.log10(1 / (length - 1))

            if cutoff is not None:
                if vu_dist > cutoff:
                    continue
            if u in dist:
                u_dist = dist[u]
                if vu_dist < u_dist:
                    raise ValueError("Contradictory paths found:", "negative weights?")
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths[u] = paths[v] + [u]  # noqa

    paths = {k: v for k, v in paths.items() if len(v) > 1}
    dist = {k: 10**-v for k, v in dist.items() if k in paths}
    return paths, dist


def read_network(fname, full_output=False):
    """
    Read a network file and return a DataFrame.

    Subsets the graph to interactions if given, else to the number top interactions given by edges.
    If both interactions and edges are none, return the whole graph.

    Peak memory usage is about ~5GB per million edges (tested with 0.1m, 1m and 10m edges)
    """
    # read GRN files
    data_columns = ["tf_target", "prob"]
    if full_output:
        data_columns = [
            "tf_target",
            "prob",
            "tf_expression",
            "target_expression",
            "weighted_binding",
            "activity",
        ]
    # read the GRN file
    rnet = pd.read_csv(
        fname,
        sep="\t",
        usecols=data_columns,
        dtype="float64",
        converters={"tf_target": str},
        index_col="tf_target",
    )

    return rnet


def read_network_to_graph(
    fname,
    edges: Union[int, None] = 100_000,
    interactions=None,
    sort_by="prob",
    full_output=False,
):
    """
    Read a network file and return a networkx DiGraph.

    Subsets the graph to interactions if given, else to the number top interactions given by edges.
    If both interactions and edges are none, return the whole graph.

    Peak memory usage is about ~5GB per million edges (tested with 0.1m, 1m and 10m edges)
    """
    rnet = read_network(fname, full_output)
    if interactions is not None:
        rnet = rnet[rnet.index.isin(interactions)]
    elif edges is not None:
        rnet = rnet.sort_values(sort_by).tail(edges)

    # split the transcription factor and target gene into 2 columns
    rnet[["source", "target"]] = rnet.index.to_series().str.split(
        SEPARATOR, expand=True
    )
    rnet.reset_index(drop=True, inplace=True)

    # rename the columns
    rnet.rename(columns=GRN_COLUMNS, inplace=True)

    # load into a network with TFs and TGs as nodes, and the interaction scores as edges
    grn = nx.from_pandas_edgelist(rnet, edge_attr=True, create_using=nx.DiGraph)

    return grn


def difference(
    grn_source,
    grn_target,
    outfile=None,
    edges=100_000,
    sort_by="prob",
    full_output=False,
    select_after_join=False,
):
    """
    Calculate the network different between two GRNs.

    First take the nodes from both networks, and add
    edges from the target network that are missing in the source network.
    Then add edges present in both but with a higher interaction
    score in the target network.
    """

    # read GRN files
    logger.info("Loading source network.")
    source = read_network(grn_source, full_output)
    if edges and not select_after_join:
        source = source.sort_values(sort_by).tail(edges)
    source.rename(columns=GRN_COLUMNS, inplace=True)

    logger.info("Loading target network.")
    target = read_network(grn_target, full_output)
    if edges and not select_after_join:
        logger.info(f"    Selecting top {edges} edges before calculating difference")
        target = target.sort_values(sort_by).tail(edges)
    target.rename(columns=GRN_COLUMNS, inplace=True)

    # Calculate difference
    logger.info("Calculating differential network.")
    # Fill weight not present in source network with 0
    diff_network = target.join(source, lsuffix="_target", rsuffix="_source").fillna(0)
    diff_network["weight"] = (
        diff_network["weight_target"] - diff_network["weight_source"]
    )

    # Only keep edges that are higher in target network
    diff_network = diff_network[diff_network["weight"] > 0]
    if len(diff_network) == 0:
        logger.error("No differences between networks!")
        sys.exit(1)

    if not full_output:
        diff_network.drop(columns=["weight_target", "weight_source"], inplace=True)

    # Only keep top edges
    if edges and select_after_join:
        logger.info(f"    Selecting top {edges} edges after calculating difference")
        sort_by = GRN_COLUMNS.get(sort_by, sort_by)
        if sort_by != "weight":
            diff_network[sort_by] = (
                diff_network[f"{sort_by}_target"] - diff_network[f"{sort_by}_source"]
            )
        diff_network = diff_network.sort_values(sort_by).tail(edges)

    # split the transcription factor and target gene into 2 columns, make sure they end up
    # in the first columns
    source_target = (
        diff_network.index.to_series()
        .str.split(SEPARATOR, expand=True)
        .rename(columns={0: "source", 1: "target"})
    )
    diff_network = pd.concat((source_target, diff_network), axis=1)
    diff_network.reset_index(drop=True, inplace=True)

    if outfile:
        logger.info("Saving differential network.")
        diff_network.to_csv(outfile, sep="\t", index=False)

    # load into a network with TFs and TGs as nodes, and the interaction scores as edges
    grn = nx.from_pandas_edgelist(diff_network, edge_attr=True, create_using=nx.DiGraph)

    return grn


def target_score(expression_change, targets):
    """
    Calculate the target score, as (mostly) explained in equation 5:
    https://academic.oup.com/nar/article/49/14/7966/6318498#M5
    """
    ts = 0
    for target, weight in targets.items():
        # g: expression score of the target
        g = expression_change[target].score
        # weight: cumulative probability normalized by the length
        score = g * weight
        ts += score
    return ts


def influence_scores(node, grn, expression_change, de_genes, max_steps=2):
    """
    Calculate the influence scores of a transcription factor.

    Parameters
    ----------
    node : str
        Transcription factor name, present in grn as a node
    grn : nx.DiGraph
        A network with gene names as nodes and interaction scores as weights
    expression_change : dict
        A dictionary with interaction scores and log fold changes per transcription factor
    de_genes : list or set or dict
        A list-like with genes present in expression_change that have a score > 0
    max_steps : int
        The maximum number of steps between the TF and the target gene
        (example with 2 steps: TF -> intermediate TF -> target gene)

    Returns
    -------
    tuple
        interaction data of the given transcription factor
    """
    # sum target scores for all genes that are
    # - up to 'max_steps' away from the TF
    # - differentially expressed
    sub_grn = nx.generators.ego_graph(grn, node, radius=max_steps)
    # dijkstra_prob_length cutoff between 0.25 to 0.32 yields the same targets
    paths, weights = dijkstra_prob_length(sub_grn, node, "weight")
    de_targets = {k: v for k, v in weights.items() if k in de_genes}
    targetscore = target_score(expression_change, de_targets)

    pval, target_fc_diff = fold_change_scores(node, grn, expression_change)
    factor_fc = expression_change[node].absfc if node in expression_change else 0
    return (
        node,  # factor
        grn.out_degree(node),  # noqa. direct_targets
        len(paths),  # total_targets
        targetscore,  # target_score
        expression_change[node].score,  # G_score
        factor_fc,  # factor_fc
        pval,  # pval
        target_fc_diff,  # target_fc
    )


def fold_change_scores(node, grn, expression_change):
    """
    Get the Mann-Whitney U p-value of direct targets vs. non-direct targets,
    as well as the difference of the mean fold changes.
    """
    direct_targets = set(grn[node]) & set(expression_change)
    if len(direct_targets) == 0:
        return np.NAN, np.NAN
    non_direct_targets = (set(grn.nodes) & set(expression_change)) - direct_targets
    if len(non_direct_targets) == 0:
        return np.NAN, np.NAN

    target_fc = [expression_change[t].absfc for t in direct_targets]
    non_target_fc = [expression_change[t].absfc for t in non_direct_targets]
    try:
        # auto method prevents recursion errors.
        pval = mannwhitneyu(target_fc, non_target_fc, method="auto")[1]
    except (RecursionError, ValueError) as e:
        pval = np.NAN
        logger.warning(e)
        # logger.warning(
        #     f"Could not calculate p-val (target vs non-target fold-change) for {node}, "
        #     f"targets = {len(target_fc)}, non-target = {len(non_target_fc)}."
        # )
        # logger.warning(f"targets = {target_fc[0:min(5, len(target_fc))]}...")
        # logger.warning(f"non_target = {non_target_fc[0:min(5, len(non_target_fc))]}...")
    target_fc_diff = np.mean(target_fc) - np.mean(non_target_fc)
    return pval, target_fc_diff


def filter_tf(scores_df, network=None, tpmfile=None, tpm=20, overlap=0.98):
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


class Influence(object):
    def __init__(
        self,
        outfile,
        degenes,
        gene_gtf=None,
        grn_source_file=None,
        grn_target_file=None,
        filter_tfs=False,  # variable not exposed in CLI
        edges=100_000,
        ncore=1,
        sort_by="prob",
        padj_cutoff=0.05,
        full_output=False,
        select_after_join=False,
    ):
        self.ncore = ncore
        self.gene_gtf = gene_gtf
        self.full_output = full_output
        self.outfile = outfile
        self.filter_tfs = filter_tfs

        # Load GRNs
        if grn_target_file is None:
            logger.error("You should provide at least an ANANSE target network file!")
            sys.exit(1)

        logger.info(f"Loading network data, using the top {edges} edges")
        if edges and not full_output and sort_by != "prob":
            logger.error(
                f"Sorting by column '{sort_by}' is not possible without the full output!"
            )
            sys.exit(1)
        if grn_source_file is None:
            self.grn = read_network_to_graph(
                grn_target_file, edges=edges, sort_by=sort_by
            )
            logger.warning("You only provided the target network!")
        else:
            outfile = os.path.splitext(self.outfile)[0] + "_diffnetwork.tsv"
            self.grn = difference(
                grn_source_file,
                grn_target_file,
                outfile,
                edges,
                sort_by,
                full_output,
                select_after_join,
            )
            logger.info(f"    Differential network has {len(self.grn.edges)} edges.")

        # Load expression file
        self.expression_change = self.read_expression(degenes, padj_cutoff)

    def read_expression(self, fname, padj_cutoff=0.05):
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

        padj_cutoff: float, optional
            cutoff below which genes are flagged as differential, default is 0.05

        Returns
        -------
        dict
            namedtuples of scores, absolute fold change and "real" (directional) fold change.
        """
        cutoff = 0.6  # fraction of overlap that is "good enough"
        logger.info(
            f"Loading expression data, using genes with an adjusted p-value below {padj_cutoff}"
        )

        df = pd.read_table(
            fname,
            index_col=0,
            header=0,
            dtype=str,
        )
        for col in ["log2FoldChange", "padj"]:
            if col not in df.columns:
                logger.error(
                    f"Column '{col}' not in differential gene expression file!"
                )
                sys.exit(1)
        df = df[["log2FoldChange", "padj"]].astype(float)

        # convert to gene names if overlap is poor
        network_genes = set(self.grn.nodes)
        df_genes = set(df.index)
        pct_overlap = len(network_genes & df_genes) / min(
            len(network_genes), len(df_genes)
        )
        logger.debug(
            f"{int(100 * pct_overlap)}% of genes found in DE genes and network(s)"
        )
        if pct_overlap < cutoff and self.gene_gtf is not None:
            logger.warning(
                "Converting genes in differential expression table to HGNC symbols"
            )
            backup_pct_overlap = pct_overlap
            backup_df = df.copy()

            gp = genomepy.Annotation(self.gene_gtf, quiet=True)
            tid2gid = gp.gtf_dict("transcript_id", "gene_id")
            tid2name = gp.gtf_dict("transcript_id", "gene_name")
            gid2name = gp.gtf_dict("gene_id", "gene_name")
            df = (
                df.rename(index=tid2name)
                .rename(index=tid2gid)
                .rename(index=gid2name)
                .reset_index()
            )
            # take the most significant gene per duplicate (if applicable)
            df = df.groupby("index").min("padj")

            df_genes = set(df.index)
            pct_overlap = len(network_genes & df_genes) / min(
                len(network_genes), len(df_genes)
            )
            logger.debug(
                f"{int(100 * pct_overlap)}% of genes found in DE genes and network(s)"
            )
            if pct_overlap <= backup_pct_overlap:
                df = backup_df

        # unnamed genes cannot be matched
        df.dropna(inplace=True)
        # merge duplicate genes
        dup_df = df[df.index.duplicated()]
        if len(dup_df) > 0:
            logger.warning(
                "Duplicated gene names detected in differential expression file e.g. "
                f"'{str(dup_df.index[0])}'. Averaging values for duplicated genes..."
            )
            df = df.groupby(by=df.index, dropna=True).mean(0)

        overlap = len(network_genes & set(df.index))
        if overlap == 0:
            logger.error(
                "Gene names don't overlap between the "
                "differential gene expression file and network file(s)!"
            )
            if self.gene_gtf is None:
                logger.info(
                    "If you provide a GTF file we can try to convert genes to HGNC symbols"
                )
            sys.exit(1)
        logger.debug(
            f"{overlap} genes overlap between the "
            "differential expression file and the network file(s)"
        )

        # absolute fold change
        df["fc"] = df["log2FoldChange"].abs()

        # get the gscore (absolute fold change if significantly differential)
        df["score"] = df["fc"] * (df["padj"] < padj_cutoff)

        expression_change = dict()
        for k, row in df.iterrows():
            expression_change[k] = Expression(
                score=row.score, absfc=row.fc, realfc=row.log2FoldChange
            )
        return expression_change

    def run_target_score(self):
        """Run target score for all TFs."""

        tfs = [
            node for node in self.grn.nodes() if self.grn.out_degree(node) > 0  # noqa
        ]
        logger.info(f"Differential network contains {len(tfs)} transcription factors.")

        # differentially expressed TFs
        de_tfs = [tf for tf in tfs if tf in self.expression_change]
        if len(de_tfs) == 0:
            logger.error(
                "No overlapping transcription factors found between the network file(s) "
                "(-s/--source, -t/--target) and the differential expression data (-d/--degenes)!"
            )
            sys.exit(1)

        # TODO: should 'realfc' be 'score' (padj<cutoff), or even 'absfc'?
        de_tfs = set(tf for tf in de_tfs if self.expression_change[tf].realfc > 0)
        if len(de_tfs) == 0:
            # expression_change[tf].score > 0 == differentially expressed
            logger.error("No increasingly expressed TFs found!")
            sys.exit(1)
        else:
            logger.info(f"    Out of these, {len(de_tfs)} are upregulated.")

        # differentially expressed genes
        genes = self.grn.nodes
        logger.info(f"Differential network contains {len(genes)} genes.")
        de_genes = set(
            g
            for g in self.expression_change
            if g in genes and self.expression_change[g].score > 0
        )
        logger.info(f"    Out of these, {len(de_genes)} are differentially expressed.")

        tmpdir = mytmpdir()
        tmpfile = os.path.join(tmpdir, os.path.basename(self.outfile))
        influence_file = open(tmpfile, "w")
        influence_file.write(
            "factor\tdirect_targets\ttotal_targets\ttarget_score\tG_score\tfactor_fc\tpval\ttarget_fc\n"
        )

        try:
            if self.ncore > 1:
                pool = mp.Pool(self.ncore)
                jobs = []
                for tf in de_tfs:
                    jobs.append(
                        pool.apply_async(
                            influence_scores,
                            (tf, self.grn, self.expression_change, de_genes),
                        )
                    )
                pool.close()
                with tqdm(total=len(jobs)) as pbar:
                    for j in jobs:
                        print(*j.get(), file=influence_file, sep="\t")
                        pbar.update(1)
                pool.join()

            else:
                for tf in tqdm(de_tfs):
                    line = influence_scores(
                        tf, self.grn, self.expression_change, de_genes
                    )
                    print(*line, file=influence_file, sep="\t")

            influence_file.close()
            shutil.move(tmpfile, self.outfile)

        except Exception as e:
            pool = None  # noqa: force garbage collection on orphaned workers
            if "multiprocessing" in e.__repr__():
                msgs = [
                    str(e),
                    "The error seems to be related to multiprocessing.",
                    "In some cases running `ananse influence` with `-n 1` will solve this issue.",
                    "If it doesn't, please file a bug report (with the output of the command run with `-n 1`) at:",
                    "https://github.com/vanheeringen-lab/ANANSE/issues",
                ]
                _ = [logger.error(msg) for msg in msgs]
                sys.exit(1)
            raise e

    def run_influence_score(self, influence_file, fin_expression=None):
        """Calculate influence score from target score and gscore"""

        scores_df = pd.read_table(influence_file, index_col="factor")
        scores_df["target_score_scaled"] = minmax_scale(
            rankdata(scores_df["target_score"], method="dense")
        )
        scores_df["G_score_scaled"] = minmax_scale(
            rankdata(scores_df["G_score"], method="dense")
        )
        scores_df["influence_score_raw"] = scores_df.target_score + scores_df.G_score
        scores_df["influence_score"] = minmax_scale(
            rankdata(
                scores_df.target_score_scaled + scores_df.G_score_scaled, method="dense"
            )
        )
        scores_df.sort_values("influence_score", inplace=True, ascending=False)
        scores_df = scores_df[
            [
                "influence_score",
                "influence_score_raw",
                "target_score",
                "target_score_scaled",
                "G_score",
                "G_score_scaled",
                "direct_targets",
                "factor_fc",
            ]
        ]

        scores_df.to_csv(self.outfile, sep="\t")

        if self.filter_tfs:
            scores_df2 = filter_tf(
                network=self.grn, scores_df=scores_df, tpmfile=fin_expression
            )
            scores_df2.to_csv(
                ".".join(self.outfile.split(".")[:-1]) + "_filtered.txt", sep="\t"
            )

    def run_influence(self, fin_expression=None):

        logger.info("Calculating target scores.")
        self.run_target_score()

        logger.info("Calculating influence scores.")
        self.run_influence_score(self.outfile, fin_expression)
