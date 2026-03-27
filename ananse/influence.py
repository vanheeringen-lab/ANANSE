"""Predict TF influence score"""
import multiprocessing as mp
import os
import shutil
import sys
import warnings
from collections import namedtuple
from typing import Union

import genomepy
import networkx as nx
import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import mannwhitneyu, rankdata
from sklearn.preprocessing import minmax_scale
from tqdm.auto import tqdm

from ananse import SEPARATOR
from ananse.nx import dijkstra_prob_length
from ananse.utils import load_whitelist, mytmpdir

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


def _read_network(fname, full_output=False):
    """
    Read a network file and return a DataFrame.

    Peak memory usage is about ~5GB per million edges
    (tested with 0.1m, 1m and 10m edges).
    """
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
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=pd.errors.ParserWarning)
        rnet = pd.read_csv(
            fname,
            sep="\t",
            usecols=data_columns,
            dtype="float64",
            converters={"tf_target": str},
            index_col="tf_target",
        )
    rnet.rename(columns=GRN_COLUMNS, inplace=True)
    return rnet


def read_network_to_graph(
    fname,
    sort_by="prob",
    edges: Union[int, None] = 100_000,
    whitelist=None,
    interactions=None,
    full_output=False,
):
    """
    Read a network file and return a networkx DiGraph.

    Subsets the graph to interactions if given, else to the number top interactions given by edges.
    Optionally accepts a list of genes or TF—target interactions to include with the top edges.

    If both interactions and edges are none, return the whole graph.
    """
    rnet = _read_network(fname, full_output)
    if interactions is not None:
        rnet = rnet[rnet.index.isin(interactions)]
    if edges is not None:
        rnet = _filter_network_edges(rnet, sort_by, edges, whitelist)

    rnet = _separate_index(rnet)

    # load into a network with TFs and TGs as nodes, and the interaction scores as edges
    grn = nx.from_pandas_edgelist(rnet, edge_attr=True, create_using=nx.DiGraph)

    return grn


def difference(
    grn_source,
    grn_target,
    sort_by="prob",
    edges=100_000,
    select_after_join=False,
    whitelist=None,
    full_output=False,
    outfile=None,
):
    """
    Calculate the network differences between two GRNs.

    First take the nodes from both networks, and add
    edges from the target network that are missing in the source network.
    Then add edges present in both but with a higher interaction
    score in the target network.
    """
    # read GRN files
    logger.info("Loading source network.")
    source = _read_network(grn_source, full_output)
    if edges and not select_after_join:
        source = _filter_network_edges(source, sort_by, edges, whitelist)

    logger.info("Loading target network.")
    target = _read_network(grn_target, full_output)
    if edges and not select_after_join:
        target = _filter_network_edges(target, sort_by, edges, whitelist)
        n = max(len(source), len(target))
        logger.info(f"    Selected top {n} edges before calculating difference")

    # Calculate difference
    logger.info("Calculating differential network.")
    # Fill weight not present in source network with 0
    diff_network = target.join(source, lsuffix="_target", rsuffix="_source").fillna(0)
    diff_network["weight"] = (
        diff_network["weight_target"] - diff_network["weight_source"]
    )

    # Only keep edges that are higher in target network
    diff_network = diff_network[diff_network["weight"] > 0]
    # if whitelist is None:
    #     diff_network = diff_network[diff_network["weight"] > 0]
    # else:
    #     tfs = set(diff_network.index.str.startswith(whitelist).index)
    #     targets = set(diff_network.index.str.endswith(whitelist).index)
    #     weights = set(diff_network[diff_network["weight"] > 0].index)
    #     diff_network = diff_network[diff_network.index.isin(tfs | targets | weights)]
    #
    if len(diff_network) == 0:
        logger.error("No differences between networks!")
        sys.exit(1)

    if not full_output:
        diff_network.drop(columns=["weight_target", "weight_source"], inplace=True)

    # Only keep top edges
    if edges and select_after_join:
        if sort_by not in ["weight", "prob"]:
            sort_by = GRN_COLUMNS.get(sort_by, sort_by)
            diff_network[sort_by] = (
                diff_network[f"{sort_by}_target"] - diff_network[f"{sort_by}_source"]
            )
        diff_network = _filter_network_edges(diff_network, sort_by, edges, whitelist)
        logger.info(
            f"    Selected top {len(diff_network)} edges after calculating difference"
        )

    diff_network = _separate_index(diff_network)

    if outfile:
        logger.info("Saving differential network.")
        diff_network.to_csv(outfile, sep="\t", index=False)

    # load into a network with TFs and TGs as nodes, and the interaction scores as edges
    grn = nx.from_pandas_edgelist(diff_network, edge_attr=True, create_using=nx.DiGraph)

    return grn


def _filter_network_edges(df, sort_by: str, n_edges: int, whitelist: tuple = None):
    """
    sort a dataframe by a column and filter to a number of edges.
    Optionally accepts a list of TFs, targets or TF—target interactions
    (present in the index) to include in the final network.
    """
    sort_by = GRN_COLUMNS.get(sort_by, sort_by)
    if whitelist is None:
        df = df.sort_values(sort_by).tail(n_edges)
    else:
        df.sort_values(sort_by, inplace=True)
        tfs = set(df[df.index.str.startswith(whitelist)].index)
        targets = set(df[df.index.str.endswith(whitelist)].index)
        tail = set(df.tail(n_edges).index)
        df = df[df.index.isin(tfs | targets | tail)]
    return df


def _separate_index(df):
    """
    split the df index into 2 columns (source and target),
    returns a df starting with these columns
    """
    source_target = (
        df.index.to_series()
        .str.split(SEPARATOR, expand=True)
        .rename(columns={0: "source", 1: "target"})
    )
    df = pd.concat((source_target, df), axis=1)
    df.reset_index(drop=True, inplace=True)
    return df


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
    target_fc_diff = np.mean(target_fc) - np.mean(non_target_fc)
    return pval, target_fc_diff


# def filter_tf(scores_df, network=None, tpmfile=None, tpm=20, overlap=0.98):
#     """Filter TFs:
#     1) it have high expression in origin cell type;
#     2) 98% of its target genes are also regulated by previous TFs.
#     """
#
#     tpmscore = {}
#     with open(tpmfile) as tpf:
#         next(tpf)
#         for line in tpf:
#             tpmscore[line.split()[0]] = float(line.split()[1])
#
#     tftarget = {}
#     for tf in scores_df.index:
#         tftarget[tf] = set(network[tf]) if tf in network else set()
#
#     ltf = list(scores_df.index)
#
#     keeptf = []
#     for i in ltf:
#         passtf = []
#         if len(tftarget[i]) > 0:
#             for j in ltf[: ltf.index(i)]:
#                 if len(tftarget[i] & tftarget[j]) / len(tftarget[i]) > overlap:
#                     break
#                 else:
#                     passtf.append(j)
#             if passtf == ltf[: ltf.index(i)] and i in tpmscore and tpmscore[i] < tpm:
#                 keeptf.append(i)
#     scores_df = scores_df.loc[keeptf]
#     scores_df.sort_values("sumScaled", inplace=True, ascending=False)
#     return scores_df


class Influence(object):
    grn = None
    expression_change = None

    def __init__(
        self,
        outfile,
        degenes,
        gene_gtf=None,
        grn_source_file=None,
        grn_target_file=None,
        sort_by="prob",
        edges=100_000,
        whitelist=None,
        select_after_join=False,
        padj_cutoff=0.05,
        ncore=1,
        full_output=False,
        # filter_tfs=False,  # variable not exposed in CLI
    ):
        self.ncore = ncore
        self.full_output = full_output
        self.outfile = outfile
        # self.filter_tfs = filter_tfs

        if grn_target_file is None:
            logger.error("You should provide at least an ANANSE target network file!")
            sys.exit(1)

        if edges and not full_output and sort_by != "prob":
            logger.error(
                f"Sorting by column '{sort_by}' is not possible without the full output!"
            )
            sys.exit(1)

        # Load self.grn
        self.read_networks(
            grn_source_file,
            grn_target_file,
            sort_by,
            edges,
            select_after_join,
            whitelist,
        )

        # Load self.expression_change
        self.read_expression(degenes, padj_cutoff, gene_gtf)

    def read_networks(
        self,
        grn_source_file,
        grn_target_file,
        sort_by,
        edges,
        select_after_join,
        whitelist,
    ):
        if whitelist is not None:
            whitelist = load_whitelist(whitelist)

        logger.info(f"Loading network data, using the top {edges} edges")
        if grn_source_file is None:
            self.grn = read_network_to_graph(
                grn_target_file,
                sort_by,
                edges,
                whitelist,
            )
            logger.warning("You only provided the target network!")
        else:
            outfile = os.path.splitext(self.outfile)[0] + "_diffnetwork.tsv"
            self.grn = difference(
                grn_source_file,
                grn_target_file,
                sort_by,
                edges,
                select_after_join,
                whitelist,
                self.full_output,
                outfile,
            )
            logger.info(f"    Differential network has {len(self.grn.edges)} edges.")

    def read_expression(self, fname, padj_cutoff=0.05, gene_gtf=None):
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

        gene_gtf: str, optional
            GTF file used to convert gene IDs to gene names.
            Only used if the overlap between DE genes and the network genes is low.

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
        if pct_overlap < cutoff and gene_gtf is not None:
            logger.warning(
                "Converting genes in differential expression table to HGNC symbols"
            )
            backup_pct_overlap = pct_overlap
            backup_df = df.copy()

            gp = genomepy.Annotation(gene_gtf, quiet=True)
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
            if gene_gtf is None:
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
        self.expression_change = expression_change

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

    def run_influence_score(self):  # , influence_file, fin_expression=None):
        """Calculate influence score from target score and gscore"""

        df = pd.read_table(self.outfile, index_col="factor")  # influence_file
        df["target_score_scaled"] = minmax_scale(
            rankdata(df["target_score"], method="dense")
        )
        df["G_score_scaled"] = minmax_scale(rankdata(df["G_score"], method="dense"))
        df["influence_score_raw"] = df.target_score + df.G_score
        df["influence_score"] = minmax_scale(
            rankdata(df.target_score_scaled + df.G_score_scaled, method="dense")
        )

        df.sort_values("influence_score", inplace=True, ascending=False)
        df = df[
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
        df.to_csv(self.outfile, sep="\t")

        # if self.filter_tfs:
        #     df2 = filter_tf(
        #         network=self.grn, scores_df=df, tpmfile=fin_expression
        #     )
        #     df2.to_csv(
        #         ".".join(self.outfile.split(".")[:-1]) + "_filtered.txt", sep="\t"
        #     )

    def run_influence(self):  # , fin_expression=None):
        logger.info("Calculating target scores.")
        self.run_target_score()

        logger.info("Calculating influence scores.")
        self.run_influence_score()  # self.outfile, fin_expression)
