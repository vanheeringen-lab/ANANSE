"""Predict TF influence score"""
import os
import shutil
import sys
import warnings
import genomepy
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


def read_network(
    fname, edges=100_000, interactions=None, sort_by="prob", full_output=False
):
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
    fname, edges=100_000, interactions=None, sort_by="prob", full_output=False
):
    """
    Read a network file and return a networkx DiGraph.

    Subsets the graph to interactions if given, else to the number top interactions given by edges.
    If both interactions and edges are none, return the whole graph.

    Peak memory usage is about ~5GB per million edges (tested with 0.1m, 1m and 10m edges)
    """
    rnet = read_network(
        fname,
        edges=edges,
        interactions=interactions,
        sort_by=sort_by,
        full_output=full_output,
    )
    # sort on selection variable
    rnet.sort_values(sort_by, ascending=False, inplace=True)
    if interactions is not None:
        rnet = rnet[rnet.index.isin(interactions)]
    elif edges is not None:
        rnet = rnet.head(edges)

    # split the transcription factor and target gene into 2 columns
    rnet[["source", "target"]] = rnet.index.to_series().str.split(
        SEPARATOR, expand=True
    )
    rnet.reset_index(drop=True, inplace=True)

    # rename the columns
    rnet.rename(
        columns={
            "prob": "weight",
            # "weighted_binding": "weighted_binding",
            # "tf_expression": "tf_expression",
            "target_expression": "tg_expression",
            "activity": "tf_activity",
        },
        inplace=True,
    )

    # load into a network with TFs and TGs as nodes, and the interaction scores as edges
    grn = nx.from_pandas_edgelist(rnet, edge_attr=True, create_using=nx.DiGraph)

    return grn

def difference(
    grn_source, grn_target, outfile, edges=100_000, column="prob", full_output=False, select_after_join=False
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
    source = read_network(grn_source, edges=None, full_output=full_output)
    logger.info("Loading target network.")
    target = read_network(grn_target, edges=None, full_output=full_output)

    if edges and not select_after_join:
        logger.info(f"Selecting top {edges} edges before calculating difference")
        source = source.sort_values(column).tail(edges)
        target = target.sort_values(column).tail(edges)


    # Calculate difference
    logger.info("Calculating differential network.")
    # Fill weight not present in source network with 0
    diff_network = target.join(source, lsuffix="_target", rsuffix="_source").fillna(0)
    diff_network["weight"] = (
        diff_network[f"{column}_target"] - diff_network[f"{column}_source"]
    )

    # Only keep edges that are higher in target network
    diff_network = diff_network[diff_network["weight"] > 0]

    if not full_output:
        diff_network.drop(
            columns=[f"{column}_target", f"{column}_source"], inplace=True
        )

    # Only keep top edges
    if edges and select_after_join:
        logger.info(f"Selecting top {edges} edges after calculating difference")
        diff_network = diff_network.sort_values("weight").tail(edges)

    logger.info("Saving differential network.")
    diff_network.to_csv(outfile, sep="\t")

    # split the transcription factor and target gene into 2 columns
    diff_network[["source", "target"]] = diff_network.index.to_series().str.split(
        SEPARATOR, expand=True
    )
    diff_network.reset_index(drop=True, inplace=True)
    # load into a network with TFs and TGs as nodes, and the interaction scores as edges
    grn = nx.from_pandas_edgelist(diff_network, edge_attr=True, create_using=nx.DiGraph)

    return grn


def get_weight(grn, path):
    weight = np.cumprod([grn[s][t]["weight"] for s, t in zip(path, path[1:])])[-1]
    # Add weight correction for the length of the path
    weight = weight / (len(path) - 1)
    return weight


def target_score(node, grn, expression_change, targets):
    """
    Calculate the target score, as (mostly) explained in equation 5:
    https://academic.oup.com/nar/article/49/14/7966/6318498#M5
    """
    ts = 0
    for target in targets:
        # Calculate all paths from TF to target to select to path with the highest multiplied weight
        all_paths = {}
        for path in nx.all_simple_paths(grn, node, target, cutoff=2):
            all_paths[tuple(path)] = get_weight(grn, path)
        if len(all_paths) > 0:
            path, weight = sorted(all_paths.items(), key=lambda pw: pw[1])[-1]
            # the level (or the number of steps) that gene is away from transcription factor
            # TODO: this is the number of nodes. 1 higher than the number of steps!
            pathlen = len(path)
            # expression score of the target
            g = expression_change[target].score if target in expression_change else 0
            # TODO: why divide by path length AGAIN? (already happens in the weight function)
            # TODO: maybe this was left in by accident?
            # TODO: see https://github.com/vanheeringen-lab/ANANSE/commit/ba67ebb8e7bafd1df13fb439485b6a590482e924
            score = g / pathlen * weight
            ts += score
    return ts


# TODO: this version is ~O(n^2) faster, and returns near identical scores
# def target_score(grn, expression_change, targets):
#     ts = 0
#     for target in targets:
#         path = targets[target]
#         weight = get_weight(grn, path)
#         score = expression_change[target].score / len(path) * weight
#         ts += score
#     return ts


def influence_scores(node, grn, expression_change, de_genes):
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

    Returns
    -------
    tuple
        interaction data of the given transcription factor
    """
    # get all genes that are
    # - direct or indirectly targeted by the TF (cutoff + weight)
    # - not the TF itself (because we divide by len(path)-1 elsewhere)
    # - differentially expressed
    targets = nx.single_source_dijkstra(grn, node, cutoff=2, weight=None)[1]  # noqa
    _ = targets.pop(node, None)

    de_targets = {k: v for k, v in targets.items() if k in de_genes}
    targetscore = target_score(node, grn, expression_change, de_targets)

    pval, target_fc_diff = fold_change_scores(node, grn, expression_change)
    factor_fc = expression_change[node].absfc if node in expression_change else 0
    return (
        node,  # factor
        grn.out_degree(node),  # noqa. directTargets
        len(targets),  # totalTargets
        targetscore,  # targetscore
        expression_change[node].score,  # Gscore
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
        # asymptotic method prevents recursion errors.
        # TODO: review when scipy closes https://github.com/scipy/scipy/issues/14622
        pval = mannwhitneyu(target_fc, non_target_fc, method="asymptotic")[1]
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
        filter_tfs=False,  # TODO: variable not exposed in CLI
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
                "out_diffnetwork.txt",
                edges=edges,
                column=sort_by,
                full_output=full_output,
                select_after_join=select_after_join,
            )

            if len(self.grn.edges) == 0:
                logger.error("No differences between networks!")
                sys.exit(1)
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
        df = df[["log2FoldChange", "padj"]].dropna()  # removes unneeded data
        df = df.astype(float)

        network_genes = set(self.grn.nodes)
        pct_overlap = len(network_genes & set(df.index)) / len(network_genes)
        logger.debug(
            f"{int(100 * pct_overlap)}% of genes found in DE genes and network(s)"
        )
        if pct_overlap < cutoff and self.gene_gtf is not None:
            logger.warning(
                "Converting genes in differential expression table to HGNC symbols"
            )
            backup_pct_overlap = pct_overlap
            backup_df = df.copy()

            gp = genomepy.Annotation(self.gene_gtf)
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

            pct_overlap = len(network_genes & set(df.index)) / len(network_genes)
            logger.debug(
                f"{int(100 * pct_overlap)}% of genes found in DE genes and network(s)"
            )
            if pct_overlap <= backup_pct_overlap:
                df = backup_df

        overlap = len(network_genes & set(df.index))
        if overlap == 0:
            logger.error(
                "Gene names don't overlap between the "
                "differential gene expression file and network file(s)!"
            )
            if self.gene_gtf is None:
                logger.info(
                    "If you provide a GTF file we can try to convert "
                    "the gene expression file to HGNC names."
                )
            sys.exit(1)
        logger.debug(
            f"{overlap} genes overlap between the "
            "differential expression file and the network file(s)"
        )

        # check for duplicated index rows and return an error (but continue running)
        dup_df = df[df.index.duplicated()]
        if len(dup_df) > 0:
            dupped_gene = str(dup_df.index[0])
            logger.warning(
                f"Duplicated gene names detected in differential expression file e.g. '{dupped_gene}'. "
                "Averaging values for duplicated genes..."
            )
            df = df.groupby(by=df.index, dropna=True).mean(0)

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

        # TODO: should 'realfc' this not be 'score' (padj<cutoff), or even 'absfc'?
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
            "factor\tdirect_targets\ttotalTargets\ttarget_score\tG_score\tfactor_fc\tpval\ttarget_fc\n"
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

        scores_df = pd.read_table(influence_file, index_col=0)
        scores_df["target_score_scaled"] = minmax_scale(
            rankdata(scores_df["target_score"], method="dense")
        )
        scores_df["G_score"] = scores_df["G_score"]
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
