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


def read_top_interactions(fname, edges=100_000, sort_by="prob"):
    """
    Read network file and return the top interactions.

    Incorporates the combined metric (probability score) by default to select the most
    certain interactions, however other sorting options such as the weighted binding score
    can be applied too.
    """
    # read the GRN file
    rnet = pd.read_csv(fname, usecols=["tf_target", sort_by], sep="\t")
    # sort based on the selected value (default probability score)
    rnet.sort_values(sort_by, ascending=False, inplace=True)
    # takes all edges if there are less than requested
    rnet = rnet.head(edges)
    top_int = set(rnet["tf_target"])
    return top_int


def read_network(
    fname, edges=100_000, interactions=None, sort_by="prob", full_output=False
):
    """Read network file and return networkx DiGraph"""

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
    rnet = pd.read_csv(
        fname,
        sep="\t",
        usecols=data_columns,
        dtype="float64",
        converters={"tf_target": str},
    )  # read the GRN file
    # sort on selection variable
    rnet.sort_values(sort_by, ascending=False, inplace=True)
    if interactions is None:
        rnet = rnet.head(edges)  # no interactions so take the top head of interactions
    else:
        rnet = rnet[rnet.tf_target.isin(interactions)]

    grn = nx.DiGraph()  # initiate empty network
    for _, row in rnet.iterrows():
        source, target = row["tf_target"].split(SEPARATOR, 1)
        try:
            if full_output:
                grn.add_edge(
                    source,
                    target,
                    weight=row["prob"],
                    weighted_binding=row["weighted_binding"],
                    tf_expression=row["tf_expression"],
                    tg_expression=row["target_expression"],
                    tf_activity=row["activity"],
                )
            else:
                grn.add_edge(source, target, weight=row["prob"])  # , n=1)
        except Exception:
            logger.error(f"Could not parse edge weight {source}:{target}")
            raise
    return grn


def difference(grn_source, grn_target, full_output=False):
    """
    Calculate the network different between two GRNs.

    First take the nodes from both networks, and add
    edges from the target network that are missing in the source network.
    Then add edges present in both but with a higher interaction
    score in the target network.
    """
    grn_diff = nx.create_empty_copy(grn_source)  # copying source GRN nodes
    nodes_target = nx.create_empty_copy(grn_target)  # copying target GRN nodes
    # merge the nodes if there are differences
    # (which can happen) when taking the header instead of the --union-interactions flag
    grn_diff.add_nodes_from(list(nodes_target.nodes))

    # lets check if the full GRN output is loaded and if so output all attributes to the diffnetwork:
    if full_output:
        logger.info("Calculating differential GRN with full output")
        # lets load all the  edges into the diffnetwork
        for u, v, ddict in grn_target.edges(data=True):
            # u = source node, v = target node, ddict = dictionary of edge attributes
            # calculate the weight difference and output all attributes
            source_weight = grn_source.edges[u, v]["weight"]
            target_weight = grn_target.edges[u, v]["weight"]
            diff_weight = target_weight - source_weight
            if diff_weight > 0:
                # if the interaction probability is higher in the target than in the
                # source, add the interaction to the diffnetwork:
                grn_diff.add_edge(
                    u,
                    v,
                    weight=diff_weight,
                    source_weight=source_weight,
                    target_weight=target_weight,
                    # neglogweight=-np.log(diff_weight),
                    # n=1,
                    tf_expr_diff=(
                        (ddict["tf_expression"]) - (grn_source[u][v]["tf_expression"])
                    ),
                    tf_expr_target=ddict["tf_expression"],
                    tf_expr_source=grn_source[u][v]["tf_expression"],
                    tg_expr_diff=(
                        (ddict["tg_expression"]) - (grn_source[u][v]["tg_expression"])
                    ),
                    tg_expr_target=ddict["tg_expression"],
                    tg_expr_source=grn_source[u][v]["tg_expression"],
                    wb_diff=(
                        (ddict["weighted_binding"])
                        - (grn_source[u][v]["weighted_binding"])
                    ),
                    target_wb=ddict["weighted_binding"],
                    source_wb=grn_source[u][v]["weighted_binding"],
                    TF_act_diff=(
                        (ddict["tf_activity"]) - (grn_source[u][v]["tf_activity"])
                    ),
                    TF_act_target=ddict["tf_activity"],
                    TF_act_source=grn_source[u][v]["tf_activity"],
                )
    else:
        logger.info("Calculating differential GRN")
        # if only the weight is loaded, lets load only that in the diffnetwork:
        for (u, v) in grn_target.edges():
            source_weight = grn_source.edges[u, v]["weight"]
            target_weight = grn_target.edges[u, v]["weight"]
            diff_weight = target_weight - source_weight
            if diff_weight > 0:
                grn_diff.add_edge(
                    u,
                    v,
                    weight=diff_weight,
                    # n=1,  # neglogweight=-np.log(diff_weight)
                )
    return grn_diff


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
        pval = mannwhitneyu(target_fc, non_target_fc)[1]
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
    ):
        self.ncore = ncore
        self.gene_gtf = gene_gtf
        self.full_output = full_output
        self.outfile = outfile
        self.filter_tfs = filter_tfs

        # Load GRNs
        if grn_target_file is None and grn_source_file is None:
            logger.error("You should provide at least one ANANSE network file!")
            sys.exit(1)
        logger.info(f"Loading network data, using the top {edges} edges")
        if grn_source_file is None:
            self.grn = read_network(grn_target_file, edges=edges)
            logger.warning("You only provided the target network!")
        elif grn_target_file is None:
            self.grn = read_network(grn_source_file, edges=edges)
            logger.warning("You only provided the source network!")
        else:
            top_int_source = read_top_interactions(
                grn_source_file, edges=edges, sort_by=sort_by
            )
            top_int_target = read_top_interactions(
                grn_target_file, edges=edges, sort_by=sort_by
            )
            # combine the top edges of each input network
            top_int = set.union(top_int_source, top_int_target)
            grn_source = read_network(
                grn_source_file,
                interactions=top_int,
                sort_by=sort_by,
                full_output=full_output,
            )
            grn_target = read_network(
                grn_target_file,
                interactions=top_int,
                sort_by=sort_by,
                full_output=full_output,
            )
            self.grn = difference(grn_source, grn_target, full_output=full_output)
            if len(self.grn.edges) == 0:
                raise ValueError("No differences between networks!")
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
                raise ValueError(
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

    def save_reg_network(self, filename, full_output=False):
        """Save the network difference between two cell types to a file."""
        # check if all keys are found
        n = list(self.grn.edges)[0]
        keys = self.grn.edges[n].keys()
        with open(filename, "w") as nw:
            if full_output and "source_wb" in keys:
                logger.info("output full diff network")
                header = [
                    "tf",
                    "target",
                    "weight",
                    "weight_source",
                    "weight_target",
                    "tf_expr_diff",
                    "tf_expr_source",
                    "tf_expr_target",
                    "tg_expr_diff",
                    "tg_expr_source",
                    "tg_expr_target",
                    "wb_diff",
                    "wb_source",
                    "wb_target",
                    "tf_act_diff",
                    "source_tf_act",
                    "target_tf_act",
                ]
                nw.write("\t".join(header) + "\n")
                for (u, v, ddict) in self.grn.edges(data=True):
                    cols = [
                        u,
                        v,
                        ddict["weight"],
                        ddict["source_weight"],
                        ddict["target_weight"],
                        ddict["tf_expr_diff"],
                        ddict["tf_expr_source"],
                        ddict["tf_expr_target"],
                        ddict["tg_expr_diff"],
                        ddict["tg_expr_source"],
                        ddict["tg_expr_target"],
                        ddict["wb_diff"],
                        ddict["source_wb"],
                        ddict["target_wb"],
                        ddict["TF_act_diff"],
                        ddict["TF_act_source"],
                        ddict["TF_act_target"],
                    ]
                    nw.write("\t".join(str(v) for v in cols) + "\n")
            else:
                nw.write("tf\ttarget\tweight\n")
                for (u, v, ddict) in self.grn.edges(data=True):
                    nw.write(u + "\t" + v + "\t" + str(ddict["weight"]) + "\n")

    def run_target_score(self):
        """Run target score for all TFs."""

        tfs = [node for node in self.grn.nodes() if self.grn.out_degree(node) > 0]
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
            logger.info(f"    Out of these, {len(de_tfs)} are increasingly expressed.")

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
            "factor\tdirectTargets\ttotalTargets\ttargetscore\tGscore\tfactor_fc\tpval\ttarget_fc\n"
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
        scores_df["targetScaled"] = minmax_scale(
            rankdata(scores_df["targetscore"], method="dense")
        )
        scores_df["Gscore"] = scores_df["Gscore"]
        scores_df["GscoreScaled"] = minmax_scale(
            rankdata(scores_df["Gscore"], method="dense")
        )
        scores_df["sum"] = scores_df.targetscore + scores_df.Gscore
        scores_df["sumScaled"] = minmax_scale(
            rankdata(scores_df.targetScaled + scores_df.GscoreScaled, method="dense")
        )
        scores_df.sort_values("sumScaled", inplace=True, ascending=False)
        scores_df = scores_df[
            [
                "targetscore",
                "targetScaled",
                "Gscore",
                "GscoreScaled",
                "sum",
                "sumScaled",
                "directTargets",
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
        logger.info("Saving differential network.")
        diffnetwork = os.path.splitext(self.outfile)[0] + "_diffnetwork.tsv"
        self.save_reg_network(diffnetwork, self.full_output)

        logger.info("Calculating target scores.")
        self.run_target_score()

        logger.info("Calculating influence scores.")
        self.run_influence_score(self.outfile, fin_expression)
