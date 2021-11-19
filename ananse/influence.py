"""Predict TF influence score"""
import sys
import warnings
import genomepy
from collections import namedtuple
from loguru import logger
from tqdm import tqdm
import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from sklearn.preprocessing import minmax_scale
from scipy.stats import rankdata, mannwhitneyu

from . import SEPARATOR


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
    rnet.sort_values(GRNsort_column, ascending=False, inplace=True)
    if interactions is None:
        rnet = rnet.head(edges)  # no interactions so take the top head of interactions
    else:
        rnet = rnet[rnet.tf_target.isin(interactions)]

    G = nx.DiGraph()  # initiate empty network
    for _, row in rnet.iterrows():
        source, target = row[0].split(SEPARATOR, 1)
        try:
            if full_output:
                G.add_edge(
                    source,
                    target,
                    weight=row["prob"],
                    weighted_binding=row["weighted_binding"],
                    tf_expression=row["tf_expression"],
                    tg_expression=row["target_expression"],
                    tf_activity=row["activity"],
                )
            else:
                G.add_edge(source, target, weight=row["prob"], n=1)
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
    # merge the nodes if there are differences
    # (which can happen) when taking the header instead of the --union-interactions flag
    DIF.add_nodes_from(list(nodes_target.nodes))

    # lets check if the full GRN output is loaded and if so output all attributes to the diffnetwork:
    if full_output:
        logger.info("Calculating differential GRN with full output")
        # lets load all the  edges into the diffnetwork
        for u, v, ddict in GRN_target.edges(data=True):
            # u = source node, v = target node, ddict = dictionary of edge attributes
            # calculate the weight difference and output all attributes
            source_weight = GRN_source.edges[u, v]["weight"]
            target_weight = GRN_target.edges[u, v]["weight"]
            diff_weight = target_weight - source_weight
            if diff_weight > 0:
                # if the interaction probability is higher in the target than in the
                # source, add the interaction to the diffnetwork:
                DIF.add_edge(
                    u,
                    v,
                    weight=diff_weight,
                    source_weight=source_weight,
                    target_weight=target_weight,
                    neglogweight=-np.log(diff_weight),
                    n=1,
                    tf_expr_diff=(
                        (ddict["tf_expression"]) - (GRN_source[u][v]["tf_expression"])
                    ),
                    tf_expr_target=ddict["tf_expression"],
                    tf_expr_source=GRN_source[u][v]["tf_expression"],
                    tg_expr_diff=(
                        (ddict["tg_expression"]) - (GRN_source[u][v]["tg_expression"])
                    ),
                    tg_expr_target=ddict["tg_expression"],
                    tg_expr_source=GRN_source[u][v]["tg_expression"],
                    wb_diff=(
                        (ddict["weighted_binding"])
                        - (GRN_source[u][v]["weighted_binding"])
                    ),
                    target_wb=ddict["weighted_binding"],
                    source_wb=GRN_source[u][v]["weighted_binding"],
                    TF_act_diff=(
                        (ddict["tf_activity"]) - (GRN_source[u][v]["tf_activity"])
                    ),
                    TF_act_target=ddict["tf_activity"],
                    TF_act_source=GRN_source[u][v]["tf_activity"],
                )
    else:
        logger.info("Calculating differential GRN")
        # if only the weight is loaded, lets load only that in the diffnetwork:
        for (u, v) in GRN_target.edges():
            source_weight = GRN_source.edges[u, v]["weight"]
            target_weight = GRN_target.edges[u, v]["weight"]
            diff_weight = target_weight - source_weight
            if diff_weight > 0:
                DIF.add_edge(
                    u, v, weight=diff_weight, n=1, neglogweight=-np.log(diff_weight)
                )
    return DIF


def targetScore(node, G, expression_change, max_degree=3):
    """Calculate the influence score."""
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
            # the level (or the number of steps) that gene is away from transcription factor
            pathlen = len(path)
            # expression score of the target
            g = expression_change[target].score if target in expression_change else 0
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
    except (RecursionError, ValueError):
        pval = np.NAN
        logger.warning(
            f"Could not calculate p-val (target vs non-target fold-change) for {node}, "
            f"targets = {len(target_fc)}, non-target = {len(non_target_fc)}."
        )
        logger.warning(f"targets = {target_fc[0:min(5, len(target_fc))]}...")
        logger.warning(f"non_target = {non_target_fc[0:min(5, len(non_target_fc))]}...")
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


class Influence(object):
    def __init__(
        self,
        outfile,
        degenes,
        gene_gtf=None,
        GRN_source_file=None,
        GRN_target_file=None,
        filter=False,  # TODO: variable not exposed in CLI
        edges=100_000,
        ncore=1,
        GRNsort_column="prob",
        padj_cutoff=0.05,
        full_output=False,
    ):
        self.ncore = ncore
        self.gene_gtf = gene_gtf
        self.full_output = full_output
        self.GRNsort_column = GRNsort_column

        # Load GRNs
        if GRN_target_file is None and GRN_source_file is None:
            logger.error("You should provide at least one ANANSE network file!")
            sys.exit(1)
        logger.info(f"Loading network data, using the top {edges} edges")
        if GRN_source_file is None:
            self.G = read_network(GRN_target_file, edges=edges)
            logger.warning("You only provided the target network!")
        elif GRN_target_file is None:
            self.G = read_network(GRN_source_file, edges=edges)
            logger.warning("You only provided the source network!")
        else:
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
            if len(self.G.edges) == 0:
                raise ValueError("No differences between networks!")
            logger.info(f"    Differential network has {len(self.G.edges)} edges.")

        # Load expression file
        self.expression_change = self.read_expression(degenes, padj_cutoff)
        self.outfile = outfile
        # Filter TFs
        self.filter = filter

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

        network_genes = set(self.G.nodes)
        overlap = len(network_genes & set(df.index))
        if overlap == 0 and self.gene_gtf is not None:
            logger.warning(
                "Converting genes in differential expression table to HGNC symbols"
            )

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
        n = list(self.G.edges)[0]
        keys = self.G.edges[n].keys()
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
                for (u, v, ddict) in self.G.edges(data=True):
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
                for (u, v, ddict) in self.G.edges(data=True):
                    nw.write(u + "\t" + v + "\t" + str(ddict["weight"]) + "\n")

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
            "factor\tdirectTargets\ttotalTargets\ttargetscore\tGscore\tfactor_fc\tpval\ttarget_fc\n"
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
            rankdata(scores_df["targetscore"], method="dense")
        )
        scores_df.sort_values("targetScaled", inplace=True, ascending=False)

        return self.outfile

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

        if self.filter:
            scores_df2 = filter_TF(
                network=self.G, scores_df=scores_df, tpmfile=fin_expression
            )
            scores_df2.to_csv(
                ".".join(self.outfile.split(".")[:-1]) + "_filtered.txt", sep="\t"
            )

    def run_influence(self, fin_expression=None):
        logger.info("Saving differential network.")
        self.save_reg_network(
            (".".join(self.outfile.split(".")[:-1]) + "_diffnetwork.tsv"),
            full_output=self.full_output,
        )

        logger.info("Calculating target scores.")
        influence_file = self.run_target_score()

        logger.info("Calculating influence scores.")
        self.run_influence_score(influence_file, fin_expression=fin_expression)
