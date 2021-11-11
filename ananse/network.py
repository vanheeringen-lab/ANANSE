"""Build gene regulatory network"""
import os
import math
import re
import shutil
import sys
import warnings
from typing import Union

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
import genomepy
import dask.dataframe as dd
from tempfile import NamedTemporaryFile, mkdtemp
from dask.distributed import progress
from loguru import logger
from pandas import HDFStore
from tqdm.auto import tqdm
import pyranges as pr

from . import SEPARATOR


warnings.filterwarnings("ignore")
PACKAGE_DIR = os.path.dirname(__file__)


class Network(object):
    _tmp_files = []

    def __init__(
        self,
        ncore=1,
        genome="hg38",
        gene_bed=None,
        include_promoter=False,
        include_enhancer=True,
        full_output=False,
    ):
        """
        infer cell type-specific gene regulatory network

        Parameters
        ----------
            ncore : int
                Specifies the number of threads to use during analysis. (default: 1)
            genome : str
                The genome that is used for the gene annotation and the enhancer location. (default: "hg38")
            gene_bed : str, optional
                Gene annotation for the genome specified with -g as a 12 column BED file. (default: None)
            include_promoter : bool
                Include or exclude promoter peaks (<= TSS +/- 2kb) in network inference. (default: False)
            include_enhancer : bool
                Include or exclude enhancer peaks (> TSS +/- 2kb) in network inference. (default: True)
            full_output : bool
                export all variables to the GRN file or by default only the TF_target + prob score
        """
        self.ncore = ncore
        self.genome = genome
        self.gene_bed = get_bed(gene_bed, genome)
        self.include_promoter = include_promoter
        self.include_enhancer = include_enhancer
        self.full_output = full_output

    @staticmethod
    def unique_enhancers(fname):
        """Extract a list of unique enhancers.

        Parameters
        ----------
        fname : str
            File name of a tab-separated file that contains an 'enhancer' column.

        Returns
        -------
            PyRanges object with enhancers
        """
        logger.info("reading enhancers")

        # Read enhancers from binding file
        header = pd.read_table(fname, nrows=0)
        idx = header.columns.get_loc("enhancer")
        skiprows = 1
        chunksize = 2_000_000
        enhancers = np.array([])
        while True:
            try:
                tmp = pd.read_table(
                    fname,
                    usecols=[idx],
                    header=None,
                    nrows=chunksize,
                    skiprows=skiprows,
                )
            except pd.errors.EmptyDataError:
                break
            if tmp.shape[0] == 0 or tmp.iloc[0, 0] in enhancers:
                break

            skiprows += chunksize
            enhancers = np.hstack((enhancers, tmp.iloc[:, 0].unique()))
        enhancers = np.unique(enhancers)

        # Split into columns and create PyRanges object
        p = re.compile("[:-]")
        enhancers = pr.PyRanges(
            pd.DataFrame(
                [re.split(p, e) for e in enhancers],
                columns=["Chromosome", "Start", "End"],
            )
        )
        return enhancers

    @staticmethod
    def distance_weight(
        include_promoter=False,
        include_enhancer=True,
        alpha=1e4,
        maximum_distance=100_000,
        full_weight_region=5000,
        promoter_region=2000,
    ):
        """Build weight distribution based on distance to TSS.

        The basic idea is similar to Wang et al. [1], with some modifications.
        The resulting weight ranges from 0 (far from the TSS) to 1 (near the
        TSS) and is based on several different variables.

        If `include_promoter` is `True`, then distances smaller than
        `promoter_region` are included, otherwise they are excluded, the weight
        is set to 0.
        The `full_weight_region` parameters determines the region where
        the weight will be 1, regardless of distance. The `maximum_distance`
        parameter sets the maximum distance to consider. The weight decays with
        an increasing distance, starting from 1 at `full_weight_region` to 0
        at `maximum_distance`. The `alpha` parameters controls the decay.

        Parameters
        ----------
        include_promoter : bool, optional
            Include promoter regions. Default is False.
        include_enhancer : bool, optional
            Include enhancer regions, ie. regions that are distal to the
            promoter.
        alpha : float, optional
            Controls weight decay, default is 1e4.
        maximum_distance : int, optional
            Maximum distance from TSS to consider. Default is 100kb.
        full_weight_region : int, optional
            Distance where regions will receive the full weight. Default
            is 5kb.
        promoter_region : int, optional
            Promoter region, default is 2kb.

        Returns
        -------
        DataFrame with two columns: distance and weight.

        References
        ----------
        ..[1] Wang S, Zang C, Xiao T, Fan J, Mei S, Qin Q, Wu Q, Li X, Xu K,
        He HH, Brown M, Meyer CA, Liu XS. "Modeling cis-regulation with a
        compendium of genome-wide histone H3K27ac profiles." Genome Res.
        2016 Oct;26(10):1417-1429. doi: 10.1101/gr.201574.115. PMID: 27466232
        """
        u = -math.log(1.0 / 3.0) * 1e5 / alpha

        promoter_weight = int(include_promoter)
        enhancer_weight = int(include_enhancer)

        weight1 = pd.DataFrame(
            {
                "weight": [promoter_weight for _ in range(0, promoter_region + 1)],
                "dist": range(0, promoter_region + 1),
            }
        )

        weight2 = pd.DataFrame(
            {
                "weight": [
                    enhancer_weight
                    for _ in range(promoter_region + 1, full_weight_region + 1)
                ],
                "dist": range(promoter_region + 1, full_weight_region + 1),
            }
        )

        weight3 = pd.DataFrame(
            {
                "weight": [
                    enhancer_weight
                    * 2.0
                    * math.exp(-u * math.fabs(z) / 1e5)
                    / (1.0 + math.exp(-u * math.fabs(z) / 1e5))
                    for z in range(1, maximum_distance - full_weight_region + 1)
                ],
                "dist": range(full_weight_region + 1, maximum_distance + 1),
            }
        )

        weight = pd.concat([weight1, weight2, weight3])
        return weight

    def enhancer2gene(
        self,
        peak_pr,
        up=100_000,
        down=100_000,
        alpha=1e4,
        promoter=2000,
        full_weight_region=5000,
    ):
        """Couple enhancers to genes.

        Parameters
        ----------
        peak_pr : PyRanges object
            PyRanges object with enhancer regions.
        up : int, optional
            Upstream maximum distance, by default 100kb.
        down : int, optional
            Upstream maximum distabce, by default 100kb.
        alpha : float, optional
            Parameter to control weight decay, by default 1e4.
        promoter : int, optional
            Promoter region, by default 2000.
        full_weight_region : int, optional
            Region that will receive full weight, by default 5000.

        Returns
        -------
        pandas.DataFrame
            DataFrame with enhancer regions, gene names, distance and weight.
        """
        genes = region_gene_overlap(peak_pr, self.gene_bed)
        if genes.empty:
            return pd.DataFrame()

        # Get the distance from center of enhancer to TSS
        # Correct for extension
        genes["dist"] = (
            (genes["Start_b"] + genes["End_b"]) / 2 - genes["Start"]
        ).astype(int)
        genes.loc[genes["Strand"] == "+", "dist"] -= up
        genes.loc[genes["Strand"] == "-", "dist"] -= down
        genes["dist"] = np.abs(genes["dist"])

        # Create region in chr:start:end format
        genes["loc"] = (
            genes["Chromosome"].astype(str)
            + ":"
            + genes["Start_b"].astype(str)
            + "-"
            + genes["End_b"].astype(str)
        )

        # Keep the gene-enhancer combination with the smallest distance
        genes = genes.sort_values("dist").drop_duplicates(
            subset=["loc", "Name"], keep="first"
        )

        # Return the right stuff
        genes = genes.set_index("loc")[["Name", "dist"]].rename(
            columns={"Name": "gene"}
        )

        # Get distance-based wight
        weight = self.distance_weight(
            include_promoter=self.include_promoter,
            include_enhancer=self.include_enhancer,
            alpha=alpha,
            promoter_region=promoter,
            full_weight_region=full_weight_region,
        ).set_index("dist")
        genes = genes.join(weight, on="dist")

        return genes

    def aggregate_binding(
        self,
        binding_fname,
        tfs=None,
        up=1e5,
        down=1e5,
        alpha=None,
        promoter=2000,
        full_weight_region=5000,
        combine_function="sum",
    ):
        """Summarize all binding signal per gene per TF.

        Return a dask delayed computation object.

        Parameters
        ----------
        binding_fname : str
            Filename of binding network.
        tfs : list, optional
            List of transcription factor names, by default None, which means
            that all TFs will be used.
        up : int, optional
            Maximum upstream region to include, by default 1e5
        down : [type], optional
            Maximum downstream region to include, by default 1e5
        alpha : float, optional
            Distance at which the weight will be half, by default None
        promoter : int, optional
            Promoter region, by default 2000
        full_weight_region : int, optional
            Region that will receive full weight, regardless of distance, by
            default 5000.
        combine_function : str, optional
            How to combine signal of weighted enhancers, by default "sum".
            Valid options are "sum", "mean" or "max".

        Returns
        -------
        dask.DataFrame
            DataFrame with delayed computations.
        """
        if not os.path.exists(binding_fname):
            raise ValueError(f"File {binding_fname} does not exist!")

        if combine_function not in ["mean", "max", "sum"]:
            raise NotImplementedError(
                "Unknown combine function, valid options are: mean, max, sum"
            )

        maximum_distance = max(up, down)
        if alpha is None:
            alpha = maximum_distance / 10

        if promoter > maximum_distance:
            raise ValueError(
                "promoter region is larger than the maximum distance to use"
            )

        hdf = HDFStore(binding_fname, "r")

        # TODO: This is hacky (depending on "_"), however the hdf.keys() method is
        # much slower. Currently all TF names do *not* start with "_"
        all_tfs = [x for x in dir(hdf.root) if not x.startswith("_")]
        logger.info(f"Binding file contains {len(all_tfs)} TFs.")
        if tfs is None:
            tfs = all_tfs
        else:
            not_valid = set(all_tfs) - set(tfs)
            if len(not_valid) > 1:
                logger.warning(
                    f"The following TFs are found in {binding_fname}, but do not seem to be TFs:"
                )
                logger.warning(", ".join(not_valid))
            tfs = set(tfs) & set(all_tfs)
            logger.info(f"Using {len(tfs)} TFs.")

        # Read enhancer index from hdf5 file
        enhancers = hdf.get(key="_index")
        chroms = enhancers.index.to_series().str.replace(":.*", "").unique()

        tmpdir = mkdtemp()
        self._tmp_files.append(tmpdir)  # mark for deletion later

        # Summarize enhancers per gene, per chromosome. In principle this could
        # also be done at once, however, the memory usage of dask is very finicky.
        # This is a pragmatic solution, that seems to work well, does not use a
        # lot of memory and is not too slow (~50 seconds per chromosome).
        for chrom in chroms:
            logger.info(f"Aggregating binding for genes on {chrom}")

            # Get the index of all enhancers for this specific chromosome
            idx = enhancers.index.str.contains(f"^{chrom}:")
            idx_i = np.arange(enhancers.shape[0])[idx]

            # Create a pyranges object
            enhancer_pr = pr.PyRanges(
                enhancers[idx]
                .index.to_series()
                .str.split(r"[:-]", expand=True)
                .rename(columns={0: "Chromosome", 1: "Start", 2: "End"})
            )

            # Link enhancers to genes on basis of distance to annotated TSS
            gene_df = self.enhancer2gene(
                enhancer_pr,
                up=up,
                down=down,
                alpha=alpha,
                promoter=promoter,
                full_weight_region=full_weight_region,
            )
            gene_df.dropna(inplace=True)
            if gene_df.empty:
                logger.debug(f"No genes found on {chrom}")
                continue

            bp = pd.DataFrame(index=enhancers[idx].index)

            for tf in tqdm(
                tfs, total=len(tfs), desc="Aggregating", unit_scale=1, unit=" TFs"
            ):
                # Load TF binding data for this chromosome.
                # hdf.get() is *much* faster here than pd.read_hdf()
                bp[tf] = hdf.get(key=tf)[idx_i].values

            # Skipping everything with weight 0, as it won't be counted anyway.
            gene_df = gene_df[gene_df["weight"] > 0]

            # Make sure binding score and enhancers match up (i.e. same enhancer
            # is used for multiple genes)
            gene_df = gene_df.join(bp).dropna()
            bp = gene_df[tfs]
            gene_df = gene_df[["gene", "weight"]]

            # Multiply binding score by weight
            bp = bp.mul(gene_df["weight"], axis=0)

            # Summarize weighted score per gene
            bp["gene"] = gene_df["gene"]
            tmp = bp.groupby("gene")
            if combine_function == "mean":
                tmp = tmp.mean()
            elif combine_function == "max":
                tmp = tmp.max()
            elif combine_function == "sum":
                tmp = tmp.sum()

            # Go from wide to long format, to be able to merge with other
            # information later
            tmp = tmp.reset_index().melt(
                id_vars=tmp.index.name, var_name="tf", value_name="weighted_binding"
            )

            # Create dataframe with two columns: tf_gene and weighted_binding score
            tmp["tf_target"] = tmp["tf"] + SEPARATOR + tmp["gene"]
            tmp[["tf_target", "weighted_binding"]].to_csv(
                os.path.join(tmpdir, f"{chrom}.csv"), index=False
            )

        hdf.close()

        ddf = dd.read_csv(os.path.join(tmpdir, "*.csv")).set_index("tf_target")
        return ddf

    def _save_temp_expression(self, df, name):
        tmp = df.rename(columns={"tpm": f"{name}_expression"})
        tmp[f"{name}_expression"] = minmax_scale(tmp[f"{name}_expression"].rank())
        tmp.index.rename(name, inplace=True)
        tmp["key"] = 0
        fname = NamedTemporaryFile(
            prefix="ananse.", suffix=f".{name}.parquet", delete=False
        ).name
        self._tmp_files.append(fname)
        tmp.reset_index().to_parquet(fname, index=False)
        return fname

    def create_expression_network(self, expression, tfs, column="tpm"):
        """Create a gene expression based network.

        Based on file(s) with gene expression levels (a TPM column), a
        dask DataFrame is generated with the combined expression levels
        of the tf and the target gene. By default, the expression levels
        are ranked and subsequently scaled between 0 and 1.

        Parameters
        ----------
        expression : pd.DataFrame
            gene expression data.
            First column should contain the gene names.

        tfs : list
            List of TF gene names. All TFs will be used by default.

        column : str, optional
            Column name that contains gene expression, 'tpm' by default (case insensitive).

        Returns
        -------
            Dask DataFrame with gene expression based values.
        """
        # log transform expression values
        expression[column] = np.log2(expression[column] + 1e-5)

        tf_expr = expression[expression.index.isin(tfs)]
        tf_fname = self._save_temp_expression(tf_expr, "tf")
        target_fname = self._save_temp_expression(expression, "target")
        # Read files (delayed) and merge on 'key' to create a Cartesian product
        # combining all TFs with all target genes.
        a = dd.read_parquet(tf_fname)
        b = dd.read_parquet(target_fname)
        network = a.merge(b, how="outer")

        # Use one-column index that contains TF and target genes.
        # This is necessary for dask, as dask cannot merge on a MultiIndex.
        # Otherwise this would be an inefficient and unnecessary step.
        network["tf_target"] = network["tf"] + SEPARATOR + network["target"]
        network = network[
            ["tf", "target", "tf_target", "tf_expression", "target_expression"]
        ]

        return network

    def run_network(
        self,
        binding,
        fin_expression=None,
        tfs=None,
        outfile=None,
        up=1e5,
        down=1e5,
        alpha=None,
        promoter=2000,
        full_weight_region=5000,
    ):

        """Create network.

        Generates a gene-regulatory network with a TF-target gene interaction "prob" score based on the mean rank of:
        1. Binding score based on ATAC/H3K27ac data of nearby enhancers.
        2. TF Activity (single score per TF based on general TF motif behaviour in the trained dataset)
        3. TF expression score
        4. Target expression score

        Parameters
        ----------
        binding : str
            Filename with binding information. Should contain at least three
            columns: "factor", "enhancer" and "binding".
        fin_expression : str or list, optional
            Filename of list of filenames with expression information.
        tfs : list, optional
            List of transcription factors to use, by default None, which means
            all TFs will be used.
        outfile : str, optional
            Output file. If None, returns a dataframe.
        up : int, optional
            Upstream maximum distance, by default 100kb.
        down : int, optional
            Upstream maximum distabce, by default 100kb.
        alpha : float, optional
            Parameter to control weight decay, by default 1e4.
        promoter : int, optional
            Promoter region, by default 2000.
        full_weight_region : int, optional
            Region that will receive full weight, by default 5000.
        """
        column = "tpm"  # TODO: expose as argument

        # get all TFs from binding (or database),
        # then intersect with given TFs (if any)
        # COMPATIBILITY: validates binding-tf (or database-tf) overlap
        tfs = get_factors(binding, tfs)

        # expression dataframe with transcription factors as index
        expression = combine_expression_files(fin_expression, column)

        # check for sufficient overlap in gene/transcript names/identifiers
        # attempt to fix issues if a genomepy assembly was used
        # COMPATIBILITY: validates/fixes tf-expression and tf-gene_bed overlap
        expression = self.gene_overlap(expression, tfs)

        # Expression base network
        df_expression = self.create_expression_network(expression, tfs, column)

        # Use a version of the binding network:
        # either promoter-based, enhancer-based or both.
        if self.include_promoter or self.include_enhancer:
            df_binding = self.aggregate_binding(
                binding,
                tfs=tfs,
                up=up,
                down=down,
                alpha=alpha,
                promoter=promoter,
                full_weight_region=full_weight_region,
                combine_function="sum",
            )

            try:
                act = pd.read_hdf(binding, key="_factor_activity")
                if "factor" in act.columns:  # noqa
                    act = act.set_index("factor")  # noqa
                logger.info("Reading factor activity")
                act.index.name = "tf"
                act["activity"] = minmax_scale(rankdata(act["activity"], method="min"))
                df_expression = df_expression.merge(
                    act, right_index=True, left_on="tf", how="left"
                ).fillna(0.5)
            except KeyError:
                pass

            df_expression = df_expression.drop(columns=["tf"])

            # This is where the heavy lifting of all delayed computations gets done
            if fin_expression is not None:
                result = df_expression.merge(
                    df_binding, right_index=True, left_on="tf_target", how="left"
                )
                result = result.persist()
                result = result.fillna(0)
                logger.info("Processing expression+binding network")
                progress(result)
                result = result.compute()
            else:
                result = df_binding

            result["weighted_binding"] = minmax_scale(
                rankdata(result["weighted_binding"], method="min")
            )
            columns = [
                "tf_expression",
                "target_expression",
                "weighted_binding",
                "activity",
            ]
            columns = [col for col in columns if col in result]
            logger.info(f"Using {', '.join(columns)}")
            # Combine the individual scores
            result["prob"] = result[columns].mean(1)

        else:
            result = df_expression
            result["prob"] = result[["tf_expression", "target_expression"]].mean(1)
            logger.info("Processing expression network")
            result = result.compute()

        output_cols = ["tf_target", "prob"]
        if self.full_output:
            columns = [
                "tf_target",
                "prob",
                "tf_expression",
                "target_expression",
                "weighted_binding",
                "activity",
            ]
            output_cols = [col for col in columns if col in result]
        if outfile:
            logger.info("Writing network")
            out_dir = os.path.abspath(os.path.dirname(outfile))
            os.makedirs(out_dir, exist_ok=True)
            result[output_cols].to_csv(outfile, sep="\t", index=False)
        else:
            return result[output_cols]

    def gene_overlap(self, expression: pd.DataFrame, tfs: list):
        """
        For ANANSE Network to run properly, we need overlap between 3 files/sets with gene names:
        - TFs (from motif2factors.txt/binding file)
        - expression file(s)
        - gene annotation BED file

        if the overlap is low, attempt to convert the names to those of the TFs with genomepy.
        """
        cutoff = 0.6  # fraction of overlap that is "good enough"
        tfs = set(tfs)
        gp = genomepy.Annotation(self.gene_bed)
        bed_genes = gp.genes()

        overlapping_tfs = set(expression.index).intersection(tfs)
        overlap_tf_exp = len(overlapping_tfs) / len(tfs)
        logger.debug(
            f"{int(100 * overlap_tf_exp)}% of TFs found in the expression file(s)"
        )

        overlapping_tfs = set(bed_genes).intersection(tfs)
        overlap_tf_bed = len(overlapping_tfs) / len(tfs)
        logger.debug(f"{int(100 * overlap_tf_bed)}% of TFs found in the BED file")

        overlapping_tfs = tfs.intersection(set(expression.index), set(bed_genes))
        overlap_total = len(overlapping_tfs) / len(tfs)
        logger.debug(
            f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
        )

        if overlap_total > cutoff:
            logger.info(
                f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
            )
            return expression
        if gp.annotation_gtf_file is None:
            if overlap_total == 0:
                incompatible_gene_error()
            logger.warning(
                "Little overlap between genes! "
                "Are TFs, expression and gene_bed using the same symbols?"
            )
            logger.info(
                f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
            )
            return expression

        logger.warning(
            "Little overlap between genes! "
            "Converting genes in expression table and BED to HGNC symbols"
        )
        # assumption: you used gimme motif2factors on the GTF file of this genome
        tid2gid = gp.gtf_dict("transcript_id", "gene_id")
        tid2name = gp.gtf_dict("transcript_id", "gene_name")
        gid2name = gp.gtf_dict("gene_id", "gene_name")

        if overlap_tf_exp <= cutoff:
            expression = (
                expression.rename(index=tid2name)
                .rename(index=tid2gid)
                .rename(index=gid2name)
            )

            # metrics
            overlapping_tfs = set(expression.index).intersection(tfs)
            overlap_tf_exp = len(overlapping_tfs) / len(tfs)
            logger.debug(
                f"{int(100 * overlap_tf_exp)}% of TFs found in the expression file(s)"
            )

            # merge duplicate genes
            expression = expression.groupby(by=expression.index).sum()

        if overlap_tf_bed <= cutoff:
            bed = (
                gp.bed.copy()
                .set_index("name")
                .rename(index=tid2name)
                .rename(index=tid2gid)
                .rename(index=gid2name)
                .reset_index()
            )

            # metrics
            bed_genes = bed.name
            overlapping_tfs = set(bed_genes).intersection(tfs)
            overlap_tf_bed = len(overlapping_tfs) / len(tfs)
            logger.debug(f"{int(100 * overlap_tf_bed)}% of TFs found in the BED file")

            # merge duplicate genes
            group = bed.groupby("name")
            bed["start"] = group["start"].transform("min")
            bed["end"] = group["end"].transform("max")
            cols = set(bed.columns) - {"name", "chrom", "start", "end", "strand"}
            for col in cols:
                bed[col] = 0
            bed.drop_duplicates(inplace=True, ignore_index=True)
            tpm_bed = NamedTemporaryFile(
                prefix="ananse.", suffix=f".annotation.bed", delete=False
            ).name
            genomepy.annotation.utils.write_annot(bed, tpm_bed)
            self._tmp_files.append(tpm_bed)
            self.gene_bed = tpm_bed

        overlapping_tfs = tfs.intersection(set(expression.index), set(bed_genes))
        overlap_total = len(overlapping_tfs) / len(tfs)
        # logger.debug(
        #     f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
        # )

        if overlap_total > 0:
            if overlap_total <= cutoff:
                logger.warning(
                    "Little overlap between genes! "
                    "Are TFs, expression and gene_bed using the same symbols?"
                )
            logger.info(
                f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
            )
            return expression
        incompatible_gene_error()

    def __del__(self):
        if not hasattr(self, "_tmp_files"):
            return

        for fname in self._tmp_files:
            if os.path.exists(fname):
                shutil.rmtree(fname, ignore_errors=True)


def incompatible_gene_error():
    msg = [
        "None of the transcription factors are found in your expression file(s)/BED file.",
        "If you have human data, please make sure you use HGNC symbols (gene names).",
        "If you have non-human data, you have to create a custom motif to gene mapping.",
        "See this link for one possibility to create this file: ",
        "https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors",
        "If you use a custom motif mapping, you will also have (re)run `gimme binding` with this file.",
    ]
    for line in msg:
        logger.error(line)
    sys.exit(1)


def get_bed(gene_bed, genome):
    out_bed = gene_bed
    if out_bed is None:
        if genome in ["hg38", "hg19"]:
            out_bed = os.path.join(PACKAGE_DIR, "db", f"{genome}.genes.bed")
        else:
            gp = genomepy.Genome(genome)  # can raise descriptive FileNotFoundError
            out_bed = gp.annotation_bed_file
    if out_bed is None:
        raise TypeError("Please provide a gene bed file with the -a argument.")
    if not os.path.exists(out_bed):
        raise FileNotFoundError(f"Could not find gene bed file {out_bed}.")
    return out_bed


def get_factors(bindingfile: str, tfs=None):
    """
    Return a list of transcription factors in the bindingfile or the database.
    If TFs is given, this is used to filter the output list.

    Returns
    -------
    list
        of unique transcription factors
    """
    # Create the TF list, based on valid transcription factors
    try:
        act = pd.read_hdf(bindingfile, key="_factor_activity")
        if "factor" in act.columns:  # noqa
            act = act.set_index("factor")  # noqa
        out_tfs = act.index
    except KeyError:
        tffile = os.path.join(PACKAGE_DIR, "db", "tfs.txt")
        out_tfs = pd.read_csv(tffile, header=None)[0]
        logger.warning("No TFs found in binding.h5 file. Using database file")
    if tfs is not None:
        out_tfs = set(out_tfs).intersection(set(tfs))
        if len(out_tfs) == 0:
            raise ValueError(
                "No transcription factor overlap between given TFs and binding file"
            )
    return list(set(out_tfs))


def region_gene_overlap(
    region_pr,
    gene_bed,
    up=100_000,
    down=100_000,
):
    """
    Couple enhancers to genes.

    Parameters
    ----------
    region_pr : PyRanges object
        PyRanges object with enhancer regions.
    gene_bed : str
        gene_bed
    up : int, optional
        Upstream maximum distance, by default 100kb.
    down : int, optional
        Upstream maximum distance, by default 100kb.

    Returns
    -------
    pandas.DataFrame
        DataFrame with enhancer regions, gene names, distance and weight.
    """
    # TODO: Ensembl chrom MT interpreted as number (ONLY WITH MY CONVERTED BED)
    genes = genomepy.Annotation(gene_bed).bed  # pr.read_bed(gene_bed)
    genes.columns = [col.capitalize() for col in genes.columns]
    # Convert to DataFrame & we don't need intron/exon information
    genes = genes.as_df().iloc[:, :6]

    # Get the TSS only
    genes.loc[genes["Strand"] == "+", "End"] = genes.loc[
        genes["Strand"] == "+", "Start"
    ]
    genes.loc[genes["Strand"] == "-", "Start"] = genes.loc[
        genes["Strand"] == "-", "End"
    ]

    # Extend up and down
    genes.loc[genes["Strand"] == "+", "Start"] -= up
    genes.loc[genes["Strand"] == "+", "End"] += down
    genes.loc[genes["Strand"] == "-", "Start"] -= down
    genes.loc[genes["Strand"] == "-", "End"] += up

    # Perform the overlap
    genes = pr.PyRanges(genes)
    genes = genes.join(region_pr).as_df()

    return genes


def combine_expression_files(fin_expression: Union[str, list], column="tpm"):
    """
    Extract the index and expression column from one or more expression files.
    We expect the index to be gene/transcript names/identifiers, and the expression to be TPMs.

    Within each expression file, duplicate indexes are summed.
    Between expression files, indexes are averaged.
    NAs are dropped.

    Parameters
    ----------
    fin_expression : str or list
        One of more files that contains gene expression data.
        First column should contain the gene/transcript names/identifiers.

    column : str, optional
        Column name that contains gene expression, 'tpm' by default (case insensitive).

    Returns
    -------
    pd.DataFrame
        a dataframe with one expression column and all unique indexes
    """
    logger.info("Loading expression")

    # Convert to a list of filename(s)
    if isinstance(fin_expression, str):
        fin_expression = [fin_expression]

    # Read all expression input files and take the mean expression per gene
    re_column = re.compile(fr"^{column}$", re.IGNORECASE)
    expression = pd.DataFrame()
    for f in fin_expression:
        # keep only the name and expression columns
        subdf = pd.read_table(f, index_col=0).filter(regex=re_column)
        # combine expression values for duplicate gene names
        # also removes NaNs from the index
        subdf = subdf.groupby(by=subdf.index, dropna=True).sum()
        subdf.dropna(inplace=True)
        expression = pd.concat([expression, subdf], axis=1)
    expression = pd.DataFrame(expression.mean(1), columns=[column])
    return expression
