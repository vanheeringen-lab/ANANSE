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
from tqdm.auto import tqdm
import pyranges as pr

from ananse.utils import cleanpath
from ananse.view import get_binding_tfs
from . import SEPARATOR, PACKAGE_DIR


warnings.filterwarnings("ignore")


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

    # @staticmethod
    # def unique_enhancers(fname):
    #     """Extract a list of unique enhancers.
    #
    #     Parameters
    #     ----------
    #     fname : str
    #         File name of a tab-separated file that contains an 'enhancer' column.
    #
    #     Returns
    #     -------
    #         PyRanges object with enhancers
    #     """
    #     logger.info("reading enhancers")
    #
    #     # Read enhancer regions from binding file
    #     header = pd.read_table(fname, nrows=0)
    #     idx = header.columns.get_loc("enhancer")
    #     skiprows = 1
    #     chunksize = 2_000_000
    #     enhancers = np.array([])
    #     while True:
    #         try:
    #             tmp = pd.read_table(
    #                 fname,
    #                 usecols=[idx],
    #                 header=None,
    #                 nrows=chunksize,
    #                 skiprows=skiprows,
    #             )
    #         except pd.errors.EmptyDataError:
    #             break
    #         if tmp.shape[0] == 0 or tmp.iloc[0, 0] in enhancers:
    #             break
    #
    #         skiprows += chunksize
    #         enhancers = np.hstack((enhancers, tmp.iloc[:, 0].unique()))
    #     enhancers = np.unique(enhancers)
    #
    #     # Split into columns and create PyRanges object
    #     p = re.compile("[:-]")
    #     enhancers = pr.PyRanges(
    #         pd.DataFrame(
    #             [re.split(p, e) for e in enhancers],
    #             columns=["Chromosome", "Start", "End"],
    #         )
    #     )
    #     return enhancers

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
        genes = region_gene_overlap(peak_pr, self.gene_bed, up, down)
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
        binding,
        tfs=None,
        regions=None,
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
        binding : str
            Filename of binding network.
        tfs : list, optional
            List of transcription factor names, by default None, which means
            that all TFs will be used.
        up : int, optional
            Maximum upstream region to include, by default 1e5
        regions : list, optional
            List of regions to aggregate on.
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
        logger.info("Loading binding data")

        if not os.path.exists(cleanpath(binding)):
            raise ValueError(f"File {binding} does not exist!")

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

        # Read TFs
        all_tfs = get_binding_tfs(binding)
        if tfs is None:
            tfs = all_tfs
        if len(tfs) == 0:
            raise ValueError(
                "No transcription factor overlap between requested TFs and binding file! "
                f"Use `ananse view --list-tfs {binding}` to inspect the TFs in the file."
            )
        if len(tfs) == len(all_tfs):
            logger.info(f"Using all {len(set(tfs))} TFs.")
        else:
            logger.info(f"Using {len(tfs)} of {len(set(all_tfs))} TFs.")

        hdf = pd.HDFStore(binding, "r")

        # Read enhancer regions
        enhancers = hdf.get(key="_index")
        chroms = set(enhancers.index.str.replace(":.*", ""))
        if regions is None:
            logger.info(f"Using all {len(set(enhancers.index))} regions.")
        else:
            regions = set(regions) & set(enhancers.index)
            chroms = [region.split(":")[0] for region in regions]
            logger.info(f"Using {len(regions)} of {len(set(enhancers.index))} regions.")
        if len(chroms) == 0:
            raise ValueError(
                "No regions in the binding file overlap with given regions! "
                f"Use `ananse view --list-regions {binding}` to inspect the regions in the file."
            )

        tmpdir = mkdtemp()
        self._tmp_files.append(tmpdir)  # mark for deletion later

        # Summarize enhancers per gene, per chromosome. In principle this could
        # also be done at once, however, the memory usage of dask is very finicky.
        # This is a pragmatic solution, that seems to work well, does not use a
        # lot of memory and is not too slow (~50 seconds per chromosome).
        for chrom in chroms:
            # Get the index of all enhancers for this specific chromosome
            idx = enhancers.index.str.startswith(f"{chrom}:")
            if regions:
                region_idx = enhancers.index.isin(regions)
                idx = list(np.array(idx) * np.array(region_idx))

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
                tfs,
                total=len(tfs),
                desc=f"Aggregating on {chrom}",
                unit_scale=1,
                unit=" TFs",
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
        binding=None,
        fin_expression=None,
        tfs=None,
        regions=None,
        outfile=None,
        column="tpm",
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
        binding : str, optional
            Filename with binding information. Should contain at least four
            columns: "factor", "enhancer", "binding" and "_factor_activity".
            Required if no expression file(s) are provided.
        fin_expression : str or list, optional
            Filename of list of filenames with expression information.
            Required if no binding file is provided.
        tfs : list, optional
            List of transcription factors to use, by default None, which means
            all TFs in the binding file will be used.
            Required if no binding file is provided.
        regions : list, optional
            List of regions to limit the binding network to.
        outfile : str, optional
            Output file. If None, returns a dataframe.
        column : string, optional
            Name of the column containing the expression data in the expression file(s).
            Defaults to 'tpm' (case insensitive).
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
        # get all TFs from binding (or database),
        # then intersect with given TFs (if any)
        # COMPATIBILITY: validates binding-tf (or database-tf) overlap
        tfs = get_factors(binding, tfs)

        # create the expression network
        df_expression = None
        if fin_expression is not None:
            # expression dataframe with transcription factors as index
            expression = combine_expression_files(fin_expression, column)

            # check for sufficient overlap in gene/transcript names/identifiers
            # attempt to fix issues if a genomepy assembly was used
            # COMPATIBILITY: validates/fixes tf-expression and tf-gene_bed overlap
            expression = self.gene_overlap(expression, tfs)

            df_expression = self.create_expression_network(expression, tfs, column)
            if binding is not None:
                logger.debug("Loading tf binding activity data")
                act = pd.read_hdf(binding, key="_factor_activity")
                act = act.set_index("factor")
                act.index.name = "tf"
                act["activity"] = minmax_scale(rankdata(act["activity"], method="min"))
                df_expression = df_expression.merge(
                    act, right_index=True, left_on="tf", how="left"
                ).fillna(0.5)
            df_expression = df_expression.drop(columns=["tf"])

        # create the binding network
        # promoter-based, enhancer-based or both.
        df_binding = None
        if binding is not None and (self.include_promoter or self.include_enhancer):
            df_binding = self.aggregate_binding(
                binding,
                tfs=tfs,
                regions=regions,
                up=up,
                down=down,
                alpha=alpha,
                promoter=promoter,
                full_weight_region=full_weight_region,
                combine_function="sum",
            )

        # (combine networks and) compute delayed operations
        if df_binding is None:
            if df_expression is None:
                raise ValueError(
                    "Networks are based on at least one of "
                    "expression, promoter binding or enhancer binding!"
                )
            else:
                logger.info("Processing expression network")
                result = df_expression
        else:
            if df_expression is None:
                logger.info("Processing binding network")
                result = df_binding
            else:
                logger.info("Processing expression-binding network")
                result = df_expression.merge(
                    df_binding, right_index=True, left_on="tf_target", how="left"
                )
                result = result.persist()
                result = result.fillna(0)
                # TODO: move client from commands.network into this module
                if "client" in globals():
                    progress(result)
        # This is where the heavy lifting of all delayed computations gets done
        result = result.compute()

        if "weighted_binding" in result:
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
        logger.debug(f"Using {', '.join(columns)}")
        # Combine the individual scores
        result["prob"] = result[columns].mean(1)

        # filter output
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

        overlapping_tfs = set(expression.index) & tfs
        overlap_tf_exp = len(overlapping_tfs) / len(tfs)
        logger.debug(
            f"{int(100 * overlap_tf_exp)}% of TFs found in the expression file(s)"
        )

        overlapping_tfs = set(bed_genes) & tfs
        overlap_tf_bed = len(overlapping_tfs) / len(tfs)
        logger.debug(f"{int(100 * overlap_tf_bed)}% of TFs found in the BED file")

        overlapping_tfs = tfs.intersection(set(expression.index), set(bed_genes))
        overlap_total = len(overlapping_tfs) / len(tfs)
        logger.info(
            f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
        )

        if overlap_total > cutoff:
            return expression
        if gp.annotation_gtf_file is None:
            if overlap_total == 0:
                incompatible_gene_error()
            logger.warning("Are TFs, expression and gene_bed using the same symbols?")
            return expression

        logger.warning("Converting genes in expression table and BED to HGNC symbols")
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
            overlapping_tfs = set(expression.index) & tfs
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
            overlapping_tfs = set(bed_genes) & tfs
            overlap_tf_bed = len(overlapping_tfs) / len(tfs)
            logger.debug(f"{int(100 * overlap_tf_bed)}% of TFs found in the BED file")

            # merge duplicate genes
            group = bed.groupby("name")
            bed["start"] = group["start"].transform("min")
            bed["end"] = group["end"].transform("max")
            drop_cols = set(bed.columns) - {"name", "chrom", "start", "end", "strand"}
            for col in drop_cols:
                bed[col] = 0
            bed.drop_duplicates(inplace=True, ignore_index=True)
            tpm_bed = NamedTemporaryFile(
                prefix=f"ananse.{gp.name}", suffix=".annotation.bed", delete=False
            ).name
            cols = genomepy.annotation.utils.BED12_FORMAT  # fixes column order
            genomepy.annotation.utils.write_annot(bed[cols], tpm_bed)
            self._tmp_files.append(tpm_bed)
            self.gene_bed = tpm_bed

        overlapping_tfs = tfs.intersection(set(expression.index), set(bed_genes))
        overlap_total = len(overlapping_tfs) / len(tfs)
        logger.info(
            f"{int(100 * overlap_total)}% of TFs found in both BED and expression file(s)"
        )

        if overlap_total > 0:
            if overlap_total <= cutoff:
                logger.warning(
                    "Are TFs, expression and gene_bed using the same symbols?"
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
        "If you have human data, please make sure you use HGNC symbols (gene names) in the BED and expression file(s).",
        "If you have non-human data, you have to create a custom motif to gene mapping.",
        "  See this link for one possibility to create this file: ",
        "  https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors",
        "If you use a custom motif mapping, you will also have (re)run `ananse binding` with this file.",
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


def get_factors(binding: str = None, tfs: list = None):
    """
    Return a list of transcription factors in the binding or the database.
    If TFs is given, this is used to filter the output list.

    Returns
    -------
    list
        of unique transcription factors
    """
    if binding is None and tfs is None:
        raise ValueError(
            "A binding file, a list of transcription factors, or both is required!"
        )

    out_tfs = tfs
    if binding is not None:
        out_tfs = get_binding_tfs(binding)
        if tfs is not None:
            not_valid = set(tfs) - set(out_tfs)
            if len(not_valid) > 0:
                logger.warning(
                    f"The following TFs are requested, but not found in {binding}:"
                )
                logger.warning(", ".join(not_valid))
                logger.warning(
                    "Perhaps an error occurred for these TFs in ananse binding."
                )
            out_tfs = set(tfs) & set(out_tfs)

    if len(out_tfs) == 0:
        raise ValueError(
            "No transcription factor overlap between requested TFs and binding file! "
            f"Use `ananse view --list-tfs {binding}` to inspect the TFs in the file."
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
    genes = pr.read_bed(gene_bed)
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
    logger.info("Loading expression data")

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
