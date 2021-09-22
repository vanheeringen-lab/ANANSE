#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

"""Build gene regulatory network"""

# Python imports
import os
import math
import re
import shutil
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
import dask.dataframe as dd
from tempfile import NamedTemporaryFile, mkdtemp
from dask.distributed import progress
from loguru import logger
from pandas import HDFStore
from tqdm.auto import tqdm

import pyranges as pr

warnings.filterwarnings("ignore")
PACKAGE_DIR = os.path.dirname(__file__)


class Network(object):
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
        self._tmp_files = []

        # # Motif information file
        # if pfmfile is None:
        #     self.pfmfile = "../data/gimme.vertebrate.v5.1.pfm"
        # else:
        #     self.pfmfile = pfmfile

        # self.motifs2factors = self.pfmfile.replace(".pfm", ".motif2factors.txt")
        # self.factortable = self.pfmfile.replace(".pfm", ".factortable.txt")

        # Gene information file
        self.gene_bed = gene_bed
        if gene_bed is None:
            if self.genome in ["hg38", "hg19"]:
                self.gene_bed = os.path.join(
                    PACKAGE_DIR, "db", f"{self.genome}.genes.bed"
                )
            else:
                raise TypeError("Please provide a gene bed file with -a argument.")
        if not os.path.exists(self.gene_bed):
            raise FileNotFoundError(
                f"Could not find the gene bed file {self.gene_bed}."
            )

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
            tmp["tf_target"] = tmp["tf"] + "_" + tmp["gene"]
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

    def create_expression_network(
        self, fin_expression, column="tpm", tfs=None, bindingfile=None
    ):
        """Create a gene expression based network.

        Based on file(s) with gene expression levels (a TPM column), a
        dask DataFrame is generated with the combined expression levels
        of the tf and the target gene. By default, the expression levels
        are ranked and subsequently scaled between 0 and 1.

        Parameters
        ----------
        fin_expression : str or list
            One of more files that contains gene expression data.
            First column should contain the gene names in HGNC symbols.

        column : str, optional
            Column name that contains gene expression, 'tpm' by default (case insensitive).

        tfs : list, optional
            List of TF gene names. All TFs will be used by default.

        bindingfile : str, optional
            Output file from ANANSE binding.

        Returns
        -------
            Dask DataFrame with gene expression based values.
        """
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
        # log transform expression values
        expression[column] = np.log2(expression[column] + 1e-5)

        genes = pd.read_table(
            self.gene_bed, usecols=[3], comment="#", names=["name"], index_col=0
        )
        overlap = len(genes.index.intersection(expression.index))
        if overlap / expression.shape[0] < 0.1:
            logger.error(
                "gene annotation identifiers do not seem to match between annotation and expression files!"
            )
            sample_exp = ", ".join(expression.sample(5).index.values)
            sample_gene = ", ".join(genes.sample(5).index.values)
            logger.error(f"expression sample: {sample_exp}")
            logger.error(f"annotation sample: {sample_gene}")
            sys.exit(1)

        # Create the TF list, based on valid transcription factors
        if tfs is None:
            try:
                act = pd.read_hdf(bindingfile, key="_factor_activity")
                if "factor" in act.columns:
                    act = act.set_index("factor")
                tfs = list(set(act.index.tolist()))
            except KeyError:
                tffile = os.path.join(PACKAGE_DIR, "db", "tfs.txt")
                tfs = pd.read_csv(tffile, header=None)[0].tolist()

        # Save TFs and targets as temporary files
        idx = expression.index[expression.index.isin(tfs)]
        tmp = expression.loc[idx]
        if tmp.shape[0] == 0:
            logger.error(
                "None of the transcription factors are found in your expression file."
            )
            logger.error(
                "If you have human data, please make sure you use HGNC symbols (gene names)."
            )
            logger.error(
                "If you have non-human data, you have to create a custom motif to gene mapping."
            )
            logger.error("See this link for one possibility to create this file: ")
            logger.error(
                "https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors"
            )
            logger.error(
                "If you use a custom motif mapping, you will also have (re)run `gimme binding` with this file."
            )
            sys.exit(1)

        tf_fname = self._save_temp_expression(tmp, "tf")
        target_fname = self._save_temp_expression(expression, "target")
        # Read files (delayed) and merge on 'key' to create a Cartesian product
        # combining all TFs with all target genes.
        a = dd.read_parquet(tf_fname)
        b = dd.read_parquet(target_fname)
        network = a.merge(b, how="outer")

        # Use one-column index that contains TF and target genes.
        # This is necessary for dask, as dask cannot merge on a MultiIndex.
        # Otherwise this would be an inefficient and unnecessary step.
        network["tf_target"] = network["tf"] + "_" + network["target"]
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
        full_output=False,
    ):

        """Create network.

        Generates a gene-regulatory network with a TF-target gene interaction "prob" score based on the mean rank of:
        1. Binding score based on ATAC/H3K27ac data of neirby enhancers.
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
        # Expression base network
        logger.info("Loading expression")
        df_expression = self.create_expression_network(
            fin_expression, tfs=tfs, bindingfile=binding
        )

        # Use a version of the binding network, either promoter-based, enhancer-based
        # or both.
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
                if "factor" in act.columns:
                    act = act.set_index("factor")
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
            # logger.info("Computing network")
            if fin_expression is not None:
                result = df_expression.merge(
                    df_binding, right_index=True, left_on="tf_target", how="left"
                )
                result = result.persist()
                result = result.fillna(0)
                logger.info("Computing network")
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
            result = result.compute()

        if outfile:
            logger.info("Writing network")
            out_dir = os.path.abspath(os.path.dirname(outfile))
            os.makedirs(out_dir, exist_ok=True)
            if self.full_output:
                result[
                    [
                        "tf_target",
                        "prob",
                        "tf_expression",
                        "target_expression",
                        "weighted_binding",
                        "activity",
                    ]
                ].to_csv(outfile, sep="\t", index=False)
            else:
                result[["tf_target", "prob"]].to_csv(outfile, sep="\t", index=False)

        else:
            if self.full_output:
                return result[
                    [
                        "tf_target",
                        "prob",
                        "tf_expression",
                        "target_expression",
                        "weighted_binding",
                        "activity",
                    ]
                ]
            return result[["tf_target", "prob"]]

    def __del__(self):
        if not hasattr(self, "_tmp_files"):
            return

        for fname in self._tmp_files:
            if os.path.exists(fname):
                shutil.rmtree(fname, ignore_errors=True)


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
