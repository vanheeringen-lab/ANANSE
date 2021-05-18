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
import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from loguru import logger

from genomepy import Genome
import pyranges as pr

import ananse
from ananse.utils import region_gene_overlap

warnings.filterwarnings("ignore")


class Network(object):
    def __init__(
        self,
        ncore=1,
        genome="hg38",
        gene_bed=None,
        include_promoter=False,
        include_enhancer=True,
    ):
        """[infer cell type-specific gene regulatory network]

        Arguments:
            object {[type]} -- [description]

        Keyword Arguments:
            ncore {int} -- [Specifies the number of threads to use during analysis.] (default: {1})
            genome {str} -- [The genome that is used for the gene annotation and the enhancer location.] (default: {"hg38"})
            gene_bed {[type]} -- [Gene annotation for the genome specified with -g as a 12 column BED file.] (default: {None})
            include_promoter {bool} -- [Include or exclude promoter peaks (<= TSS +/- 2kb) in network inference.] (default: {False})
            include_enhancer {bool} -- [Include or exclude enhancer peaks (> TSS +/- 2kb) in network inference.] (default: {True})

        Raises:
            TypeError: [description]
        """
        self.ncore = ncore
        self.genome = genome
        g = Genome(self.genome)
        self.gsize = g.sizes_file

        # # Motif information file
        # if pfmfile is None:
        #     self.pfmfile = "../data/gimme.vertebrate.v5.1.pfm"
        # else:
        #     self.pfmfile = pfmfile

        # self.motifs2factors = self.pfmfile.replace(".pfm", ".motif2factors.txt")
        # self.factortable = self.pfmfile.replace(".pfm", ".factortable.txt")

        package_dir = os.path.dirname(ananse.__file__)

        # Gene information file
        if self.genome == "hg38":
            if gene_bed is None:
                self.gene_bed = os.path.join(package_dir, "db", "hg38.genes.bed")
            else:
                self.gene_bed = gene_bed
        elif self.genome == "hg19":
            if gene_bed is None:
                self.gene_bed = os.path.join(package_dir, "db", "hg19.genes.bed")
            else:
                self.gene_bed = gene_bed
        else:
            if gene_bed is None:
                raise TypeError("Please provide a gene bed file with -a argument.")
            else:
                self.gene_bed = gene_bed

        # self.promoter = promoter
        self.include_promoter = include_promoter

        self.include_enhancer = include_enhancer

    def unique_enhancers(self, fname, chrom=None):
        """Extract a list of unique enhancers.

        Parameters
        ----------
        fname : str
            File name of a tab-separated file that contains an 'enhancer' column.

        chrom : str, optional
            Only return enhancers on this chromosome.

        Returns
        -------
            PyRanges object with enhancers
        """
        p = re.compile("[:-]")
        logger.info("reading enhancers")

        # Read enhancers from binding file
        # This is relatively slow for a large file. May need some optimization.
        enhancers = pd.read_table(fname, usecols=["enhancer"])["enhancer"]
        enhancers = enhancers.unique()

        # Split into columns and create PyRanges object
        p = re.compile("[:-]")
        enhancers = pr.PyRanges(
            pd.DataFrame(
                [re.split(p, e) for e in enhancers],
                columns=["Chromosome", "Start", "End"],
            )
        )
        return enhancers

    def distance_weight(
        self,
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
        include_promoer : bool, optional
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
                "weight": [promoter_weight for z in range(0, promoter_region + 1)],
                "dist": range(0, promoter_region + 1),
            }
        )

        weight2 = pd.DataFrame(
            {
                "weight": [
                    enhancer_weight
                    for z in range(promoter_region + 1, full_weight_region + 1)
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

        if combine_function not in ["mean", "max", "sum"]:
            raise ValueError(
                "Unknown combine function, valid options are: mean, max, sum"
            )

        maximum_distance = max(up, down)
        if alpha is None:
            alpha = maximum_distance / 10

        if promoter > maximum_distance:
            raise ValueError(
                "promoter region is larger than the maximum distance to use"
            )

        # Get list of unique enhancers from the binding file
        enhancer_pr = self.unique_enhancers(binding_fname)

        # Link enhancers to genes on basis of distance to annotated TSS
        gene_df = self.enhancer2gene(
            enhancer_pr,
            up=up,
            down=down,
            alpha=alpha,
            promoter=promoter,
            full_weight_region=full_weight_region,
        )

        # print(gene_df)
        logger.info("Reading binding file...")
        ddf = dd.read_csv(
            binding_fname,
            sep="\t",
            usecols=["factor", "enhancer", "binding"],
        )
        if tfs is not None:
            ddf = ddf[ddf["factor"].isin(tfs)]

        # Merge binding information with gene information.
        # This may be faster than creating index on enhancer first, but need to check!
        tmp = ddf.merge(gene_df, left_on="enhancer", right_index=True)

        # Remove everything with weight 0
        tmp = tmp[tmp["weight"] > 0]

        # Modify the binding by the weight, which is based on distance to TSS
        tmp["weighted_binding"] = tmp["weight"] * tmp["binding"]

        logger.info("Grouping by tf and target gene...")
        # Using one column that combines TF and target as dask cannot handle MultiIndex
        tmp["tf_target"] = tmp["factor"] + "_" + tmp["gene"]

        tmp = tmp.groupby("tf_target")[["weighted_binding"]]

        if combine_function == "mean":
            tmp = tmp.mean()
        elif combine_function == "max":
            tmp = tmp.max()
        elif combine_function == "sum":
            tmp = tmp.sum()

        logger.info("Done grouping...")
        return tmp

    def create_expression_network(
        self, fin_expression, column="tpm", tfs=None, rank=True, bindingfile=None
    ):
        """Create a gene expression based network.

        Based on a file with gene expression levels (a TPM column), a
        dask DataFrame is generated with the combined expression levels
        of the tf and the target gene. By default, the expresison levels
        are ranked and subsequently scaled between 0 and 1.

        Parameters
        ----------
        fin_expression : str or list
            Filename of file that contains gene expression data (TPM), or a
            list of filenames. First column should contain the gene name.

        column : str, optional
            Column name that contains gene expression, 'tpm' by default.

        tfs : list, optional
            List of TF gene names. All TFs will be used by default.

        rank : bool, optional
            Rank expression levels before scaling.

        bindingfile : str, optional
            Filename with binding information.

        Returns
        -------
            Dask DataFrame with gene expression based values.
        """
        # Convert it to a list if it's not a list of files, but a single file name
        if isinstance(fin_expression, str):
            fin_expression = [fin_expression]

        # Read all expression input files and take the mean expression per gene
        expression = pd.DataFrame(
            pd.concat(
                [pd.read_table(f, index_col=0)[[column]] for f in fin_expression],
                axis=1,
            ).mean(1),
            columns=[column],
        )
        expression[column] = np.log2(expression[column] + 1e-5)

        # Create the target gene list, based on all genes
        expression.index.rename("target", inplace=True)
        expression = expression.reset_index()
        expression = expression.rename(columns={"tpm": "target_expression"})

        # Create the TF list, based on valid transcription factors
        if tfs is None:
            activity_fname = bindingfile.replace("binding.tsv", "factor_activity.tsv")
            if os.path.exists(activity_fname):
                tfs = list(set(pd.read_table(activity_fname, index_col=0).index.tolist()))
            else:
                package_dir = os.path.dirname(ananse.__file__)
                tffile = os.path.join(package_dir, "db", "tfs.txt")
                tfs = pd.read_csv(tffile, header=None)[0].tolist()

        tfs = expression[expression.target.isin(tfs)]
        tfs = tfs.reset_index()
        tfs = tfs.drop(columns=["index"])
        tfs.rename(
            columns={"target": "tf", "target_expression": "tf_expression"}, inplace=True
        )

        expression["key"] = 0
        tfs["key"] = 0

        # Merge TF and target gene expression information
        network = expression.merge(tfs, how="outer")
        network = network[["tf", "target", "tf_expression", "target_expression"]]

        # Rank and scale
        for col in ["tf_expression", "target_expression"]:
            if rank:
                network[col] = rankdata(network[col])
            network[col] = minmax_scale(network[col])

        # Use one-column index that contains TF and target genes.
        # This is necessary for dask, as dask cannot merge on a MultiIndex.
        # Otherwise this would be an inefficient and unnecessary step.
        network["tf_target"] = network["tf"] + "_" + network["target"]
        network = network.set_index("tf_target").drop(columns=["target"])

        # Convert to a dask DataFrame.
        logger.info("creating expression dataframe")
        network = dd.from_pandas(network, npartitions=30)

        return network

    def run_network(
        self,
        binding,
        fin_expression=None,
        tfs=None,
        corrfiles=None,
        outfile=None,
        up=1e5,
        down=1e5,
        alpha=None,
        promoter=2000,
        full_weight_region=5000,
    ):
        """Create network.

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
        corrfiles : [type], optional
            Correlation files by default None. CURRENTLY UNUSED.
        outfile : str, optional
            Output file.
        up : int, optional
            Upstream maximum distance, by default 100kb.
        down : int, optional
            Upstream maximum distabce, by default 100kb.
        alpha : float, optional
            Parameter to control weight decay, by default 1e4.
        promoter : int, optional
            Promoter region, by default 2000.
        full_weight_region : int, optional
            Region that will receive full weight, by default 5000."""
        # Expression base network
        logger.info("Loading expression")
        df_expression = self.create_expression_network(
            fin_expression, tfs=tfs, rank=True, bindingfile=binding
        )

        # Use a version of the binding network, either promoter-based, enhancer-based
        # or both.
        if self.include_promoter or self.include_enhancer:
            logger.info("Aggregate binding")
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

            activity_fname = binding.replace("binding.tsv", "factor_activity.tsv")
            if os.path.exists(activity_fname):
                logger.info("Reading factor activity")
                act = pd.read_table(activity_fname, index_col=0)
                act.index.name = "tf"
                act["activity"] = minmax_scale(rankdata(act["activity"], method="min"))
                df_expression = df_expression.merge(
                    act, right_index=True, left_on="tf", how="left"
                ).fillna(0.5)
            df_expression = df_expression.drop(columns=["tf"])

            # This is where the heavy lifting of all delayed computations gets done
            logger.info("Computing network")
            if fin_expression is not None:
                with ProgressBar():
                    result = df_expression.join(df_binding)
                    result = result.compute()
                result = result.fillna(0)
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

        logger.info("Saving file")
        result[["prob"]].to_csv(outfile, sep="\t")
