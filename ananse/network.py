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
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from loguru import logger

from pybedtools import BedTool
from genomepy import Genome

from ananse import mytmpdir
import ananse

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

    def get_gene_dataframe(
        self,
        peak_bed,
        up=100000,
        down=100000,
        alpha=1e4,
        promoter=2000,
        full_weight_region=5000,
    ):
        # all overlap Enh-TSS(100000-tss-100000) pair with distance
        logger.info("relating enhancers to genes")
        peaks = BedTool(peak_bed)
        b = BedTool(self.gene_bed)
        b = b.flank(l=1, r=0, s=True, g=self.gsize).slop(  # noqa: E741
            l=up, r=down, g=self.gsize, s=True  # noqa: E741
        )
        # bedtools flank  -r 0 -l 1 -i b.bed -g
        # #all gene upstream 1bp position (TSS), Chr01 12800   12801   in Chr01    4170    12800   Xetrov90000001m.g   0   -
        # bedtools slop  -r down -l up -i b.bed -g
        # # |100000--TSS--100000|

        vals = []
        # for f in b.intersect(peaks, wo=True, nonamecheck=True):
        for f in b.intersect(peaks, wo=True):
            # bedtools intersect -wo -nonamecheck -b peaks.bed -a b.bed
            chrom = f[0]
            strand = f[5]
            if strand == "+":
                tss = f.start + up
            else:
                tss = f.start + down
            gene = f[3]
            peak_start, peak_end = int(f[-3]), int(f[-2])
            vals.append([chrom, tss, gene, peak_start, peak_end])
        p = pd.DataFrame(
            vals, columns=["chrom", "tss", "gene", "peak_start", "peak_end"]
        )
        p["peak"] = [int(i) for i in (p["peak_start"] + p["peak_end"]) / 2]
        # peak with int function, let distance int
        p["dist"] = np.abs(p["tss"] - p["peak"])
        p["loc"] = (
            p["chrom"]
            + ":"
            + p["peak_start"].astype(str)
            + "-"
            + p["peak_end"].astype(str)
        )
        p = p.sort_values("dist").drop_duplicates(["loc", "gene"], keep="first")[
            ["gene", "loc", "dist"]
        ]

        # remove distance more than 100k interaction, for weight calculate
        p = p[p["dist"] < up - 1]

        weight = self.distance_weight(
            include_promoter=self.include_promoter,
            include_enhancer=self.include_enhancer,
            remove=promoter,
            keep1=full_weight_region,
        ).set_index("dist")
        p = p.join(weight, on="dist")

        p.gene = [i.upper() for i in list(p.gene)]
        p = p.set_index("loc")
        logger.info("done")
        return p

    def distance_weight(
        self,
        include_promoter=False,
        include_enhancer=True,
        alpha=1e4,
        padding=100000,
        keep1=5000,
        remove=2000,
    ):
        """
        Built weight distribution from TSS.
        """
        # alpha is half site, default setting is 1e4, which means at 1e4 position weight is 0.5
        # padding is the full range we used
        # remove is promoter removed range
        # keep1 is keep full binding score range

        u = -math.log(1.0 / 3.0) * 1e5 / alpha

        promoter_weight = int(include_promoter)
        enhancer_weight = int(include_enhancer)

        weight1 = pd.DataFrame(
            {
                "weight": [promoter_weight for z in range(0, remove + 1)],
                "dist": range(0, remove + 1),
            }
        )

        weight2 = pd.DataFrame(
            {
                "weight": [enhancer_weight for z in range(remove + 1, keep1 + 1)],
                "dist": range(remove + 1, keep1 + 1),
            }
        )

        weight3 = pd.DataFrame(
            {
                "weight": [
                    enhancer_weight
                    * 2.0
                    * math.exp(-u * math.fabs(z) / 1e5)
                    / (1.0 + math.exp(-u * math.fabs(z) / 1e5))
                    for z in range(1, padding - keep1 + 1)
                ],
                "dist": range(keep1 + 1, padding + 1),
            }
        )

        weight = pd.concat([weight1, weight2, weight3])
        return weight

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

        if combine_function not in ["mean", "max", "sum"]:
            raise ValueError(
                "Unknown combine function, valid options are: mean, max, sum"
            )

        padding = max(up, down)
        if alpha is None:
            alpha = padding / 10

        if promoter > padding:
            raise ValueError("promoter region is > region to use")

        # Get list of unique enhancers from the binding file
        # This is relatively slow for a large file. May need some optimization.
        enhancer_bed = self.unique_enhancers(binding_fname)

        # Link enhancers to genes on basis of distance to annotated TSS
        gene_df = self.get_gene_dataframe(
            enhancer_bed,
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
        )[["factor", "enhancer", "binding"]]
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
            tmp = tmp.mean()
        elif combine_function == "sum":
            tmp = tmp.sum()

        logger.info("Done grouping...")
        return tmp

    def create_expression_network(
        self, fin_expression, column="tpm", tfs=None, rank=True
    ):
        expression = pd.DataFrame(
            pd.concat(
                [pd.read_table(f, index_col=0)[[column]] for f in fin_expression],
                axis=1,
            ).mean(1),
            columns=[column],
        )
        expression[column] = np.log2(expression[column] + 1e-5)

        expression.index.rename("target", inplace=True)
        expression = expression.reset_index()
        expression = expression.rename(columns={"tpm": "target_expression"})

        if tfs is None:
            package_dir = os.path.dirname(ananse.__file__)
            tffile = os.path.join(package_dir, "db", "tfs.txt")
            tfs = pd.read_csv(tffile, header=None)[0].tolist()

        # print(expression)
        tfs = expression[expression.target.isin(tfs)]
        tfs = tfs.reset_index()
        tfs = tfs.drop(columns=["index"])
        tfs.rename(
            columns={"target": "tf", "target_expression": "tf_expression"}, inplace=True
        )

        expression["key"] = 0
        tfs["key"] = 0

        network = expression.merge(tfs, how="outer")
        network = network[["tf", "target", "tf_expression", "target_expression"]]

        for col in ["tf_expression", "target_expression"]:
            if rank:
                network[col] = rankdata(network[col])
            network[col] = minmax_scale(network[col])

        network["tf_target"] = network["tf"] + "_" + network["target"]
        network = network.set_index("tf_target").drop(columns=["tf", "target"])
        logger.info("creating expression dataframe")
        network = dd.from_pandas(network, npartitions=30)
        # print(network)
        return network

    def unique_enhancers(self, binding_fname, chrom=None):
        p = re.compile("[:-]")
        logger.info("reading enhancers")
        enhancers = pd.read_table(binding_fname, usecols=["enhancer"])["enhancer"]
        # enhancers = pd.read_table("/ceph/rimlsfnwi/data/moldevbio/heeringen/heeringen/atac_ananse/models/saved/H3K27ac.qnorm.ref.txt", usecols=[0]).iloc[:,0]

        if chrom:
            logger.info(f"filtering for {chrom}")
            enhancers = enhancers[enhancers.str.contains(f"^{chrom}:")]

        enhancers = enhancers.unique()
        logger.info("writing temporary enhancers")
        f = NamedTemporaryFile(mode="w+", dir=mytmpdir(), delete=False)
        for e in enhancers:
            print("{}\t{}\t{}".format(*re.split(p, e)), file=f)

        # print(open(f.name).readlines()[:10])

        return f.name

    def run_network(
        self,
        binding,
        tfs=None,
        fin_expression=None,
        corrfiles=None,
        outfile=None,
        up=1e5,
        down=1e5,
        alpha=None,
        promoter=2000,
        full_weight_region=5000,
    ):
        # Expression base network
        logger.info("Loading expression")
        df_expression = self.create_expression_network(
            fin_expression, tfs=tfs, rank=True
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

            # This is where the heavy lifting of all delayed computations gets done
            logger.info("Computing network")
            with ProgressBar():
                result = df_expression.join(df_binding)
                result = result.compute()
                result = result.fillna(0)
                result["binding"] = minmax_scale(
                    rankdata(result["weighted_binding"], method="min")
                )
                # Combine binding score with expression score
                result["binding"] = result[
                    ["tf_expression", "target_expression", "binding"]
                ].mean(1)
        else:
            result = df_expression
            result["binding"] = result[["tf_expression", "target_expression"]].mean(1)

        logger.info("Saving file")
        result[["binding"]].to_csv(outfile, sep="\t")
