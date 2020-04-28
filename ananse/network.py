#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

"""Built gene regulatory network"""

# Python imports
import os
import math
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
    # def __init__(self, ncore=1, genome="hg38", gene_bed=None, pfmfile=None, promoter=False):
    def __init__(self, ncore=1, genome="hg38", gene_bed=None):

        self.ncore = ncore
        self.genome = genome
        g = Genome(self.genome)
        self.gsize = g.props["sizes"]["sizes"]

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
                self.gene_bed =  os.path.join(package_dir, "db", "hg38_genes.bed")
            else:
                self.gene_bed = gene_bed
        elif self.genome == "hg19":
            if gene_bed is None:
                self.gene_bed =  os.path.join(package_dir, "db", "hg19_genes.bed")
            else:
                self.gene_bed = gene_bed            
        else:
            if gene_bed is None:
                raise TypeError("Please provide a gene bed file with -a argument.")
            else:
                self.gene_bed = gene_bed
        
        # self.promoter = promoter


    def clear_peak(self, ddf):
        """
        Filter the enhancer peaks in promoter range.
        """
        ddf = ddf.compute()

        global alltfs
        alltfs = list(set(ddf.factor))

        enhancerbed = pd.DataFrame(set(ddf.enhancer))
        enhancerbed[["chr","site"]]=enhancerbed[0].str.split(":",expand=True)
        enhancerbed[["start","end"]]=enhancerbed.site.str.split("-",expand=True)
        enhancerbed.drop(columns=[0,"site"], inplace=True)

        enhancerfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        enhancerbed.to_csv(enhancerfile, sep="\t", header=False, index=False)
        # print(enhancerfile.name)
        return enhancerfile.name


    def get_promoter_dataframe(self, peak_bed, up=2000, down=2000):
        # all overlap Enh-TSS(up8000 to down2000) pair
        peaks = BedTool(peak_bed)
        b = BedTool(self.gene_bed)
        b = b.flank(l=1, r=0, s=True, g=self.gsize).slop(  # noqa: E741
            l=up, r=down, g=self.gsize, s=True  # noqa: E741
        )
        vals = []
        # for f in b.intersect(peaks, wo=True, nonamecheck=True):
        for f in b.intersect(peaks, wo=True):
            chrom = f[0]
            gene = f[3]
            peak_start, peak_end = int(f[13]), int(f[14])
            vals.append([chrom, gene, peak_start, peak_end])
        prom = pd.DataFrame(vals, columns=["chrom", "gene", "peak_start", "peak_end"])
        prom["loc"] = (
            prom["chrom"]
            + ":"
            + prom["peak_start"].astype(str)
            + "-"
            + prom["peak_end"].astype(str)
        )
        prom.gene = [i.upper() for i in list(prom.gene)]
        return prom

    def get_gene_dataframe(self, peak_bed, up=100000, down=100000):
        # all overlap Enh-TSS(100000-tss-100000) pair distance
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
            peak_start, peak_end = int(f[13]), int(f[14])
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

        p = p[p["dist"] < up - 1]
        # remove distance more than 100k interaction, for weight calculate
        p.gene = [i.upper() for i in list(p.gene)]

        return p

    def distance_weight(self, alpha=1e5, padding=100000, keep1=5000, remove=2000):
        """
        Built weight distribution from TSS.
        """
        # alpha is half site, default setting is 1e5, which means at 1e5 position weight is 0.5
        # padding is the full range we used
        # remove is promoter removed range
        # keep1 is keep full binding score range

        u = -math.log(1.0 / 3.0) * 1e5 / alpha
        weight1 = pd.DataFrame(
            {"weight": [0 for z in range(1, remove + 1)], "dist": range(1, remove + 1)}
        )
        weight2 = pd.DataFrame(
            {
                "weight": [1 for z in range(remove + 1, keep1 + 1)],
                "dist": range(remove + 1, keep1 + 1),
            }
        )
        weight3 = pd.DataFrame(
            {
                "weight": [
                    2.0
                    * math.exp(-u * math.fabs(z) / 1e5)
                    / (1.0 + math.exp(-u * math.fabs(z) / 1e5))
                    for z in range(1, padding - keep1 + 1)
                ],
                "dist": range(keep1 + 1, padding + 1),
            }
        )
        weight = pd.concat([weight1, weight2, weight3])

        weightfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        weight.to_csv(weightfile)
        return weightfile.name

    def aggregate_binding(self, ddf, prom, p, weight):

        # ddf = dd.read_hdf(binding, key="/binding")[["factor", "enhancer", "binding"]]
        # ddf = dd.read_csv(binding, sep="\t")[["factor", "enhancer", "binding"]]
        prom_table = ddf.merge(prom, left_on="enhancer", right_on="loc")
        prom_table = prom_table.groupby(["factor", "gene"])[["binding"]].max()
        prom_table = prom_table.rename(columns={"binding": "max_binding_in_promoter"})
        prom_table = prom_table.reset_index()
        prom_table["source_target"] = (
            prom_table["factor"].map(str) + "_" + prom_table["gene"].map(str)
        )

        f_table = ddf.merge(p, left_on="enhancer", right_on="loc")
        sum_enh = f_table.groupby(["factor", "gene"])[["binding"]].count()
        f_table["sum_weighted_logodds"] = (
            f_table["binding"]
            .div(f_table["binding"].mean())
            .apply(np.log, meta=("binding", np.float64))
            .rmul(50000)
            .div(f_table["dist"])
        )
        f_table["sum_logodds"] = (
            f_table["binding"]
            .div(f_table["binding"].mean())
            .apply(np.log, meta=("binding", np.float64))
        )

        weight = dd.read_csv(weight)
        f_table = f_table.merge(weight, how="left", on="dist")
        f_table["sum_dist_weight"] = f_table["binding"] * f_table["weight"]

        f_table_sum = f_table.groupby(["factor", "gene"]).sum()[
            ["sum_weighted_logodds", "sum_logodds", "binding", "sum_dist_weight"]
        ]
        f_table_max = f_table.groupby(["factor", "gene"])[
            ["binding", "sum_dist_weight"]
        ].max()

        f_table_sum = f_table_sum.rename(columns={"binding": "sum_binding"})
        f_table_max = f_table_max.rename(columns={"binding": "max_binding"})
        f_table_max = f_table_max.rename(
            columns={"sum_dist_weight": "max_sum_dist_weight"}
        )

        sum_enh = sum_enh.rename(columns={"binding": "enhancers"})
        f_table_sum = f_table_sum.reset_index()
        f_table_max = f_table_max.reset_index()

        f_table = f_table.reset_index()
        sum_enh = sum_enh.reset_index()

        f_table_sum["source_target"] = f_table_sum["factor"] + "_" + f_table_sum["gene"]
        f_table_max["source_target"] = f_table_max["factor"] + "_" + f_table_max["gene"]

        f_table["source_target"] = f_table["factor"] + "_" + f_table["gene"]
        sum_enh["source_target"] = sum_enh["factor"] + "_" + sum_enh["gene"]
        f_table_max = f_table_max.rename(columns={"factor": "factor2"})
        f_table_max = f_table_max.rename(columns={"gene": "gene2"})

        f_table = f_table_sum.merge(
            f_table_max, left_on="source_target", right_on="source_target", how="outer"
        )
        f_table = f_table.merge(
            sum_enh, left_on="source_target", right_on="source_target", how="outer"
        )
        f_table = f_table.merge(
            prom_table, left_on="source_target", right_on="source_target", how="outer"
        )
        f_table = f_table[
            [
                "source_target",
                "factor",
                "gene",
                "sum_weighted_logodds",
                "sum_dist_weight",
                "sum_logodds",
                "sum_binding",
                "enhancers",
                "max_binding_in_promoter",
                "max_binding",
                "max_sum_dist_weight",
            ]
        ]
        f_table["log_sum_binding"] = (
            f_table["sum_binding"]
            .add(1e-5)
            .apply(np.log, meta=("sum_binding", np.float64))
        )
        f_table["log_enhancers"] = (
            f_table["enhancers"].add(1).apply(np.log, meta=("enhancers", np.float64))
        )
        f_table["factor"] = f_table["source_target"].str.replace("_.*", "")
        f_table["gene"] = f_table["source_target"].str.replace(".*_", "")
        f_table["max_binding_in_promoter"] = f_table["max_binding_in_promoter"].fillna(
            0
        )

        features_file = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        # print("computing, output file {}".format(features_file.name))
        with ProgressBar():
            f_table.compute(num_workers = self.ncore)
        f_table.to_hdf(features_file.name, key="/features")
        return features_file.name

        # return f_table

    def get_expression(self, fin_expression, features, min_tpm=1e-10, column="tpm"):
        # df = dd.read_hdf(features)
        # print(features.head())
        # features = dd.from_pandas(features, chunksize=100000)
        df = pd.read_hdf(features, key="/features", 
                    columns=["source_target", "factor", "gene"])
        # df = features
        # df = df[["source_target", "factor", "gene"]]
        df.source_target = [i.upper() for i in list(df.source_target)]
        df.gene = [i.upper() for i in list(df.gene)]
        df = df.set_index("source_target")
        # fa2name={}
        # fa2=open("/home/qxu/projects/regulatoryNetwork/run20180716/scripts/data/gene2name.txt","r")

        # #switch the Factor name to gene name
        # for i in fa2:
        #     a=i.split()
        #     if a[0].startswith("gene"):
        #         fa2name[a[1]]=a[0]

        # flist=[]
        # for f in list(df["factor"]):
        #     if str.lower(f) in fa2name:
        #         flist.append(fa2name[str.lower(f)])
        #     elif f in fa2name:
        #         flist.append(fa2name[f])
        #     else:
        #         flist.append("")
        # df["gfactor"]=flist

        # Take mean of all TPMs
        expression = pd.DataFrame(
            pd.concat(
                [pd.read_table(f, index_col=0)[[column]] for f in fin_expression],
                axis=1,
            ).mean(1),
            columns=[column],
        )
        expression.index = [i.upper() for i in list(expression.index)]
        # print(expression)
        expression[column] = np.log2(expression[column] + 1e-5)
        df = df.join(expression, on="factor")
        df = df.rename(columns={column: "factor_expression"})
        df = df.join(expression, on="gene")
        df = df.rename(columns={column: "target_expression"})

        df = df.dropna()

        for col in ["factor_expression", "target_expression"]:
            df[col + ".scale"] = minmax_scale(df[col])
            df[col + ".rank.scale"] = minmax_scale(rankdata(df[col]))

        # with ProgressBar():
        #     df.compute(num_workers = self.ncore)
        expression_file = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        # outfile = os.path.join(outdir, "expression.txt")
        df.to_csv(expression_file, sep="\t")
        return expression_file.name

    def get_factorExpression(self, fin_expression):
        import numpy as np
        import pandas as pd
        from scipy.stats import rankdata
        from sklearn import preprocessing
        import warnings

        warnings.filterwarnings("ignore")
        factorsExpression = {}

        # for line in open(self.motifs2factors):
        #     if not line.split("\t")[1].strip().split(",") == [""]:
        #         for factor in line.split("\t")[1].strip().split(","):
        #             factorsExpression[factor.upper()] = []
        
        for tf in alltfs:
            factorsExpression[tf] = []

        for f in fin_expression:
            with open(f) as fa:
                for line in fa:
                    if not line.startswith("target_id"):
                        gene = line.split("\t")[0].upper()
                        expression = float(line.split("\t")[1])
                        if gene in factorsExpression:
                            if expression < 1e-10:
                                expression = 1e-10
                            factorsExpression[gene].append(np.log10(expression))

        # for line in open(fin_b):
        #     if not line.startswith('target_id'):
        #         gene = line.split('\t')[0].upper()
        #         expression = float(line.split('\t')[4])
        #         if gene in factorsExpression:
        #             if expression < 1e-10:
        #                 expression = 1e-10
        #             factorsExpression[gene].append(np.log10(expression))

        factors_expression_file = NamedTemporaryFile(
            mode="w", dir=mytmpdir(), delete=False
        )
        factors_expression_file.write("#factor\tfactorExpression\n")
        for factor in factorsExpression:
            if len(factorsExpression[factor]) == 0:
                factors_expression_file.write(
                    "{}\t{}\n".format(factor, np.log10(1e-10))
                )
            else:
                factors_expression_file.write(
                    "{}\t{}\n".format(factor, np.mean(factorsExpression[factor]))
                )

        scores_df = pd.read_table(factors_expression_file.name, sep="\t", index_col=0)
        # scores_df['factorExpressionRank'] = preprocessing.MinMaxScaler().fit_transform(rankdata(scores_df['factorExpression'], method='average'))
        scores_df["factorExpressionRank"] = preprocessing.MinMaxScaler().fit_transform(
            rankdata(scores_df["factorExpression"], method="average").reshape(-1, 1)
        )

        scores_df.to_csv(factors_expression_file.name, sep="\t")
        return factors_expression_file.name

    def get_correlation(self, corrfiles, features):
        df = pd.read_hdf(features)
        df = df[["source_target"]]
        df.source_target = [i.upper() for i in list(df.source_target)]
        df = df.set_index("source_target")

        for i, corrfile in enumerate(corrfiles):
            corr = pd.read_table(corrfile, sep="\t", index_col=0)
            corr = corr.rename(columns={corr.columns[0]: "corr_file{}".format(i + 1)})
            df = df.join(corr)

        corr_file = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        # outfile = os.path.join(outdir, "correlation.txt")
        df.to_csv(corr_file, sep="\t")
        return corr_file.name

    def join_features(self, features, other):
        network = pd.read_hdf(features, key="/features")

        for fname in other:
            df = pd.read_table(fname, sep="\t")
            for col in ["factor", "gene"]:
                if col in df.columns:
                    df = df.drop(col, 1)
            network = network.merge(
                df, left_on="source_target", right_on="source_target",
            )

        # Compute before saving, will result in an error otherwise
        # network = network.compute()
        # with ProgressBar():
            # network.compute(num_workers = self.ncore)
        # featurefile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        # print(network.head())
        # network.to_csv(featurefile, sep="\t", index=False)
        # network.to_hdf(featurefile, key="/features")
        # return featurefile.name
        return network

    def create_network(self, featurefile, outfile, impute=False):

        # network = dd.read_hdf(featurefile, key="/features")
        # network = dd.read_csv(featurefile, sep="\t")
        network = featurefile
        exclude_cols = [
            "sum_weighted_logodds",
            "enhancers",
            "log_enhancers",
            "sum_binding",
            "sum_logodds",
            "log_sum_binding",
            "factorExpression",
            "targetExpression",
            "factor",
            "gene",
            "factor_expression",
            "target_expression",
            "factor_expression.scale",
            "target_expression.scale",
            # "factor_expression.rank.scale", 
            "target_expression.rank.scale",
            "corr_file1",
            "correlation",
            "correlationRank",
            "max_binding_in_promoter",
            "max_binding",
            "max_sum_dist_weight",
            #                "sum_dist_weight"
        ]
        network = network[[c for c in network.columns if c not in exclude_cols]]
        network = network.set_index("source_target")
        network["binding"] = minmax_scale(
            rankdata(network["sum_dist_weight"], method="dense")
        )
        network.drop(["sum_dist_weight"], axis=1, inplace=True)

        bp = network.mean(axis=1)
        bpd = pd.DataFrame(bp)
        bpd = bpd.rename(columns={0: "binding"})
        bpd["prob"] = minmax_scale(rankdata(bpd["binding"], method="dense"))

        # with ProgressBar():
        #     bpd.compute(num_workers = self.ncore)
        bpd.to_csv(outfile, sep="\t")

    def create_promoter_network(self, featurefile, outfile, impute=False):

        # network = pd.read_hdf(featurefile, key="/features")
        network = pd.read_csv(featurefile, sep="\t")

        exclude_cols = [
            "sum_weighted_logodds",
            "enhancers",
            "log_enhancers",
            "sum_binding",
            "sum_logodds",
            "log_sum_binding",
            "factorExpression",
            "targetExpression",
            "factor",
            "gene",
            "factor_expression",
            "target_expression",
            "factor_expression.scale",
            "target_expression.scale",
            "factor_expression.rank.scale", "target_expression.rank.scale",
            "corr_file1",
            "correlation",
            "correlationRank",
            # "max_binding_in_promoter",
            "max_binding",
            "max_sum_dist_weight",
            "sum_dist_weight"
        ]
        network = network[[c for c in network.columns if c not in exclude_cols]]
        network = network.set_index("source_target")
        network["binding"] = minmax_scale(
            rankdata(network["max_binding_in_promoter"], method="dense")
        )
        network.drop(["max_binding_in_promoter"], axis=1, inplace=True)

        bp = network.mean(axis=1)
        bpd = pd.DataFrame(bp)
        bpd = bpd.rename(columns={0: "binding"})
        bpd["prob"] = minmax_scale(rankdata(bpd["binding"], method="dense"))

        bpd.to_csv(outfile, sep="\t")

    def create_expression_network(self, featurefile, outfile, impute=False):

        # network = pd.read_hdf(featurefile, key="/features")
        network = pd.read_csv(featurefile, sep="\t")

        exclude_cols = [
            "sum_weighted_logodds",
            "enhancers",
            "log_enhancers",
            "sum_binding",
            "sum_logodds",
            "log_sum_binding",
            "factorExpression",
            "targetExpression",
            "factor",
            "gene",
            "factor_expression",
            "target_expression",
            "factor_expression.scale",
            "target_expression.scale",
            # "factor_expression.rank.scale", "target_expression.rank.scale",
            "corr_file1",
            "correlation",
            "correlationRank",
            "max_binding_in_promoter",
            "max_binding",
            "max_sum_dist_weight",
            "sum_dist_weight"
        ]
        network = network[[c for c in network.columns if c not in exclude_cols]]
        network = network.set_index("source_target")
        # network["binding"] = minmax_scale(
        #     rankdata(network["sum_dist_weight"], method="dense")
        # )
        # network.drop(["sum_dist_weight"], axis=1, inplace=True)

        bp = network.mean(axis=1)
        bpd = pd.DataFrame(bp)
        bpd = bpd.rename(columns={0: "binding"})
        bpd["prob"] = minmax_scale(rankdata(bpd["binding"], method="dense"))

        bpd.to_csv(outfile, sep="\t")

    def create_promoter_expression_network(self, featurefile, outfile, impute=False):

        # network = pd.read_hdf(featurefile, key="/features")
        network = pd.read_csv(featurefile, sep="\t")

        exclude_cols = [
            "sum_weighted_logodds",
            "enhancers",
            "log_enhancers",
            "sum_binding",
            "sum_logodds",
            "log_sum_binding",
            "factorExpression",
            "targetExpression",
            "factor",
            "gene",
            "factor_expression",
            "target_expression",
            "factor_expression.scale",
            "target_expression.scale",
            # "factor_expression.rank.scale", "target_expression.rank.scale",
            "corr_file1",
            "correlation",
            "correlationRank",
            # "max_binding_in_promoter",
            "max_binding",
            "max_sum_dist_weight",
            "sum_dist_weight"
        ]
        network = network[[c for c in network.columns if c not in exclude_cols]]
        network = network.set_index("source_target")
        network["binding"] = minmax_scale(
            rankdata(network["max_binding_in_promoter"], method="dense")
        )
        network.drop(["max_binding_in_promoter"], axis=1, inplace=True)

        bp = network.mean(axis=1)
        bpd = pd.DataFrame(bp)
        bpd = bpd.rename(columns={0: "binding"})
        bpd["prob"] = minmax_scale(rankdata(bpd["binding"], method="dense"))

        bpd.to_csv(outfile, sep="\t")

    def run_network(self, binding, fin_expression, corrfiles, outfile):
        # gene_bed="/home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed"
        # peak_bed="data/krt_enhancer.bed"
        # pfmfile="../data/gimme.vertebrate.v5.1.pfm"
        # binding="results/binding.predicted.h5"

        # b=self.interaction.Interaction(genome=self.genome, gene_bed= self.gene_bed, pfmfile=self.pfmfile)

        # print("1, Read data")
        logger.info("Read data")

        ddf = dd.read_csv(binding, sep="\t")[["factor", "enhancer", "binding"]]

        filter_bed = self.clear_peak(ddf)

        prom = self.get_promoter_dataframe(filter_bed)
        p = self.get_gene_dataframe(filter_bed)
        weight = self.distance_weight()


        logger.info("Aggregate binding")

        features = self.aggregate_binding(ddf, prom, p, weight)

        logger.info("Join expression")

        expression_file = self.get_expression(fin_expression, features)
        # factors_expression_file = self.get_factorExpression(fin_expression)

        if corrfiles is None:
            other = [
                expression_file,
            ]
        else:
            corr_file = self.get_correlation(corrfiles, features)
            other = [
                expression_file,
                corr_file,
            ]

        # outfile = 'full_features.h5'
        logger.info("Join features")
        featurefile = self.join_features(features, other)

        logger.info("Create network")
        self.create_network(featurefile, outfile)

        # if self.promoter:
        #     self.create_promoter_network(featurefile, ".".join(outfile.split(".")[:-1])+"_promoter."+outfile.split(".")[-1])
        #     self.create_expression_network(featurefile, ".".join(outfile.split(".")[:-1])+"_expression."+outfile.split(".")[-1])
        #     self.create_promoter_expression_network(featurefile, ".".join(outfile.split(".")[:-1])+"_promoter_expression."+outfile.split(".")[-1])
