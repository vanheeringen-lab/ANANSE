import itertools
import os
import re
import sys
from glob import glob
from tempfile import NamedTemporaryFile
from urllib.request import urlretrieve

import joblib  # noqa
import networkx as nx
import numpy as np
import pandas as pd
import qnorm  # noqa
from fluff.fluffio import load_heatmap_data  # noqa
from genomepy import Genome
from gimmemotifs.moap import moap
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import scan_regionfile_to_table
from gimmemotifs.preprocessing import coverage_table
from loguru import logger
from pandas import HDFStore
from pyfaidx import FastaIndexingError  # noqa
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale, scale
from tqdm.auto import tqdm

from ananse.bed import map_counts
from ananse.utils import (
    load_tfs,
    load_regions,
    get_motif_factors,
    check_cores,
    mytmpdir,
)
from . import PACKAGE_DIR

BLACKLIST_TFS = [
    "NO ORTHOLOGS FOUND",  # gimme motif2factors artifact
]


class PeakPredictor:
    atac_data = None
    histone_data = None
    p300_data = None
    cage_data = None
    _avg = None
    _dist = None
    all_data = None
    all_data_columns = []

    factor_model_db = "custom"
    factor_models = {}

    def __init__(
        self,
        reference=None,
        genome=None,
        atac_bams=None,
        histone_bams=None,
        p300_bams=None,
        cage_tpms=None,
        columns=None,
        regions=None,
        factors=None,
        pfmfile=None,
        pfmscorefile=None,
        ncore=4,
    ):
        # the reference data directory contains TF binding models,
        # and (optionally) a distribution to quantile normalize the enhancer data to.
        if reference is None:
            reference = os.path.join(PACKAGE_DIR, "db", "default_reference")
            self.factor_model_db = "default"
        if not os.path.exists(reference):
            raise NotADirectoryError(f"Could not find {reference}")
        self.data_dir = reference

        if all(b is None for b in [atac_bams, histone_bams, p300_bams, cage_tpms]):
            raise ValueError(
                "Need either "
                "- ATAC-seq, p300 or H3K27ac ChIP-seq BAM file(s), "
                "- ATAC-seq or H3K27ac ChIP-seq coverage table(s), "
                "- a CAGE bidirectional sites TPM table, "
                "- a combination of ATAC and H3K27ac data. "
                "See the documentation for examples."
            )

        if genome is None:
            logger.warning("Assuming genome is hg38")
            genome = "hg38"
        self.genome = genome
        self.pfmfile = pfmfile
        self.ncore = ncore

        # load factor2motifs (f2m) and motif_graph
        self._load_factor2motifs(factors=factors)

        # load CAGE bidirectional sites data and
        # load regions, region_factor_df (optional: _avg)
        # load models
        if cage_tpms is not None:
            self._load_cage(cage_tpms, regions, pfmscorefile)
            self._load_models()
            return  # CAGE-seq does not combine with other data types

        # load regions, region_factor_df (optional: _avg, _dist)
        self._load_regions(regions, pfmscorefile)

        # load ATAC-seq, p300 and/or H3K27ac ChIP-seq data
        if atac_bams is not None:
            self.atac_data = self._load_enhancer_data(
                "ATAC", atac_bams, columns, window=200
            )

        if histone_bams is not None:
            self.histone_data = self._load_enhancer_data(
                "H3K27ac", histone_bams, columns, "H3K27ac", window=2000
            )
        elif p300_bams is not None:
            self.p300_data = self._load_enhancer_data(
                "p300", p300_bams, columns, window=500
            )

        # load models
        self._load_models()

    def _load_factor2motifs(self, indirect=True, factors=None):
        """Load motif-associated data.

        For now, only default motifs are supported.
        Will read factors associated to motifs, and generates a graph of
        related factors based on different factors binding to the same motif.
        This information is used to select the most appropriate TF model.

        Sets
        ----
        self.f2m : dict
        self.motif_graph : nx.Graph

        Parameters
        ----------
        indirect : bool, optional
            Include TF-motif associations that are not curated, for instance
            based on ChIP-seq motif prediction, or binding inference. This will
            greatly increase TF coverage. By default True.
        """
        if self.pfmfile is None:
            logger.info("Loading default motif file")
        else:
            logger.info(f"Loading specified motif file: {self.pfmfile}")

        # load factor2motifs
        self.f2m = read_factor2motifs(
            pfmfile=self.pfmfile, indirect=indirect, factors=factors
        )

        # The default motif db must be pruned
        if self.pfmfile is None:
            self.f2m = _prune_f2m(self.f2m, self.genome)

        n = len(self.f2m)
        logger.info(f"  Using motifs for {n} factor{'' if n == 1 else 's'}")

        # load motif_graph
        self._jaccard_motif_graph()

    def _jaccard_motif_graph(self):
        """
        Create a graph of TF nodes where edges are the Jaccard index of the motifs that they bind to.
        For instance, if TF1 binds motif A and B and TF2 binds motif B and C,
        then the edge weight of TF1 to TF2 will be 1/3.

        Sets
        ----
        self.motif_graph : nx.Graph
        """
        # load the complete (unfiltered) factor2motifs
        complete_f2m = read_factor2motifs(self.pfmfile)
        if self.pfmfile is None:
            # The default motif db must be pruned
            complete_f2m = _prune_f2m(complete_f2m, self.genome)

        # convert motif lists to sets
        for k, v in complete_f2m.items():
            complete_f2m[k] = set(v)

        # if an alternative motif2factors is used, we can use the
        # jaccard index to link the user's TFs to
        # ortholog TF models in the ananse reference database
        if self.pfmfile is not None:
            reference_orthologs = _prune_f2m(read_factor2motifs(), "hg38")
            for k, v in reference_orthologs.items():
                if k in complete_f2m:
                    complete_f2m[k].update(set(v))
                else:
                    complete_f2m[k] = set(v)

        # compute the jaccard index between each combination of TFs
        self.motif_graph = nx.Graph()
        combinations = set(itertools.combinations(complete_f2m.keys(), 2))
        for tf1, tf2 in combinations:
            i = len(complete_f2m[tf1].intersection(complete_f2m[tf2]))
            u = len(complete_f2m[tf1].union(complete_f2m[tf2]))
            if i and u:  # both > 0
                jaccard = i / u
                self.motif_graph.add_edge(tf1, tf2, weight=1 - jaccard)

    def _load_cage(self, cage_tpms, regions=None, pfmscorefile=None, window=200):
        """Load CAGE data: Bidirectional regions and TPMs

        The CAGE file needs to contain two columns:
            1. regions in chr:start-end format
            2. TPMs

        Sets
        ----
        self.regions : list
            Bidirectional regions are loaded.
        self.region_factor_df : pd.DataFrame
            Bidirectional regions are scanned for motifs.
        self.cage_data : pd.DataFrame
            CAGE TPMs are log transformed.
        self._avg : pd.DataFrame
            Bidirectional regions are used to determine ReMap 2022 ChIP-seq coverage
            (Only hg19 or hg38 supported).

        Parameters
        ----------
        cage_tpms : str
            CAGE-seq output with bidirectional regions and TPMs
        regions : list, optional
            list of regions to filter the CAGE bidirectional regions/pfmscorefile by
        pfmscorefile : str, optional
            pre-scanned gimmemotifs scores file of the CAGE bidirectional regions
        window : int, optional
            scanning window around the center of each region to scan for motifs
        """
        logger.info("Loading CAGE data")
        data = _load_cage_tpm(cage_tpms, window)

        # Get the overlap of normalized regions between
        # 1) CAGE data, 2) regions and 3) pfmscorefile.
        # Then get their motif scores
        if regions is None:
            regions = data.index.tolist()
        else:
            regions = _remove_invalid_regions(regions, data.index)

        if pfmscorefile is not None:
            # if we have a pre-scanned file,
            # use regions from the index (optionally filtered by regions)
            self.region_factor_df = self._load_prescanned_motifs(pfmscorefile, regions)
            self.regions = self.region_factor_df.index.tolist()
        else:
            # Scan for motifs
            self.regions = regions
            self.region_factor_df = self._scan_motifs()

        logger.info(f"  Using {len(self.regions)} regions")
        if len(self.regions) != len(data):
            data = data.loc[self.regions]

        # normalize & save TPMs
        logger.debug("Log transforming CAGE data")
        data = np.log1p(data)
        # data = qnorm.quantile_normalize(data)
        data.loc[:, :] = minmax_scale(data)
        self.cage_data = data
        self.all_data_columns.append("CAGE")

        # Determine average ReMap coverage (Currently only for genomes "hg19" or "hg38")
        if "hg19" in self.genome or "hg38" in self.genome:
            self._avg = _load_cage_remap_data(data, window, self.genome, self.ncore)
            self.all_data_columns.append("average")

    def _load_regions(self, regions=None, pfmscorefile=None):
        """
        sets
        ----
        self.regions : list
            a list of regions to work with
        self.region_factor_df : pd.DataFrame
            with regions as index, motifs as columns, and motif scores as values
        """
        if regions is not None:
            # unique regions
            regions = list(dict.fromkeys(regions))

        if pfmscorefile is not None:
            # if we have a pre-scanned file,
            # use regions from the index (optionally filtered by regions)
            self.region_factor_df = self._load_prescanned_motifs(pfmscorefile, regions)
            self.regions = self.region_factor_df.index.tolist()
        elif regions is not None:
            # if we have custom regions we have to scan for motifs.
            self.regions = regions
            self.region_factor_df = self._scan_motifs()
        else:
            # assume regions, region_factor_df (and maybe more) is
            # included in the custom reference data.
            self._load_custom_data()
        logger.info(f"  Using {len(self.regions)} regions")

    def _load_prescanned_motifs(self, pfmscorefile, regions=None):
        """
        Use pre-scanned gimmemotifs motif scores.

        Parameters
        ----------
        pfmscorefile : str
            pre-scanned gimmemotifs scores file
        regions : list, optional
            list of regions to filter the pfmscorefile by
        """
        logger.info("Loading pre-scanned motif scores")

        region_motif_df = pd.read_table(pfmscorefile, comment="#", index_col=0)

        if regions is not None:
            regions = _remove_invalid_regions(regions, region_motif_df.index)

            # filter pfmscorefile and regions for overlap
            overlap = [r for r in regions if r in region_motif_df.index]
            if len(overlap) < len(region_motif_df.index):
                logger.debug(
                    f"Subsetting pfmscorefile to requested {len(regions)} regions"
                )
                region_motif_df = region_motif_df.loc[overlap]

        region_factor_df = self._average_motifs_per_factor(region_motif_df)
        return region_factor_df

    def _scan_motifs(self, **kwargs):
        """Scan regions for motifs using gimmemotifs.
        Both regions and motifs have already been filtered.

        Parameters
        ----------
        kwargs : dict, optional
            arguments passed to gimmemotifs' scan_regionfile_to_table()
        """
        logger.info("Scanning regions for motifs")

        # only scan motifs for our factors
        all_motifs = read_motifs(self.pfmfile)
        our_motifs = [m for f2m_vals in self.f2m.values() for m in f2m_vals]
        motifs = [m for m in all_motifs if m.id in our_motifs]

        with NamedTemporaryFile("w") as regionfile:
            print("region", file=regionfile)
            for region in self.regions:
                print(region, file=regionfile)
            regionfile.flush()

            region_motif_df = scan_regionfile_to_table(
                regionfile.name,
                self.genome,
                scoring="score",
                pfmfile=motifs,
                ncpus=self.ncore,
                **kwargs,
            )

        region_factor_df = self._average_motifs_per_factor(region_motif_df)
        return region_factor_df

    def _load_custom_data(self):
        """
        Load custom reference data, which must contain:
        * Peak regions (self.regions)
        * Factor scores in peak regions (self.region_factor_df)
        And optionally contains:
        * The average peak coverage (self._avg)
        * The distance from the peak to nearest TSS. (self._dist)

        We built this function with our hg38-specific REMAP dataset in mind.
        Other datasets will need to mirror that dataset structure.
        """
        logger.info("Loading the custom reference data")
        fpath = os.path.join(self.data_dir, "reference.factor.feather")
        if not os.path.exists(fpath):
            raise FileNotFoundError(
                f"{fpath} not found. For hg38, download the REMAP "
                "reference dataset and specify its location with -R."
            )

        # Read factor scores
        self.region_factor_df = pd.read_feather(fpath)
        self.region_factor_df.set_index(self.region_factor_df.columns[0], inplace=True)

        # Set regions
        self.regions = self.region_factor_df.index.tolist()

        # Read average coverage
        fname = f"{self.data_dir}/reference.coverage.txt"
        if os.path.exists(fname):
            self._avg = pd.read_table(
                fname,
                comment="#",
                index_col=0,
            )
            self._avg.columns = ["average"]
            self._avg["average"] = self._avg["average"] / self._avg["average"].max()
            self._avg = self._avg.loc[self.regions]  # order
            assert self.region_factor_df.index.identical(self._avg.index)
            self.all_data_columns.append("average")

        # Read distance to TSS
        fname = f"{self.data_dir}/reference.dist_to_tss.txt"
        if os.path.exists(fname):
            self._dist = pd.read_table(
                fname,
                comment="#",
                index_col=0,
            )
            self._dist = self._dist.loc[self.regions]  # order
            assert self.region_factor_df.index.identical(self._dist.index)
            self.all_data_columns.append("dist")

    def _average_motifs_per_factor(self, region_motif_df):
        """
        Accepts a dataframe with regions as index and motifs as columns.
        Returns a dataframe with regions as index, and factors as columns.
        For each factor, average the scores of all associated motifs found in the input.
        Factors for which no motifs were found are skipped.
        """
        region_factor_df = pd.DataFrame(index=region_motif_df.index)
        motifs_found = set(region_motif_df.columns)
        # for each factor, average the score of all its motifs
        for factor, motifs in self.f2m.items():
            if len(set(motifs) & motifs_found) == 0:
                logger.warning(f"No motifs for factor '{factor}' were found")
                continue
            region_factor_df[factor] = region_motif_df[motifs].mean(1)
        return region_factor_df

    def _load_enhancer_data(
        self, dtype, infiles, columns=None, target_distribution="ATAC", window=200
    ):
        """
        Read one or more BAM files, or one counts table.
        Returns a dataframe with binned regions as index and normalized read counts as values
        """
        logger.info(f"Loading {dtype} data")

        if _istable(infiles):
            data = self._load_counts(infiles, columns, window)
        else:
            data = self._load_bams(infiles, window)

        data = self._normalize_reads(dtype, data, target_distribution)

        if self.factor_model_db == "custom":
            data = self._load_custom_enhancer_data(dtype, data, target_distribution)

        self.all_data_columns.append(dtype)
        return data

    def _load_bams(self, bams, window=200):
        """Bin BAM reads to regions within the given window

        Parameters
        ----------
        bams: list
            One or more bam files
        window : int, optional
            search window for enhancers within peaks
        """
        tmp = pd.DataFrame(index=self.regions)
        n_regions = len(self.regions)
        with NamedTemporaryFile(mode="w") as f_out:

            for region in self.regions:
                print("{}\t{}\t{}".format(*re.split("[:-]", region)), file=f_out)
            f_out.flush()

            t = tqdm(bams, total=len(bams), unit="BAM", desc="Reading")
            for bam in t:
                t.set_description(f"Reading {os.path.basename(bam)}. Overall progress")
                name, regions, r_data, _ = load_heatmap_data(
                    f_out.name,
                    bam,
                    bins=1,
                    up=window // 2,
                    down=window // 2,
                    rmdup=True,
                    rmrepeats=True,
                )
                r_data = r_data.T[0]
                if len(r_data) == n_regions:
                    tmp[name] = r_data
                else:
                    # "regions" contains overlap, "r_data.T[0]" contains scores (order is identical)
                    df = pd.DataFrame(
                        np.array(regions)[:, 0:3],
                        columns=["chrom", "start", "end"],
                        dtype=str,
                    )
                    df["region"] = df["chrom"] + ":" + df["start"] + "-" + df["end"]
                    df[name] = r_data
                    df = df[["region", name]].set_index("region")

                    # merge on regions to get scores for overlapping regions (set rest to 0)
                    tmp = tmp.join(df)
                    tmp[name].fillna(0.0, inplace=True)
        return tmp

    def _load_counts(self, table, columns=None, window=200):
        """Load counts from a seq2science raw counts table.

        Parameters
        ----------
        table : str or list
            Counts table with the number of reads under each peak per sample
        columns : list or string, optional
            Name of the columns (case-insensitive) in the counts table to use
            (default: all columns)
        window : int, optional
            search window for enhancers within peaks
        """
        # error checking
        if isinstance(table, list):
            table = table[0]
        if not os.path.exists(table):
            raise FileNotFoundError(f"Could not find {table}")

        # load & filter df
        df = pd.read_table(table, comment="#", index_col=0)
        if any(df.index.duplicated()):
            logger.info(
                f"  Averaging counts for duplicate regions in {os.path.basename(table)}"
            )
            df = df.groupby(df.index).mean()
        if len(set(self.regions) & set(df.index)) != len(self.regions):
            logger.debug("  Mapping to regions")
            df = map_counts(self.regions, df)
        elif len(self.regions) != len(df.index):
            df = df[df.index.isin(self.regions)]
        if len(df) == 0:
            raise ValueError(
                "regionsfile (or pfmscorefile) and counts table don't overlap!"
            )
        if isinstance(columns, str):
            columns = [columns]
        if isinstance(columns, list):
            logger.info(f"  Using {len(columns)} columns")
            cols = "|".join(columns)
            re_column = re.compile(rf"^{cols}$", re.IGNORECASE)
            df = df.filter(regex=re_column)  # noqa
            if len(df.columns) != len(columns):
                logger.warning(
                    f"{len(columns)} columns requested, but only {len(df.columns)} "
                    f"found in {os.path.basename(table)}"
                )
        elif columns is None:
            logger.debug(f"  Using all {len(df.columns)} columns")
        else:
            raise TypeError(
                f"columns must be a sting, list or None. Received {type(columns)}"
            )
        if len(df.columns) == 0:
            raise ValueError("Columns must contain at least one valid column name!")

        # check the average region distance
        dist = df.index.str.split(r":|-", expand=True)
        dist = dist.to_frame()
        dist["dist"] = dist[2].astype(int) - dist[1].astype(int)  # end-start
        dist = int(dist["dist"].mean(axis=0))
        if abs(1 - dist / window) > 0.1:  # allow some variation
            logger.warning(f"Expected region width is {window}, got {dist}.")
        return df

    def _normalize_reads(self, dtype, df, target_distribution):
        """
        For reproducibility, we use specific distributions to normalize the data too.
        If unavailable, log transform the values.

        Parameters
        ----------
        dtype : string
            Data type
            (options: ["ATAC", "H3K27ac", "p300])
        df : pd.DataFrame
            Dataframe with reads under peaks
        target_distribution : str
            Distribution to quantile normalize to
            (options: ["ATAC", "H3K27ac"])
        """
        fname = f"{self.data_dir}/{target_distribution}.qnorm.ref.txt.gz"
        if os.path.exists(fname):
            logger.debug(f"Quantile normalizing {dtype} data")
            qnorm_ref = pd.read_table(fname, comment="#", index_col=0)[
                "qnorm_ref"
            ].values
            if len(self.regions) != len(qnorm_ref):
                qnorm_ref = np.random.choice(
                    qnorm_ref, size=len(self.regions), replace=True
                )

            df = qnorm.quantile_normalize(df, target=qnorm_ref)
        else:
            logger.debug(f"Log transforming {dtype} data")
            df = np.log1p(df)

        # Limit memory usage by using float16
        df = df.mean(1).astype("float16").to_frame(dtype)

        # Scale the data
        df[dtype] = df[dtype] / df[dtype].max()

        return df

    def _load_custom_enhancer_data(self, dtype, df, target_distribution):
        # optional: add a column with relative activity
        fname = f"{self.data_dir}/{target_distribution}.mean.ref.txt.gz"
        if os.path.exists(fname):
            logger.debug(f"Adding relative {dtype} data")
            mean_ref = pd.read_table(fname, comment="#", index_col=0)
            mean_ref = mean_ref.loc[df.index]
            # TODO: overwrite the original column instead (benchmark)
            # TODO: test if this is OK, instead of the comment below
            df[f"{dtype}.relative"] = scale(df[dtype] - mean_ref["mean_ref"])
            self.all_data_columns.append(f"{dtype}.relative")
            # if mean_ref.shape[0] == df.shape[0]:
            #     mean_ref.index = df.index
            #     df[f"{dtype}.relative"] = scale(df[dtype] - mean_ref["mean_ref"])
            #     self.all_data_columns.append(f"{dtype}.relative")
            # else:
            #     logger.debug(f"Regions of {fname} are not the same as input regions.")
            #     logger.debug("Skipping calculation of relative values.")
        return df

    def _load_models(self):
        """Select the models to use for binding prediction.

        Different models are available depending on the amount of input data
        (e.g. ATAC, H3K27ac, ATAC and H3K27ac)

        Additional models are available in the REMAP reference for hg38.
        REMAP's reference regions will have the most information.
        """
        model_dir = self._model_input()
        logger.debug(f"Loading {model_dir} models")
        for fname in glob(os.path.join(self.data_dir, model_dir, "*.pkl")):
            factor = os.path.basename(fname).replace(".pkl", "")
            self.factor_models[factor] = joblib.load(fname)
        logger.debug(f"  Using {len(self.factor_models)} models")

    def _model_input(self):
        """
        The models accept input with a specific number and order of columns.
        The number of columns is reflected in the names of the directories,
        The order is the sorted(columns), except with p300 which uses H3K27ac models.

        default_reference model directories:
          - ATAC_motif              2 columns
          - H3K27ac_motif           2 columns
          - CAGE_average_motif      3 columns
          - ATAC_H3K27ac_motif      3 columns

        REMAP reference model directories:
          - ATAC_ATAC.relative_average_dist_motif           5 columns
          - H3K27ac_average_dist_motif                      4 columns
          - ATAC_ATAC.relative_H3K27ac_average_dist_motif   6 columns
        """
        # We use p300 as if it is H3K27ac, so the order is not exactly the same as sorted()
        column_order = [
            "ATAC",
            "ATAC.relative",
            "CAGE",
            "H3k27ac",
            "p300",
            "average",
            "dist",
            "motif",
        ]
        for col in column_order.copy():
            if col not in self.all_data_columns + ["motif"]:
                column_order.remove(col)
        self.all_data_columns = column_order

        model_dir = "_".join(column_order)
        # TODO: Branco, does this always assume the hg19/hg38-specific self._avg data?
        model_dir = model_dir.replace("CAGE_motif", "CAGE_average_motif")
        model_dir = model_dir.replace("p300", "H3K27ac")
        return model_dir

    def predict_binding_probability(self, factor, jaccard_cutoff=0.0):
        """Predict binding probability for a transcription factor.

        Result is based on the data that been loaded: ATAC-seq,
        H3K27ac ChIP-seq, p300 ChIP-seq, ATAC-seq and H3K27ac ChIP-seq
        or CAGE-seq bidirectional regions.

        Parameters
        ----------
        factor : str
            Transcription factor name.
        jaccard_cutoff : float, optional
            Cutoff for the minimum jaccard overlap between motifs of two TFs for them to be considered related.
            Related TFs can share models. Default = 0 (0.1 seems to work well based on subjective testing).

        Returns
        -------
        pandas.DataFrame
            DataFrame with binding probabilities
        """
        if factor not in self.f2m:
            raise ValueError(f"Motif not known for {factor}")

        model = self._load_factor_model(factor, jaccard_cutoff)
        data = self._load_factor_data(factor)
        probabilities = model.predict_proba(data)[:, 1]

        return pd.DataFrame(probabilities, index=self.regions)

    def _load_factor_data(self, factor):
        if self.all_data is None:
            self._aggregate_data()
        tmp = self.region_factor_df[[factor]].rename(columns={factor: "motif"})
        tmp = tmp.join(self.all_data).dropna()
        tmp = tmp[self.all_data_columns]
        return tmp

    def _aggregate_data(self):
        tmp = None
        for data in [
            self.atac_data,
            self.histone_data,
            self.p300_data,
            self.cage_data,
            self._avg,
            self._dist,
        ]:
            if data is None:
                continue
            logger.info(f"Aggregating {data.columns[0]} data")
            if tmp is None:
                tmp = data
            else:
                tmp = tmp.join(data)
        tmp = tmp.dropna()
        self.all_data = tmp

    def _load_factor_model(self, factor, jaccard_cutoff=0.0):
        """
        Load TF-binding model that is:
        1. trained for that specific TF
        2. trained on a different TF, but with a motif overlap (jaccard similarity) larger than the cutoff
        3. a general TF binding model if the other options are not available

        Parameters
        ----------
        factor : str
            Transcription factor name.
        jaccard_cutoff : float, optional
            minimum jaccard similarity score that is needed to use the model of TF1 for TF2.
            0: any shared motifs, 1: all motifs shared. Default is 0.
        """
        model = None
        # 1. trained for that specific TF
        if factor in self.factor_models:
            logger.info(f"Using {factor} model")
            model = self.factor_models[factor]

        # 2. trained on a different TF, but with a motif jaccard index larger than the cutoff
        elif factor in self.motif_graph:
            # scores are inverted due to nx's cutoff method
            tfs = nx.single_source_dijkstra_path_length(
                self.motif_graph, factor, 1 - jaccard_cutoff
            )
            # tfs are sorted best>worst. Use the first one with a trained model
            for tf in tfs:
                if tf in self.factor_models:
                    ji = round(1 - tfs[tf], 2)
                    logger.info(f"Using {tf} model for {factor} (jaccard index {ji})")
                    model = self.factor_models[tf]
                    break

        # 3. a general TF binding model if the other options are not available
        if model is None:
            logger.info(f"Using general model for {factor} (no related TF found)")
            model = self.factor_models["general"]

        return model

    def predict_factor_activity(self, nregions=50_000):
        """Predict TF activity.

        Predicted based on motif activity using ridge regression.

        Parameters
        ----------
        nregions : int, optional
            number of regions to sample from
        """
        # Run ridge regression using motif score to predict (relative) ATAC/H3K27ac/CAGE signal
        activity = pd.DataFrame()
        state = np.random.RandomState(567)  # Consistently select same regions
        for df in (self.atac_data, self.histone_data, self.p300_data, self.cage_data):
            if df is None:
                continue

            for col in df.columns:
                with NamedTemporaryFile() as f:
                    signal = df[col].astype("float32")  # float16 will give NaN's
                    signal = pd.DataFrame({col: scale(signal)}, index=df.index)
                    # Run 3 times for more stable result
                    for i in range(3):
                        logger.info(
                            f"    Motif activity prediction on {col} data, run {i+1}/3"
                        )
                        if len(df) <= nregions:
                            signal.to_csv(f.name, sep="\t")
                        else:
                            signal.sample(nregions, random_state=state).to_csv(
                                f.name, sep="\t"
                            )
                        try:
                            activity = activity.join(
                                moap(
                                    f.name,
                                    genome=self.genome,
                                    method="bayesianridge",
                                    pfmfile=self.pfmfile,
                                    ncpus=self.ncore,
                                ),
                                how="outer",
                                rsuffix=f"_{i}",
                            )
                        except Exception as e:
                            print(e)

        # Rank aggregation
        for col in activity:
            activity[col] = rankdata(activity[col])
        activity = activity.mean(1)
        activity[:] = minmax_scale(activity)

        # Take the maximum activity from the motifs of each factor
        factor_activity = []
        for factor, motifs in self.f2m.items():
            act = activity.loc[motifs].max()
            factor_activity.append([factor, act])

        factor_activity = pd.DataFrame(factor_activity, columns=["factor", "activity"])

        return factor_activity


def read_factor2motifs(pfmfile=None, indirect=True, factors=None):
    motifs = read_motifs(pfmfile, as_dict=True)
    f2m = {}
    for name, motif in motifs.items():
        for factor in get_motif_factors(motif, indirect=indirect):
            # user specified filter
            if factors is not None and factor not in factors:
                continue

            f2m.setdefault(factor, []).append(name)

    # remove blacklisted TFs
    for tf in BLACKLIST_TFS:
        if tf in f2m:
            del f2m[tf]

    if len(f2m) == 0:
        raise ValueError(
            "Zero factors remain after filtering the motif2factors.txt associated with pfmfile!"
        )
    return f2m


def _prune_f2m(f2m, genome):
    """
    The default motif db contains human and mouse factors.
    Here we filter for the correct factors, and
    give a warning if you do not (seem to) have human or mouse data.
    """
    pruned_f2m = {}
    # check if the genome is human, mouse or something else
    species = _get_species(genome)
    if species == "human":
        # 1) ignore invalid human factors and all mouse factors
        # 2) make sure all factors are upper case (should not be needed)
        # 3) extend valid human factors with valid mouse ortholog motifs
        valid = _load_human_factors()
        for k, v in f2m.items():
            k = k.upper()
            if k in valid:
                pruned_f2m.setdefault(k, []).extend(v)
    elif species == "mouse":
        # Remove human factors
        for k, v in f2m.items():
            if k[1:].islower():
                pruned_f2m[k] = v
    else:
        warnings = [
            f"The genome '{genome}' is not recognized as human or mouse.",
            "If you do have another species, the motif file likely needs to be adapted.",
            "Currently mouse and human gene names are used to link motif to TFs.",
            "If your gene symbols are different, then you will need to create a new mapping",
            "and use the `-p` argument. For a possible method to do this, see here:",
            "https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors",
        ]
        for warn in warnings:
            logger.warning(warn)
        pruned_f2m = f2m
    return pruned_f2m


def _get_species(genome):
    """Returns 'human', 'mouse' or None.
    (ANANSE has default improvements for these two species)
    """
    # Try to get taxonomy id for genomepy managed genome.
    # If there is a taxonomy id, we can be really sure about the species.
    # If genome doesn't have a tax_id, then it will be 'na'.
    try:
        tax_id = Genome(genome).tax_id
        species = None
        if tax_id == 9606:
            species = "human"
        elif tax_id == 10090:
            species = "mouse"
        if isinstance(tax_id, int):
            # tax_id converts to int, so it is valid, must be not human or mouse
            return species
    except (FileNotFoundError, FastaIndexingError):
        pass

    # backup: try to get the species from the filename
    mapping = {
        "hg38": "human",
        "hg19": "human",
        "GRCh3": "human",
        "Homo": "human",
        "mm10": "mouse",
        "mm9": "mouse",
        "GRCm3": "mouse",
        # "Mus": "mouse",
    }
    base_genome = os.path.basename(genome.strip("/"))
    for name, species in mapping.items():
        if name in base_genome:
            return species


def _load_human_factors():
    tf_xlsx = os.path.join(PACKAGE_DIR, "db", "lovering.tfs.xlsx")
    valid_factors = pd.read_excel(
        tf_xlsx,
        engine="openpyxl",
        sheet_name=1,
    )
    valid_factors = valid_factors.loc[
        valid_factors["Pseudogene"].isnull(), "HGNC approved gene symbol"
    ].values
    valid_factors = list(set(valid_factors) - {"EP300"})
    return valid_factors


def _remove_invalid_regions(regions, other):
    """
    If regions are specified, check that all are found in the enhancer data.
    Give a warning for each invalid region, but continue running with the remainder.

    If regions are unspecified, create a regions list from the enhancer data.
    """
    # warn if regions are not found in the file with enhancer data (likely a typo/later addition)
    invalid = set(regions) - set(other)
    if len(invalid) > 0:
        logger.warning(f"{len(invalid)} regions not found enhancer data:")
        logger.warning(", ".join(list(invalid)))
        logger.warning("These regions will be ignored.")
        for region in list(invalid):
            regions.remove(region)
    return regions


def _load_cage_tpm(cage_tpms, window=200):
    """load CAGE-seq TPM file, standardize regions to window size,
    and average duplicated regions."""
    data = pd.read_table(cage_tpms, comment="#", dtype=str)
    if len(data.columns) != 2:
        raise ValueError(
            "For the CAGE TPM file, please give only two columns: regions and TPMs"
        )
    if len(data) == 0:
        raise ValueError("CAGE TPM file is empty.")
    data.columns = ["regions", "CAGE"]
    data["CAGE"] = data["CAGE"].astype(np.float32)

    # Normalize regions to the window size
    regions_split = data["regions"].str.split("[:-]", expand=True)
    regions_split.columns = ["chrom", "start", "end"]
    regions_split["start"] = regions_split["start"].astype(int)
    regions_split["end"] = regions_split["end"].astype(int)
    center = (regions_split["start"] + regions_split["end"]) // 2
    regions_split["start"] = center - window // 2
    regions_split["end"] = center + window // 2

    # Merge normalized regions into chr:start-end format
    data["regions"] = (
        regions_split["chrom"]
        + ":"
        + regions_split["start"].astype(str)
        + "-"
        + regions_split["end"].astype(str)
    )
    data.set_index("regions", inplace=True)
    if any(data.index.duplicated()):
        logger.info("  Averaging TPMs for duplicate regions in CAGE file")
        data = data.groupby(data.index).mean()
    return data


def _load_cage_remap_data(data, window, genome, ncore):
    """load the REMAP CAGE-seq bigwig and determine coverage in target regions"""
    genome = "hg19" if "hg19" in genome else "hg38"
    logger.info(f"  Using REMAP coverage for {genome}!")
    coverage_bw_path = os.path.join(PACKAGE_DIR, "db", f"remap2022.{genome}.w50.bw")
    if not os.path.exists(coverage_bw_path):
        logger.info("Downloading bigwig...")
        link = f"https://zenodo.org/record/6404593/files/remap2022.{genome}.w50.bw"
        _ = urlretrieve(link, coverage_bw_path)

    # Create a regions file from CAGE input
    regions_bed = os.path.join(mytmpdir(), "regions.bed")
    data.index.to_series().str.split("[:-]", expand=True).to_csv(
        regions_bed, header=False, index=False, sep="\t"
    )

    logger.info("Determining average peak coverage")
    remap_cov = coverage_table(
        peakfile=regions_bed,
        datafiles={coverage_bw_path},
        window=window,
        ncpus=ncore,
    )
    remap_cov = remap_cov.set_index(data.index)
    # remap_cov.rename(columns={remap_cov.columns[0]: "average"}, inplace=True)
    remap_cov.columns = ["average"]
    remap_cov["average"] = remap_cov["average"] / remap_cov["average"].max()
    return remap_cov


def _istable(arg):
    as_str = isinstance(arg, str) and arg.endswith(".tsv")
    as_lst = isinstance(arg, list) and len(arg) == 1 and arg[0].endswith(".tsv")
    return as_str or as_lst


def _check_input_files(*args):
    files = []
    for arg in args:
        if arg is None:
            continue
        if isinstance(arg, list):
            files.extend(arg)
        else:
            files.append(arg)

    all_files_found = True
    for fname in files:
        if not os.path.exists(fname):
            logger.exception(f"Could not find {fname}!")
            all_files_found = False

    if not all_files_found:
        sys.exit(1)


def predict_peaks(
    outdir,
    atac_bams=None,
    histone_bams=None,
    p300_bams=None,
    cage_tpms=None,
    columns=None,
    regions=None,
    reference=None,
    factors=None,
    genome=None,
    pfmfile=None,
    pfmscorefile=None,
    jaccard_cutoff=0.0,
    ncore=4,
):
    """Predict binding in a set of genomic regions.

    Binding is predicted based on ATAC-seq and/or H3K27ac ChIP-seq data in
    combination with motif scores. The model that is used is flexible, based
    on the input data. The most accurate model will be the one that uses the
    references regions in combination with both ATAC-seq and H3K27ac ChIP-seq.

    In addition, CAGE-seq bidirectional sites (TPM) can be used as a proxy for
    cis-regulatory elements. The most accurate model uses ReMap2022 TF peak data.
    (Currently, only hg19 and hg38 have been taken along)

    The result will be saved to an output file called `binding.h5` in the
    output directory, specified by the `outdir` argument. This file wil contain
    three columns: factor, enhancer and binding. The binding columns represents
    the binding probability.

    To predict binding, `predict_peaks()` needs a set of input regions. For
    human, you have two options. You can either use the reference set of
    putative enhancer regions, as described in the ANANSE manuscript [1]. This
    is specified by the `reference` argument.
    Alternatively, you can specify one or more region files with the
    `regionfiles` argument. These are files in BED or narrowPeak format, that
    describe potential enhancers. For instance, a reference enhancer set, peaks
    from your ATAC-seq experiments or any other collection of regions. For
    accurate motif analysis, these should be as precise as possible. BroadPeaks
    from histone ChIP-seq are not really suitable. NarrowPeaks from ATAC-seq,
    DNase-seq or TF ChIP-seq will be fine.

    Parameters
    ----------
    outdir : str
        Name of output directory.
    atac_bams : list, optional
        List of BAM files
        (or one counts table with reads per peak), by default None
    histone_bams : list, optional
        List of H3K27ac ChIP-seq BAM files
        (or one counts table with reads per peak), by default None
    p300_bams : list, optional
        List of p300 ChIP-seq BAM files
        (or one counts table with reads per peak), by default None
    cage_tpms : str, optional
        table with bidirectional regions chr:start-end (columns 1)
        and TPM values (column 2) generated with CAGEfightR, by default None
    columns : list, optional
        List of count table columns to use, by default all columns are used.
    regions : str or list, optional
        BED file or text file with regions, or a list of BED, narrowPeak or
        broadPeak files If None, then the reference regions are used.
    reference : str, optional
        Directory name to a reference.
    factors : list, optional
        List of TF names or file with TFs, one per line. If None (default),
        then all TFs are used.
    genome : str, optional
        Genome name. The default is hg38.
    pfmfile : str, optional
        Motifs in PFM format, with associated motif2factors.txt file.
    pfmscorefile : str, optional
        Path to file with pre-scanned motif scores.
    jaccard_cutoff : int, optional
        minimum jaccard similarity score that is needed to use the model of TF1 for TF2.
        0: any shared motifs, 1: all motifs shared. Default is 0.
    ncore : int, optional
        Number of threads to use. Default is 4.
    """

    if (
        reference is None
        and regions is None
        and pfmscorefile is None
        and cage_tpms is None
    ):
        warnings = [
            "Need either 1) a `reference` directory, 2) one or more `regionfiles`, "
            "or 3) a `pfmscorefile` (a set of pre-scanned regions)!",
            "1) For human data, you can download the REMAP reference here: "
            "https://doi.org/10.5281/zenodo.4768075",
            "2) As regionfiles you can specify one or more BED/narrowPeak files "
            "(with potential enhancer regions, for instance, all ATAC-seq peaks "
            "from your combined experiments)",
            "3) To get pre-scanned regions, run "
            "`gimme scan -Tz --gc -g GENOME REGIONFILES > pfmscorefile.tsv` "
            "(nice if you plan to reuse these regions)",
            "Please see the docs for more info on these options.",
        ]
        for warn in warnings:
            logger.error(warn)
        sys.exit(1)

    if reference is not None and (regions is not None or pfmscorefile is not None):
        logger.error(
            "Need either a reference directory *or* a set of input regions"
            "/pfmscorefile (pfmscorefile and regionfiles can be combined)"
        )
        sys.exit(1)

    # Check if all specified BAM files exist
    _check_input_files(atac_bams, histone_bams, p300_bams, cage_tpms)

    # Check genome, will fail if it is not a correct genome name or file
    Genome(genome)

    if reference is not None:
        if not os.path.exists(reference):
            logger.error(f"Reference directory {reference} does not exist!")
            sys.exit(1)

    ncore = check_cores(ncore)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Read the factors, from a file if needed
    factors = load_tfs(factors)

    # If regions are specified, read them in, combining multiple files if necessary.
    regions = load_regions(regions, genome, outdir)

    p = PeakPredictor(
        reference=reference,
        atac_bams=atac_bams,
        columns=columns,
        histone_bams=histone_bams,
        p300_bams=p300_bams,
        cage_tpms=cage_tpms,
        regions=regions,
        genome=genome,
        pfmfile=pfmfile,
        factors=factors,
        pfmscorefile=pfmscorefile,
        ncore=ncore,
    )

    outfile = os.path.join(outdir, "binding.h5")
    with HDFStore(outfile, "w", complib="lzo", complevel=9) as hdf:

        if p.atac_data is not None:
            hdf.put(key="_atac", value=p.atac_data, format="table")

        if p.histone_data is not None:
            hdf.put(key="_h3k27ac", value=p.histone_data, format="table")

        if p.p300_data is not None:
            hdf.put(key="_p300", value=p.p300_data, format="table")

        if p.cage_data is not None:
            hdf.put(key="_cage", value=p.cage_data, format="table")

        logger.info("Predicting TF activity")
        factor_activity = p.predict_factor_activity()
        hdf.put(key="_factor_activity", value=factor_activity, format="table")

        logger.info("Predicting binding per TF:")
        proba = None
        for factor in sorted(p.f2m):
            try:
                proba = p.predict_binding_probability(
                    factor, jaccard_cutoff=jaccard_cutoff
                )
                hdf.put(
                    key=factor,
                    value=proba.iloc[:, -1].reset_index(drop=True).astype(np.float16),
                    format="table",
                )
            except ValueError as e:
                logger.debug(str(e))

        if proba is None:
            # PeakPredictor.predict_binding_probability() went wrong for *all* TFs
            logger.error(
                "Something went wrong! Please let us know by raising an issue with "
                "your ANANSE version and the log output at "
                "https://github.com/vanheeringen-lab/ANANSE/issues"
            )
            sys.exit(1)
        hdf.put(key="_index", value=proba.index.to_series(), format="table")
