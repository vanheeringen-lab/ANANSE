from glob import glob
import inspect
import os
import re
import sys
from tempfile import NamedTemporaryFile

from fluff.fluffio import load_heatmap_data
from genomepy import Genome
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import scan_regionfile_to_table
from gimmemotifs.moap import moap
import joblib
from loguru import logger
import networkx as nx
import numpy as np
import pandas as pd
from pandas import HDFStore
from sklearn.preprocessing import scale, minmax_scale
from scipy.stats import rankdata
import qnorm

import ananse
from ananse.enhancer_binding import CombineBedFiles
from ananse.utils import get_motif_factors, check_input_factors

# This motif file is not created by default
#   * f"{self.data_dir}/reference.factor.feather"


class PeakPredictor:
    def __init__(
        self,
        reference=None,
        atac_bams=None,
        histone_bams=None,
        regions=None,
        genome="hg38",
        pfmfile=None,
        factors=None,
        pfmscorefile=None,
        ncore=4,
    ):
        self.data_dir = reference

        if atac_bams is None and histone_bams is None:
            raise ValueError("Need either ATAC-seq or H3K27ac BAM file(s).")

        if genome is None:
            logger.warning("Assuming genome is hg38")
            genome = "hg38"
        self.genome = genome
        self.set_species(genome)

        if pfmfile is None and self.species not in ["human", "mouse"]:
            logger.warning(
                f"The genome '{genome}' is not recognized as human or mouse."
            )
            logger.warning(
                "If you do have another species, the motif file likely needs to be adapted."
            )
            logger.warning(
                "Currently mouse and human gene names are used to link motif to TFs."
            )
            logger.warning(
                "If your gene symbols are different, then you will need to create a new mapping"
            )
            logger.warning(
                "and use the `-p` argument. For a possible method to do this, see here:"
            )
            logger.warning(
                "https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors"
            )

        # Set basic information
        self.ncore = ncore
        self._atac_data = None
        self._histone_data = None
        self.factor_models = {}
        self.pfmfile = pfmfile
        self._load_motifs(factors=factors)

        # if the reference regions are used, we can use existing data such
        # as motif scores.
        if regions is None:
            self.region_type = "reference"
            self._load_reference_data()
        # If we have custom regions we have to scan for motifs.
        else:
            self.region_type = "custom"
            self.regions = regions
            if pfmscorefile is None:
                self._scan_motifs(regions)
            else:
                self._load_prescanned_motifs(pfmscorefile)

        # Load ATAC data
        if atac_bams is not None:
            self.load_atac(atac_bams, update_models=False)

        # Load histone ChIP-seq data
        if histone_bams is not None:
            self.load_histone(histone_bams, update_models=False)

        self._set_model_type()

    def _scan_motifs(self, regions):
        """[summary]

        Parameters
        ----------
        regions : [type]
            [description]
        """
        logger.info("Scanning regions for motifs.")
        with NamedTemporaryFile(mode="w") as f:
            print("region", file=f)
            for region in regions:
                print(region, file=f)
            f.flush()
            # TODO: we're still scanning for *all* motifs, even if we only have
            # a few factors
            motif_df = scan_regionfile_to_table(
                f.name, self.genome, "score", ncpus=self.ncore
            )

            self._motifs = pd.DataFrame(index=motif_df.index)
            for factor in self.f2m:
                # if factor not in valid_factors:
                #    continue
                self._motifs[factor] = motif_df[self.f2m[factor]].mean(1)

    def _load_prescanned_motifs(self, pfmscorefile):
        """
        Use pre-scanned gimmemotifs motif scores.

        Parameters
        ----------
        pfmscorefile : str/file
            pre-scanned gimmemotifs scores file
        """
        logger.info("loading pre-scanned motif scores.")

        motif_df = pd.read_table(pfmscorefile, comment="#", index_col=0)
        self._motifs = pd.DataFrame(index=motif_df.index)
        for factor in self.f2m:
            # if factor not in valid_factors:
            #    continue
            self._motifs[factor] = motif_df[self.f2m[factor]].mean(1)

    def _load_reference_data(self):
        """Load data for reference regions.

        Will load three types of data:
        * Motif scores.
        * The average peak coverage (self._avg)
        * The distance from the peak to nearest TSS. (self._dist)

        All of these data are only used with the reference set of regions.
        """
        # Read motifs
        logger.info("loading motifs for reference")
        self._motifs = pd.read_feather(f"{self.data_dir}/reference.factor.feather")
        self._motifs.set_index(self._motifs.columns[0], inplace=True)

        # Read average coverage
        logger.info("loading average peak coverage for reference")
        self._avg = pd.read_table(
            f"{self.data_dir}/reference.coverage.txt",
            sep="\t",
            comment="#",
            index_col=0,
        )
        self._avg.columns = ["average"]
        self._avg["average"] = self._avg["average"] / self._avg["average"].max()

        # Read distance to TSS
        logger.info("loading distance for reference")
        self._dist = pd.read_table(
            f"{self.data_dir}/reference.dist_to_tss.txt",
            sep="\t",
            comment="#",
            index_col=0,
        )

        # Set regions
        self.regions = self._avg.index

    def _load_human_factors(self):
        package_dir = os.path.dirname(ananse.__file__)
        tf_xlsx = os.path.join(package_dir, "db", "lovering.tfs.xlsx")
        valid_factors = pd.read_excel(
            tf_xlsx,
            engine="openpyxl",
            sheet_name=1,
        )
        valid_factors = valid_factors.loc[
            valid_factors["Pseudogene"].isnull(), "HGNC approved gene symbol"
        ].values
        valid_factors = list(set(valid_factors) - set(["EP300"]))
        return valid_factors

    def set_species(self, genome):
        try:
            # Try to get taxonomy id for genomepy managed genome.
            # If there is a taxonomy id, we can be really sure about the species.
            # If genome doesn't have a tax_id, then it will be 'na' and
            # fail to convert to int.
            genome = Genome(genome)
            tax_id = int(genome.tax_id)
            if tax_id == 9606:
                self.species = "human"
            elif tax_id == 10090:
                self.species = "mouse"
            else:
                # tax_id converts to int so it is valid, must be not human or mouse
                self.species = None
            return
        except Exception:
            pass

        mapping = {
            "hg38": "human",
            "hg19": "human",
            "GRCh3": "human",
            "mm10": "mouse",
            "mm9": "mouse",
            "GRCm3": "mouse",
        }

        base_genome = os.path.basename(self.genome.strip("/"))
        for name, species in mapping.items():
            if name in base_genome:
                self.species = species
                return

        self.species = None

    def factors(self):
        if self.species == "human":
            valid_factors = self._load_human_factors()
            return [f for f in self.f2m if f in valid_factors]
        if self.species == "mouse":
            # Mouse mappings are included in the default motif db.
            # Using the fact here that mouse names are not all upper-case.
            # TODO: replace with a curated set of factors.
            return [f for f in self.f2m if f[1:].islower()]
        return list(self.f2m.keys())

    def _load_factor2motifs(self, pfmfile=None, indirect=True, factors=None):
        motifs = read_motifs(pfmfile, as_dict=True)
        f2m = {}

        if self.species == "human":
            valid_factors = self._load_human_factors()

        for name, motif in motifs.items():
            for factor in get_motif_factors(motif, indirect=indirect):
                if factors is not None and factor not in factors:
                    continue

                # TODO: this is temporary, while the motif database we use
                # not very clean...
                if self.species == "human":
                    factor = factor.upper()

                if self.species == "human" and factor not in valid_factors:
                    continue

                f2m.setdefault(factor, []).append(name)
        return f2m

    def _load_motifs(self, indirect=True, factors=None):
        """Load motif-associated data.

        For now, only default motifs are supported.
        Will read factors associated to motifs, and generates a graph of
        related factors based on different factors binding to the same motif.
        This information is used to select the most appropriate TF model.

        Parameters
        ----------
        indirect : bool, optional
            Include TF-motif associations that are not curated, for instance
            based on ChIP-seq motif prediction, or binding inference. This will
            greatly increase TF coverage. By default True.
        """
        if self.pfmfile is None:
            logger.info("using default motif file")
        else:
            logger.debug(f"Motifs: {self.pfmfile}")
        self.motifs = read_motifs(self.pfmfile, as_dict=True)
        self.f2m = self._load_factor2motifs(
            pfmfile=self.pfmfile, indirect=indirect, factors=factors
        )

        if len(self.f2m) == 1:
            logger.info("using motifs for 1 factor")
        else:
            logger.info(f"using motifs for {len(self.f2m)} factors")
        # Create a graph of TFs where edges are determined by the Jaccard index
        # of the motifs that they bind to. For instance, when TF 1 binds motif
        # A and B and TF 2 binds motif B and C, the edge weight will be 0.33.
        tmp_f2m = {}
        if self.pfmfile is not None:
            logger.debug("reading default file")
            tmp_f2m = self._load_factor2motifs(indirect=True)

        for k, v in self.f2m.items():
            if k in tmp_f2m:
                tmp_f2m[k] += v
            else:
                tmp_f2m[k] = v

        self.motif_graph = nx.Graph()
        d = []
        for f1 in tmp_f2m:
            for f2 in tmp_f2m:
                jaccard = len(set(tmp_f2m[f1]).intersection(set(tmp_f2m[f2]))) / len(
                    set(tmp_f2m[f1]).union(set(tmp_f2m[f2]))
                )
                d.append([f1, f2, jaccard])
                if jaccard > 0:
                    self.motif_graph.add_edge(f1, f2, weight=1 - jaccard)

    def _load_bams(self, bams, title, window=200):
        tmp = pd.DataFrame(index=self.regions)
        with NamedTemporaryFile(mode="w") as f_out:

            for region in self.regions:
                print("{}\t{}\t{}".format(*re.split("[:-]", region)), file=f_out)
            f_out.flush()

            for bam in bams:
                result = load_heatmap_data(
                    f_out.name,
                    bam,
                    bins=1,
                    up=window // 2,
                    down=window // 2,
                    rmdup=True,
                    rmrepeats=True,
                )
                tmp[result[0]] = result[2].T[0]

        fname = f"{self.data_dir}/{title}.qnorm.ref.txt.gz"
        if os.path.exists(fname):
            logger.debug(f"quantile normalization for {title}")
            qnorm_ref = pd.read_table(fname, index_col=0)["qnorm_ref"].values
            if len(self.regions) != len(qnorm_ref):
                qnorm_ref = np.random.choice(
                    qnorm_ref, size=len(self.regions), replace=True
                )

            tmp = qnorm.quantile_normalize(tmp, target=qnorm_ref)
        else:
            tmp = np.log1p(tmp)

        # Limit memory usage by using float16
        tmp = tmp.mean(1).astype("float16").to_frame(title)

        fname = f"{self.data_dir}/{title}.mean.ref.txt.gz"
        if self.region_type == "reference" and os.path.exists(fname):
            mean_ref = pd.read_table(fname, index_col=0)
            if mean_ref.shape[0] == tmp.shape[0]:
                mean_ref.index = tmp.index
                tmp[f"{title}.relative"] = (
                    tmp[title] - mean_ref.loc[tmp.index]["mean_ref"].values
                )
                tmp[f"{title}.relative"] = scale(tmp[f"{title}.relative"])
            else:
                logger.debug(f"Regions of {fname} are not the same as input regions.")
                logger.debug("Skipping calculation of relative values.")

        tmp[title] = tmp[title] / tmp[title].max()

        return tmp

    def load_atac(self, bams, update_models=True):
        """Load ATAC-seq counts from BAM files.

        Parameters
        ----------
        bams : list
            List of file names.
        update_models : bool, optional
            Update the model used if data is loaded, by default True.
        """
        logger.info("loading ATAC data")
        self._atac_data = self._load_bams(bams, title="ATAC", window=200)
        if update_models:
            self._set_model_type()

    def load_histone(self, bams, update_models=True):
        """Load H3K27ac ChIP-seq counts from BAM files.

        Parameters
        ----------
        bams : list
            List of file names.
        update_models : bool, optional
            Update the model used if data is loaded, by default True.
        """
        logger.info("loading H3K27ac data")
        self._histone_data = self._load_bams(bams, title="H3K27ac", window=2000)
        if update_models:
            self._set_model_type()

    def _set_model_type(self):
        """Select the mode to use for binding prediction.

        Basically, this will select the columns that are available,
        based on the different types of data that are loaded.
        Reference regions will have the most information.
        """
        cols = ["motif"]
        if self._atac_data is not None:
            cols += ["ATAC"]
            if self.region_type == "reference":
                cols += ["ATAC.relative"]
        if self._histone_data is not None:
            cols += ["H3K27ac"]
        if self.region_type == "reference":
            cols += ["average", "dist"]
        cols = sorted(cols)
        self._X_columns = cols
        self._model_type = "_".join(cols)

        # Load models
        logger.info("Loading models")
        # print(os.path.join(self.data_dir, self._model_type))
        for fname in glob(os.path.join(self.data_dir, self._model_type, "*.pkl")):
            factor = fname.split("/")[-1].replace(".pkl", "")
            self.factor_models[factor] = joblib.load(fname)
        logger.info(f"{len(self.factor_models)} models found")

    def predict_proba(self, factor=None, motifs=None, jaccard_cutoff=0.0):
        """Predict binding probability.

        Predict binding probability for either a TF (factor) or a set of
        motifs. Prediction will be based on the data that been loaded,
        either ATAC-seq or H3K27ac data or both.

        Parameters
        ----------
        factor : str, optional
            Transcription factor name.
        motifs : [type], optional
            Motifs. Currently not implemented.
        jaccard_cutoff : float, optional
            Cutoff for the minimum jaccard overlap between motifs of two TFs for them to be considered related.
            Related TFs can share models. Default = 0.0 (0.1 seems to work well based on subjective testing).
        Returns
        -------
        pandas.DataFrame
            DataFrame with binding probabilities
        """
        if factor is None and motifs is None:
            raise ValueError("Need either a TF name or one or more motifs.")

        if motifs is not None:
            raise NotImplementedError("Custom motifs not yet implemented!")

        if factor not in self.f2m:
            raise ValueError(f"Motif not known for {factor}")

        model, factor = self._load_model(factor, jaccard_cutoff)

        X = self._load_data(factor)
        proba = model.predict_proba(X)[:, 1]

        return pd.DataFrame(proba, index=self.regions)

    def _load_data(self, factor):
        # if self.region_type == "reference":
        # logger.debug("Reading motif data")

        tmp = pd.DataFrame(
            {factor: self._motifs[factor]}, index=self.regions
        )  # pd.read_table(os.path.join(self.data_dir, f"{factor}.motif.txt.gz"), index_col=0)
        # else:
        tmp.columns = ["motif"]
        if self._atac_data is not None:
            tmp = tmp.join(self._atac_data)
        if self._histone_data is not None:
            tmp = tmp.join(self._histone_data)

        if self.region_type == "reference":
            tmp = tmp.join(self._avg)
            tmp = tmp.join(self._dist)
        tmp = tmp.dropna()
        # logger.debug(str(self._X_columns))
        return tmp[self._X_columns]

    def _load_model(self, factor, jaccard_cutoff=0.0):
        """Load TF-binding model that is:
        1. trained for that specific TF
        2. trained on a different TF with a motif overlap of a jacards similarity larger than the cutoff
        3. a general TF binding model if the other options are not available
            Parameters
            ----------
            jaccard_cutoff :
                minimum jacard similarity score that is needed to use the model of TFA for TFB.
        """
        model = None
        max_edge_weight = 1 - jaccard_cutoff
        if factor in self.factor_models:
            logger.info(f"Using {factor} model")
            model = self.factor_models[factor]
        elif factor in self.motif_graph:
            logger.info("Checking for alternative models based on motif overlap...")
            paths = {
                p: v
                for p, v in nx.single_source_dijkstra_path_length(
                    self.motif_graph, factor, cutoff=max_edge_weight
                ).items()
                if p in self.factor_models
            }
            try:
                sub_factor = list(paths.keys())[0]
                logger.info(f"Using {factor} motif with {sub_factor} model weights")
                model = self.factor_models[sub_factor]
                # factor = sub_factor
            except Exception:
                logger.info(f"No match for {factor} based on motifs")
        if model is None:
            logger.info(f"No related TF found for {factor}, using general model")
            model = self.factor_models["general"]

        return model, factor

    def predict_factor_activity(self, nregions=20_000):
        """Predict TF activity.

        Predicted based on motif activity using ridge regression.

        Parameters
        ----------
        """
        # Run ridge regression using motif score to predict (relative) ATAC/H3K27ac signal
        try:
            nregions = int(nregions)
        except ValueError:
            logger.warning("nregions is not an integer, using default number of 20_000")
            nregions = 20_000

        activity = pd.DataFrame()
        for df in (self._atac_data, self._histone_data):
            if df is None:
                continue

            for col in df.columns:
                with NamedTemporaryFile() as f:
                    # float16 will give NaN's
                    signal = df[col].astype("float32")
                    signal = pd.DataFrame({col: scale(signal)}, index=df.index)
                    if df.shape[0] < nregions:
                        signal.to_csv(f.name, sep="\t")
                    else:
                        signal.sample(nregions).to_csv(f.name, sep="\t")
                    try:
                        activity = activity.join(
                            moap(
                                f.name,
                                genome=self.genome,
                                method="bayesianridge",
                                pfmfile=self.pfmfile,
                            ),
                            how="outer",
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


def _check_input_regions(regionfiles, genome, outdir=".", verbose=True, force=False):
    # Load regions from BED or region text file
    if regionfiles is None:
        # Keep regions to None, use reference regions.
        return

    infile = regionfiles[0]
    if len(regionfiles) > 1:
        # merge files, assumed to be all BED
        peak_width = 200
        cbed = CombineBedFiles(genome=genome, peakfiles=regionfiles, verbose=verbose)
        combined_bed = os.path.join(outdir, "regions_combined.bed")
        cbed.run(outfile=combined_bed, width=peak_width, force=force)
        infile = combined_bed

    df = pd.read_table(infile, header=None, sep="\t", comment="#", dtype=str)
    assert df.shape[0] > 2, "regions file must have more that 2 regions."

    test = str(df.at[1, 0])
    if bool(re.match(r"^.*:\d+-\d+$", test)):
        # it's a regions list
        # or it's a Seq2science counts table
        regions = df.iloc[:, 0].tolist()

    elif df.shape[1] >= 3:
        # it's a BED file
        regions = (
            # For Ensembl genome names, make sure it's a string
            df.iloc[:, 0].astype(str)
            + ":"
            + df.iloc[:, 1].astype(str)
            + "-"
            + df.iloc[:, 2].astype(str)
        ).tolist()

    else:
        raise TypeError("Cannot identify regions file(s) type.")

    # remove the header, if any.
    header = str(regions[0])
    if not bool(re.match(r"^.*:\d+-\d+$", header)):
        regions = regions[1:]

    return regions


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
        exit(1)


def predict_peaks(
    outdir,
    atac_bams=None,
    histone_bams=None,
    regionfiles=None,
    reference=None,
    factors=None,
    genome=None,
    pfmfile=None,
    pfmscorefile=None,
    ncore=4,
):
    """Predict binding in a set of genomic regions.

    Binding is predicted based on ATAC-seq and/or H3K27ac ChIP-seq data in
    combination with motif scores. The model that is used is flexible, based
    on the input data. The most accurate model will be the one that uses the
    references regions in combination with both ATAC-seq and H3K27ac ChIP-seq.

    The result will will be saved to an outputfile called `binding.tsv` in the
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
        List of BAM files, by default None
    histone_bams : list, optional
        List of H3K27ac ChIP-seq BAM files, by default None
    regionfiles : list, optional
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
    ncore : int, optional
        Number of threads to use. Default is 4.
    """
    if reference is None and regionfiles is None:
        logger.error("Need either input regions or location of a reference set!")
        logger.error(
            "For human, you can download the REMAP reference here: https://doi.org/10.5281/zenodo.4768075 "
            "(please see the docs on how to install this)."
        )
        logger.error(
            "Otherwise you need to specify one or more BED or narrowPeak files"
        )
        logger.error(
            "with potential enhancer regions, for instance, all ATAC-seq peaks"
        )
        logger.error("from your combined experiments.")
        sys.exit(1)

    if reference is not None and regionfiles is not None:
        logger.error("Need either a reference location *or* or a set of input regions")
        sys.exit(1)

    # Check if all specified BAM files exist
    _check_input_files(atac_bams, histone_bams)

    # Read the factors, from a file if needed
    factors = check_input_factors(factors)

    # Check genome, will fail if it is not a correct genome name or file
    Genome(genome)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # If regions are specified, read them in, combining multiple files if
    # necessary.
    regions = _check_input_regions(regionfiles, genome, outdir=outdir)

    if reference is None:
        install_dir = os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))
        )
        reference = os.path.join(install_dir, "db", "default_reference")

    if reference is not None:
        if not os.path.exists(reference):
            logger.error(f"Reference directory {reference} does not exist!")
            sys.exit(1)

    p = PeakPredictor(
        reference=reference,
        atac_bams=atac_bams,
        histone_bams=histone_bams,
        regions=regions,
        genome=genome,
        pfmfile=pfmfile,
        factors=factors,
        pfmscorefile=pfmscorefile,
        ncore=ncore,
    )

    outfile = os.path.join(outdir, "binding.h5")
    # Make sure we create a new file
    with open(outfile, "w"):
        pass

    with HDFStore(outfile, complib="lzo", complevel=9) as hdf:

        if p._atac_data is not None:
            hdf.put(key="_atac", value=p._atac_data, format="table")

        if p._histone_data is not None:
            hdf.put(key="_h3k27ac", value=p._histone_data, format="table")

        logger.info("Predicting TF activity")
        factor_activity = p.predict_factor_activity()
        hdf.put(key="_factor_activity", value=factor_activity, format="table")

        for factor in p.factors():
            try:
                proba = p.predict_proba(factor)
                hdf.put(
                    key=f"{factor}",
                    value=proba.iloc[:, -1].reset_index(drop=True).astype(np.float16),
                    format="table",
                )

            except ValueError as e:
                logger.debug(str(e))

        hdf.put(key="_index", value=proba.index.to_series(), format="table")
