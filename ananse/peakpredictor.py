from glob import glob
import os
import re
from tempfile import NamedTemporaryFile

from fluff.fluffio import load_heatmap_data
from genomepy import Genome
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import scan_regionfile_to_table
import joblib
from loguru import logger
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
import qnorm

from ananse.enhancer_binding import CombineBedFiles

# Note: there are currently 3 paths hardcoded:
#  data_dir, _avg and _dist
# This should be changed to make it general
# The motif files are also not created by default
#   * f"{self.data_dir}/reference.factor.feather"
#   * f"{factor}.motif.txt.gz"


class PeakPredictor:
    def __init__(
        self,
        data_dir="/ceph/rimlsfnwi/data/moldevbio/heeringen/heeringen/atac_ananse/models/saved",
        atac_bams=None,
        histone_bams=None,
        regions=None,
        genome="hg38",
    ):
        self.data_dir = data_dir

        if atac_bams is None and histone_bams is None:
            raise ValueError("Need either ATAC-seq or H3K27ac BAM file(s).")

        if genome is None:
            logger.info("Assuming genome is hg38")
            genome = "hg38"

        # Set basic information
        self.genome = genome
        self._atac_data = None
        self._histone_data = None
        self.factor_models = {}
        self._load_motifs()

        # if the reference regions are used, we can use existing data such
        # as motif scores.
        if regions is None:
            self.region_type = "reference"
            self._load_reference_data()
        # If we have custom regions we have to scan for motifs.
        else:
            self.region_type = "custom"
            self.regions = regions
            self._scan_motifs(regions)

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

            motif_df = scan_regionfile_to_table(f.name, self.genome, "score")

            self._motifs = pd.DataFrame(index=motif_df.index)
            for factor in self.f2m:
                # if factor not in valid_factors:
                #    continue
                self._motifs[factor] = motif_df[self.f2m[factor]].mean(1)

    def _load_reference_data(self):
        """Load data for reference regions.

        Will load three types of data:
        * Motif scores.
        * The average peak coverage.
        * The distance from the peak to nearest TSS.
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

    def _load_motifs(self, indirect=True):
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
        self.motifs = read_motifs(as_dict=True)

        self.f2m = {}
        for name, motif in self.motifs.items():
            for k, factors in motif.factors.items():
                if k != "direct" and not indirect:
                    print("skip")
                    continue
                for factor in factors:
                    self.f2m.setdefault(factor.upper(), []).append(name)

        # Create a graph of TFs where edges are determined by the Jaccard index
        # of the motifs that they bind to. For instance, when TF 1 binds motif
        # A and B and TF 2 binds motif B and C, the edge weight will be 0.33.
        self.motif_graph = nx.Graph()
        d = []
        for f1 in self.f2m:
            for f2 in self.f2m:
                jaccard = len(set(self.f2m[f1]).intersection(set(self.f2m[f2]))) / len(
                    set(self.f2m[f1]).union(set(self.f2m[f2]))
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

        fname = f"{self.data_dir}/{title}.qnorm.ref.txt"
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

        tmp = tmp.mean(1).to_frame(title)

        fname = f"{self.data_dir}/{title}.mean.ref.txt"
        if self.region_type == "reference" and os.path.exists(fname):
            mean_ref = pd.read_table(fname, index_col=0)
            tmp[f"{title}.relative"] = tmp[title] - mean_ref.loc[tmp.index]["mean_ref"]
            tmp[f"{title}.relative"] = scale(tmp[f"{title}.relative"])

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
        Reference regions will have the mmost information.
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

    def predict_proba(self, factor=None, motifs=None):
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

        model, factor = self._load_model(factor)

        X = self._load_data(factor)
        proba = model.predict_proba(X)[:, 1]

        return pd.DataFrame(proba, index=self.regions)

    def _load_data(self, factor):
        # if self.region_type == "reference":
        logger.debug("Reading motif data")

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
        logger.debug(str(self._X_columns))
        return tmp[self._X_columns]

    def _load_model(self, factor):
        model = None
        if factor in self.factor_models:
            model = self.factor_models[factor]
        elif factor in self.motif_graph:
            paths = {
                p: v
                for p, v in nx.single_source_dijkstra_path_length(
                    self.motif_graph, factor
                ).items()
                if p in self.factor_models
            }
            try:
                sub_factor = list(paths.keys())[0]
                logger.info(f"Using {sub_factor} model for {factor}")
                model = self.factor_models[sub_factor]
                factor = sub_factor
            except Exception:
                logger.info(f"No match for {factor} based on motifs")
        if model is None:
            logger.info("Using general model")
            model = self.factor_models["general"]

        return model, factor


def load_default_factors():
    valid_factors = pd.read_excel(
        "https://www.biorxiv.org/content/biorxiv/early/2020/12/07/2020.10.28.359232/DC1/embed/media-1.xlsx",
        engine="openpyxl",
        sheet_name=1,
    )
    valid_factors = valid_factors.loc[
        valid_factors["Pseudogene"].isnull(), "HGNC approved gene symbol"
    ].values
    valid_factors = [f for f in valid_factors if f != "EP300"]
    return valid_factors


def _check_input_factors(factors):
    """Check factors.

    Factors can eiher be a list of transcription factors, or a filename of a
    file that containts TFs. Returns a list of factors.
    If factors is None, it will return the defaylt transcription factors.

    Returns
    -------
    list
        List of TF names.
    """
    # Load factors
    if factors is None:
        factors = load_default_factors()
    elif isinstance(factors, str) or (len(factors) == 1 and os.path.exists(factors[0])):
        factors = [line.strip() for line in open(factors[0])]
    return factors


def _check_input_regions(regionfiles, genome, outdir=".", verbose=True, force=False):
    # Load regions from BED or region text file
    if regionfiles is None:
        # Keep regions to None, use reference regions.
        return

    infile = regionfiles[0]
    if len(regionfiles) > 1:
        # merge files
        peak_width = 200
        cbed = CombineBedFiles(genome=genome, peakfiles=regionfiles, verbose=verbose)
        combined_bed = os.path.join(outdir, "regions_combined.bed")
        cbed.run(outfile=combined_bed, width=peak_width, force=force)
        infile = combined_bed

    df = pd.read_table(infile, header=None)
    if df.shape[1] >= 3:
        regions = (
            df.iloc[:, 0].astype(
                str
            )  # For Ensembl genome names, make sure it's a string
            + ":"
            + df.iloc[:, 1].astype(str)
            + "-"
            + df.iloc[:, 2].astype(str)
        ).tolist()
    else:
        regions = df.iloc[:, 0].tolist()
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
    factors=None,
    genome=None,
):
    """Predict binding in a set of genomic regions.

    Binding is predicted based on ATAC-seq and/or H3K27ac ChIP-seq data in
    combination with motif scores. The model that is used is flexible, based
    on the input data. The most accurate model will be the one that uses the
    references regions in combination with both ATAC-seq and H3K27ac ChIP-seq.

    Parameters
    ----------
    outfile : str
        Name of output file.
    atac_bams : list, optional
        List of BAM files, by default None
    histone_bams : [type], optional
        List of H3K27ac ChIP-seq BAM files, by default None
    regionfiles : list, optional
        BED file or text file with regions, or a list of BED, narrowPeak or
        broadPeak files If None, then the reference regions are used.
    factors : list, optional
        List of TF names or file with TFs, one per line. If None (default),
        then all TFs are used.
    genome : str, optional
        Genome name. The default is hg38.
    """
    # Check if all specified BAM files exist
    _check_input_files(atac_bams, histone_bams)

    # Read the factors, from a file if needed
    factors = _check_input_factors(factors)

    # Check genome, will fail if it is not a correct genome name or file
    Genome(genome)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # If regions are specified, read them in, combining multiple files if
    # necessary.
    regions = _check_input_regions(regionfiles, genome, outdir=outdir)

    p = PeakPredictor(
        atac_bams=atac_bams,
        histone_bams=histone_bams,
        regions=regions,
        genome=genome,
    )

    outfile = os.path.join(outdir, "binding.tsv")

    # Make sure we create a new file
    with open(outfile, "w") as f:
        pass

    with open(outfile, "a") as f:
        print("factor\tenhancer\tbinding", file=f)

        for factor in factors:
            try:
                proba = p.predict_proba(factor)
                proba = proba.reset_index()
                proba.columns = ["enhancer", "binding"]
                proba["factor"] = factor
                proba[["factor", "enhancer", "binding"]].to_csv(
                    f, index=False, header=False, sep="\t", float_format="%.5f"
                )
            except ValueError as e:
                logger.debug(str(e))
