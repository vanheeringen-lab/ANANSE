from loguru import logger
import pandas as pd
from fluff.fluffio import load_heatmap_data
from tempfile import NamedTemporaryFile
import re
import seaborn as sns
import qnorm
import numpy as np
from glob import glob
import os
import joblib
from sklearn.preprocessing import scale
import networkx as nx
from gimmemotifs.motif import read_motifs

# Note: there are currently 3 paths hardcoded:
#  data_dir, _avg and _dist
# This should be changed to make it general
# The motif files are also not created by default
#   * f"{self.data_dir}/reference.factor.feather"
#   * f"{factor}.motif.txt.gz"


class PeakPredictor():
    #_reference_regions = pd.read_table("/scratch/heeringen/atac_ananse/models/remap2018.enhancers.w200.bed.v1.0")
    data_dir = "/ceph/rimlsfnwi/data/moldevbio/heeringen/heeringen/atac_ananse/models/saved"
    
    def __init__(self, atac_bams=None, histone_bams=None, regions=None, genome="hg38"):
        if atac_bams is None and histone_bams is None:
            raise ValueError("Need either ATAC-seq or H3K27ac BAM file(s).")
            
        self.genome = genome
        
        self._atac_data = None
        self._histone_data = None
        self.factor_models = {}
        
        if regions is None:
            self._load_reference_data()
            self.region_type = "reference"
        else:
            self.region_type = "custom"
            self.regions = regions
            self._scan_motifs(regions)
        
        if atac_bams is not None:
            self.load_atac(atac_bams, update_models=False)
        
        if histone_bams is not None:
            self.load_histone(histone_bams, update_models=False)
        
        self._set_model_type()
        self._load_motifs()
            
    def _scan_motifs(self, regions):
        with NamedTemporaryFile() as f:
            for region in regions:
                print(region, file=f)
            f.flush()
            motif_df = scan_regionfile_to_table(f.name, genome, "score")

            self._motifs = pd.DataFrame(index=motif_df.index)
            for factor in p.f2m:
                if factor not in valid_factors:
                    continue
                self._motifs[factor] = motif_df[p.f2m[factor]].mean(1)
        
    def _load_reference_data(self):
        # Read motifs
        logger.info("loading motifs for reference")
        self._motifs = pd.read_feather(f"{self.data_dir}/reference.factor.feather")
        self._motifs.set_index(self._motifs.columns[0], inplace=True)
        
        # Read average coverage
        logger.info("loading average peak coverage for reference")
        self._avg = pd.read_table(f"{self.data_dir}/reference.coverage.txt", sep="\t", comment="#", index_col=0)
        self._avg.columns = ["average"]
        self._avg["average"] = self._avg["average"] / self._avg["average"].max()
        
        # Read distance to TSS
        logger.info("loading distance for reference")
        self._dist = pd.read_table(f"{self.data_dir}/reference.dist_to_tss.txt", sep="\t", comment="#", index_col=0)

        # Set regions
        self.regions = self._avg.index
    
    
    def _load_motifs(self, indirect=True):
        self.motifs = read_motifs(as_dict=True)

        self.f2m = {}
        for name, motif in self.motifs.items():
            for k, factors in motif.factors.items():
                if k != "direct" and not indirect:
                    print("skip")
                    continue
                for factor in factors:
                    self.f2m.setdefault(factor.upper(), []).append(name)
        
        self.motif_graph = nx.Graph()
        d = []
        for f1 in self.f2m:
            for f2 in self.f2m:
                jaccard = len(set(self.f2m[f1]).intersection(set(self.f2m[f2]))) / len(set(self.f2m[f1]).union(set(self.f2m[f2])))
                d.append([f1,f2,jaccard])
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
                qnorm_ref = np.random.choice(qnorm_ref, size=len(self.regions), replace=True)

            tmp = qnorm.quantile_normalize(tmp, target=qnorm_ref)
        else: 
            tmp = np.log1p(tmp)
            
        tmp = tmp.mean(1).to_frame(title)
        
        fname = f"{self.data_dir}/{title}.mean.ref.txt"
        if os.path.exists(fname):
            mean_ref = pd.read_table(fname, index_col=0)
            tmp[f"{title}.relative"] = tmp[title] - mean_ref.loc[tmp.index]["mean_ref"]
            tmp[f"{title}.relative"] = scale(tmp[f"{title}.relative"])
        
        tmp[title] = tmp[title] / tmp[title].max()
        
        return tmp
    
    def load_atac(self, bams, update_models=True):
        logger.info("loading ATAC data")
        self._atac_data = self._load_bams(bams, title="ATAC", window=200)
        if update_models:        
            self._set_model_type()
        
    def load_histone(self, bams, update_models=True):
        logger.info("loading H3K27ac data")
        self._histone_data = self._load_bams(bams, title="H3K27ac", window=2000)
        if update_models:
            self._set_model_type()
        
    def _set_model_type(self):
        cols = ["motif"]
        if self._atac_data is not None:
            cols += ["ATAC", "ATAC.relative"]
        if self._histone_data is not None:
            cols += ["H3K27ac"]
        if self.region_type == "reference":
            cols += ['average', 'dist']
        cols = sorted(cols)
        self._X_columns = cols
        self._model_type = "_".join(cols)
        
        # Load models
        logger.info("Loading models")
        #print(os.path.join(self.data_dir, self._model_type))
        for fname in glob(os.path.join(self.data_dir, self._model_type, "*.pkl")):
            factor = fname.split("/")[-1].replace(".pkl", "")
            self.factor_models[factor] = joblib.load(fname)
        logger.info(f"{len(self.factor_models)} models found")
        
    def predict_proba(self, factor=None, motifs=None):
        if factor is None and motifs is None:
            raise ValueError("Need either a TF name or one or more motifs.")
        
        if motifs is not None:
            raise NotImplementedError("Custom motifs not yet implemented!")
        
        if factor not in self.f2m:
            raise ValueError(f"Motif not known for {factor}")
        
        model, factor = self._load_model(factor)
        
        X = self._load_data(factor)
        proba = model.predict_proba(X)[:,1]
        
        return pd.DataFrame(proba, index=self.regions)
    
    def _load_data(self, factor):
        #if self.region_type == "reference":
        logger.info("Reading motif data")
        
        tmp = pd.DataFrame({factor:self._motifs[factor]}, index=self.regions)#pd.read_table(os.path.join(self.data_dir, f"{factor}.motif.txt.gz"), index_col=0)
        #else:
        tmp.columns = ["motif"]
        if self._atac_data is not None:
            tmp = tmp.join(self._atac_data)
        if self._histone_data is not None:
            tmp = tmp.join(self._histone_data)

        if self.region_type == "reference":
            tmp = tmp.join(self._avg)
            tmp = tmp.join(self._dist)
        tmp = tmp.dropna()
        logger.info(str(self._X_columns))
        return tmp[self._X_columns]
    
    def _load_model(self, factor):
        model = None
        if factor in self.factor_models:
            model = self.factor_models[factor]
        elif factor in self.motif_graph:
            paths = {p:v for p,v in nx.single_source_dijkstra_path_length(self.motif_graph, factor).items() if p in self.factor_models}
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
        engine='openpyxl', sheet_name=1)
    valid_factors = valid_factors.loc[valid_factors["Pseudogene"].isnull(), "HGNC approved gene symbol"].values
    valid_factors = [f for f in valid_factors if f != "EP300"]
    return valid_factors

    
def predict_peaks(outfile, atac_bams=None, histone_bams=None, regions=None, factors=None):
    
    if factors is None:
        factors = load_default_factors()
    p = PeakPredictor(
        atac_bams=atac_bams,
        histone_bams=histone_bams,
        regions=regions,
    )

    with open(outfile, "w") as f:
        pass

    with open(outfile, "a") as f:
        print("factor\tenhancer\tbinding", file=f)
    
    
    print(atac_bams)
    print(histone_bams)
    print(outfile)
    
    for factor in factors: 
        try:
            proba = p.predict_proba(factor)
            proba = proba.reset_index()
            proba.columns = ["enhancer", "binding"]
            proba["factor"] = factor
            proba[["factor", "enhancer", "binding"]].to_csv(f, index=False, header=False, sep="\t", float_format='%.5f')
        except ValueError as e:
            print(str(e))
        
