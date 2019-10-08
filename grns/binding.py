#!/usr/bin/env python
import argparse
import os
import pickle
import subprocess
import sys
import math
import ast
import warnings
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
from sklearn import preprocessing
from chest import Chest
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from pybedtools import BedTool
from genomepy import Genome

from gimmemotifs.config import MotifConfig
from gimmemotifs.scanner import Scanner
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import as_fasta
from gimmemotifs.background import create_random_genomic_bedfile
from gimmemotifs.scanner import scan_to_best_match
from gimmemotifs.background import MatchedGcFasta

from grns import mytmpdir

warnings.filterwarnings('ignore')

class Binding(object):

    def __init__(self, genome="hg19", gene_bed=None, pwmfile=None):
        
        self.genome=genome
        g = Genome(self.genome)
        self.gsize = g.props["sizes"]["sizes"]

        self.motifs2factors = pwmfile.replace(".pfm", ".motif2factors.txt")
        self.factortable= pwmfile.replace(".pfm", ".factortable.txt")

        self.gene_bed=gene_bed

        # read motifs
        with open(pwmfile) as pwm_in:
            motifs = read_motifs(pwm_in)

        self.pwmfile = NamedTemporaryFile(mode="w", dir=mytmpdir(),delete=False)
        for motif in motifs:
            if motif.factors:
                self.pwmfile.write("{}\n".format(motif.to_pwm()))
        
        self.model = "../db/dream_model.txt"

    def clear_peaks(self, peak_bed, filter_promoter=True, up=2000, down=2000):
        """Filter the enhancer peaks in promoter range.

        all overlap Enh-TSS(up2000 to down2000) pair
        """

        peaks = BedTool(peak_bed)
        b = BedTool(self.gene_bed)
        b = b.flank(l=1, r=0, s=True, g=self.gsize).slop(l=up, r=down, g=self.gsize, s=True)
        vals = []
        # print(peaks)
        # for f in b.intersect(peaks, wo=True, nonamecheck=True):
        for f in b.intersect(peaks, wo=True):
            chrom = f[0]
            gene = f[3]
            peak_start, peak_end = int(f[13]), int(f[14])
            vals.append(chrom+":"+str(peak_start)+"-"+str(peak_end))
        fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(),delete=False)
        with open(peak_bed) as pbed:
            for line in pbed:
                if filter_promoter:
                    if line.split()[0]+":"+line.split()[1]+"-"+line.split()[2] not in vals:
                        fl2.write(line)
                else:
                    fl2.write(line)
        return(fl2.name)

    def get_motif_distribution(self, normalize='gcbins', nregions=10000, length=200, pwmfile=None, force=False):
        """Calculate mean and sd of motif scores.""" 
        # if os.path.exists(outname) and not force:
        #     sys.stderr.write("File {} already exists!\n".format(outname))
        #     sys.stderr.write("Set force to True to overwrite!\n".format(outname))
        #     return

        if normalize == 'gcbins':
            # Create bed file with random regions
            config = MotifConfig()
            tmp = NamedTemporaryFile(suffix=".bed")
            create_random_genomic_bedfile(tmp.name, self.genome, length, nregions)
            
            if pwmfile is None:
                params = config.get_default_params()
                pwmfile = os.path.join(config.get_motif_dir(), params["motif_db"])
            result = scan_to_best_match(tmp.name, pwmfile, genome=self.genome, score=True)

        elif normalize == 'gc':
            bg = MatchedGcFasta(fin_regions_fa, genome=self.genome, number=nregions)
            if pwmfile is None:
                params = config.get_default_params()
                pwmfile = os.path.join(config.get_motif_dir(), params["motif_db"])
            result = scan_to_best_match(bg, pwmfile, genome=self.genome, score=True)

        fl3 = NamedTemporaryFile(mode="w", dir=mytmpdir())
        fl3.write("motif\tmean\tstd\n")
        for motif, scores in result.items():
            fl3.write("{}\t{}\t{}\n".format(motif, np.mean(scores), np.std(scores)))

    def get_peakRPKM(self, fin_rpkm):
        peaks = pd.read_table(fin_rpkm, names=["chrom", "start", "end", "peakRPKM"])
        peaks["peak"] = peaks["chrom"] + ":" + peaks["start"].astype(str) +  "-" + peaks["end"].astype(str)
        add = peaks["peakRPKM"][peaks["peakRPKM"] > 0].min()
        peaks["log10_peakRPKM"] = np.log10(peaks["peakRPKM"] + add)
        peaks["peakRPKMScale"] = minmax_scale(peaks["log10_peakRPKM"])
        peaks["peakRPKMRank"] = minmax_scale(rankdata(peaks["log10_peakRPKM"]))
        
        cols = ["peak", "peakRPKM", "log10_peakRPKM", "peakRPKMScale", "peakRPKMRank"]
        # outname = os.path.join(outdir, "peakRPKM.txt")
        peakrpkmfile = NamedTemporaryFile(mode="w", dir=mytmpdir(),delete=False)
        peaks[cols].to_csv(peakrpkmfile, sep="\t", index=False)
        return(peakrpkmfile.name)

    def get_PWMScore(self, fin_regions_fa, normalize='gcbins'):

        pwmscorefile=NamedTemporaryFile(mode="w", dir=mytmpdir(),delete=False)
        if normalize == 'gcbins':
            seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=self.genome).ids]    
            s = Scanner()
            s.set_motifs(self.pwmfile.name)
            s.set_threshold(threshold=0.0)
            s.set_genome(self.genome)
            with open(self.pwmfile.name) as f:
                motifs = read_motifs(f)
            
            chunksize = 10000
            for chunk in range(0, len(seqs), chunksize):
                chunk_seqs = seqs[chunk:chunk+chunksize]
                print(chunk, "-", chunk + chunksize)
                pwm_score = []
                it = s.best_score(chunk_seqs, zscore=True, gc=True)
                for seq,scores in zip(chunk_seqs, it):
                    for motif, score in zip(motifs, scores):
                        pwm_score.append([motif.id, seq, score])
                pwm_score = pd.DataFrame(pwm_score, columns=["motif", "enhancer", "zscore"])
                pwm_score = pwm_score.set_index("motif")
            
                print("Combine")
                pwm_score["zscoreRank"] = minmax_scale(rankdata(pwm_score["zscore"]))
                cols = ["enhancer", "zscore", "zscoreRank"]
                write_header = False
                if chunk == 0:
                    write_header = True
                pwm_score[cols].to_csv(pwmscorefile, sep='\t', header=write_header)
        elif normalize == 'gc':
            motifsdic={}
            for m in motifs:
                motifsdic[m.id]=len(m)
            # Filter for motifs that have factors assigned

            print("Motif background distribution")
            motif_stats = os.path.join(outdir, "motif_distribution.txt")
            comput_peak_background(fin_regions_fa, motif_stats, self.genome, pwmfile=pwmfile)
            motif_bg = pd.read_table(motif_stats, index_col=0)
            
            motiflen=[]
            for m in motif_bg.index:
                motiflen.append(motifsdic[m])
            motif_bg["len"]=np.log2(motiflen)
            motif_bg.to_csv(motif_stats,sep="\t")    
            print("Motif scan")
            seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=self.genome).ids]
            s = Scanner()
            s.set_motifs(pwmfile)
            s.set_threshold(threshold=0.0)
            s.set_genome(self.genome)
            with open(pwmfile) as f:
                motifs = read_motifs(f)
            
            chunksize = 10000
            for chunk in range(0, len(seqs), chunksize):
                chunk_seqs = seqs[chunk:chunk+chunksize]
                print(chunk, "-", chunk + chunksize)
                pwm_score = []
                it = s.best_score(chunk_seqs)
                for seq,scores in zip(chunk_seqs, it):
                    for motif, score in zip(motifs, scores):
                        pwm_score.append([motif.id, seq, score])
                pwm_score = pd.DataFrame(pwm_score, columns=["motif", "enhancer", "score"])
                pwm_score = pwm_score.set_index("motif")
            
                print("Combine")
                pwm_score = pwm_score.join(pd.read_table(motif_stats, index_col=0))
                # scale motif score
                pwm_score["zscore"] = (pwm_score["score"] - pwm_score["mean"]) / pwm_score["std"] 
                pwm_score["zscoreRank"] = minmax_scale(rankdata(pwm_score["zscore"]))
                cols = ["enhancer", "score", "zscore", "zscoreRank"]
                write_header = False
                if chunk == 0:
                    write_header = True
                pwm_score[cols].to_csv(pwmscorefile, sep='\t', header=write_header)
        elif normalize == 'nogc' :
            print("Motif background distribution")
            motif_stats = os.path.join(outdir, "motif_distribution.txt")
            get_motif_distribution(motif_stats, self.genome, pwmfile=pwmfile)
            motif_bg = pd.read_table(motif_stats, index_col=0)
            
            print("Motif scan")
            seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=self.genome).ids]
            
            s = Scanner()
            s.set_motifs(pwmfile)
            s.set_threshold(threshold=0.0)
            s.set_genome(self.genome)
            with open(pwmfile) as f:
                motifs = read_motifs(f)
            
            chunksize = 10000
            for chunk in range(0, len(seqs), chunksize):
                chunk_seqs = seqs[chunk:chunk+chunksize]
                print(chunk, "-", chunk + chunksize)
                pwm_score = []
                it = s.best_score(chunk_seqs)
                for seq,scores in zip(chunk_seqs, it):
                    for motif, score in zip(motifs, scores):
                        pwm_score.append([motif.id, seq, score])
        
                pwm_score = pd.DataFrame(pwm_score, columns=["motif", "enhancer", "score"])
                pwm_score = pwm_score.set_index("motif")
            
                print("Combine")
                pwm_score = pwm_score.join(pd.read_table(motif_stats, index_col=0))
                # scale motif score
                pwm_score["zscore"] = (pwm_score["score"] - pwm_score["mean"]) / pwm_score["std"]
                pwm_score["zscoreRank"] = minmax_scale(rankdata(pwm_score["zscore"]))
                cols = ["enhancer", "score", "zscore", "zscoreRank"]
                write_header = False
                if chunk == 0:
                    write_header = True
                pwm_score[cols].to_csv(pwmscorefile, sep='\t', header=write_header)
        
        return(pwmscorefile.name)

    def get_binding_score(self, pwm, peak):
        
        # Load model
        with open(self.model, "rb") as f:
            clf = pickle.load(f)

        ft = dd.read_csv(self.factortable, sep="\t") 

        r = pwm.merge(peak, left_on="enhancer", right_on="peak")[["motif", "enhancer",  "zscore", "peakRPKMScale"]]

        r = r.merge(ft, left_on="motif", right_on="motif")
        r = r.groupby(["factor", "enhancer"])[["zscore", "peakRPKMScale"]].max()
        r = r.dropna().reset_index()

        cache = Chest(available_memory=20e9)
        print("combining tables")
        table = r.compute()
        print("predict")
        table["binding"] = clf.predict_proba(table[["zscore", "peakRPKMScale"]])[:,1]
        print("save results")

        return(table)

    def run_binding(self, peak_bed, outdir):
        filter_bed = self.clear_peaks(peak_bed)

        pwm_weight = self.get_PWMScore(filter_bed)
        pwm = dd.read_csv(pwm_weight, sep="\t")

        peak_weight = self.get_peakRPKM(filter_bed)
        peak = dd.read_csv(peak_weight, sep="\t")

        table=self.get_binding_score(pwm, peak)

        outfile = os.path.join(outdir, "binding.predicted.txt")
        table.to_csv(outfile, sep="\t", index=False)

