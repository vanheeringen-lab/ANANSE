#!/usr/bin/env python
import argparse
import os
import pickle
import subprocess
import sys
import math
import ast

import warnings
warnings.filterwarnings('ignore')
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


def filter_promoter_peaks(gene_bed, peak_bed, fpomoter=True, up=2000, down=2000):
    #all overlap Enh-TSS(up8000 to down2000) pair
    g = Genome(genome)
    gsize = g.props["sizes"]["sizes"]
    
    peaks = BedTool(peak_bed)
    b = BedTool(gene_bed)
    b = b.flank(l=1, r=0, s=True, g=gsize).slop(l=up, r=down, g=gsize, s=True)
    vals = []
    for f in b.intersect(peaks, wo=True, nonamecheck=True):
        chrom = f[0]
        gene = f[3]
        peak_start, peak_end = int(f[13]), int(f[14])
        vals.append(chrom+":"+str(peak_start)+"-"+str(peak_end))
    fl2=open(os.path.join(outdir, "filtered_enahncers.txt"),"w")
    with open(peak_bed) as pbed:
        for line in pbed:
            if fpomoter:
                if line.split()[0]+":"+line.split()[1]+"-"+line.split()[2] not in vals:
                    fl2.write(line)
            else:
                fl2.write(line)

def get_motif_distribution(outname, genome, nregions=10000, length=200, pwmfile=None, force=False):
    """Calculate mean and sd of motif scores.""" 
    if os.path.exists(outname) and not force:
        sys.stderr.write("File {} already exists!\n".format(outname))
        sys.stderr.write("Set force to True to overwrite!\n".format(outname))
        return

    # Create bed file with random regions
    config = MotifConfig()
    tmp = NamedTemporaryFile(suffix=".bed")
    create_random_genomic_bedfile(tmp.name, genome, length, nregions)
    
    if pwmfile is None:
        params = config.get_default_params()
        pwmfile = os.path.join(config.get_motif_dir(), params["motif_db"])
    result = scan_to_best_match(tmp.name, pwmfile, genome=genome, score=True)
    
    with open(outname, "w") as f:
        f.write("motif\tmean\tstd\n")
        for motif, scores in result.items():
            f.write("{}\t{}\t{}\n".format(motif, np.mean(scores), np.std(scores)))

def comput_peak_background(fin_regions_fa, outname, genome, nregions=10000, length=200, pwmfile=None, force=False):
    """Calculate mean and sd of motif scores.""" 
    if os.path.exists(outname) and not force:
        sys.stderr.write("File {} already exists!\n".format(outname))
        sys.stderr.write("Set force to True to overwrite!\n".format(outname))
        return

    # seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=genome).ids]
    bg = MatchedGcFasta(fin_regions_fa, genome=genome, number=nregions)
    if pwmfile is None:
        params = config.get_default_params()
        pwmfile = os.path.join(config.get_motif_dir(), params["motif_db"])

    result = scan_to_best_match(bg, pwmfile, genome=genome, score=True)
    with open(outname, "w") as f:
        f.write("motif\tmean\tstd\n")
        for motif, scores in result.items():
            f.write("{}\t{}\t{}\n".format(motif, np.mean(scores), np.std(scores)))

def get_PWMScore(fin_regions_fa, fin_pwm, outdir, genome="hg19", gcbins=True):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # read motifs
    with open(fin_pwm) as pwm_in:
        motifs = read_motifs(pwm_in)
    motifsdic={}
    for m in motifs:
        motifsdic[m.id]=len(m)
    # Filter for motifs that have factors assigned

    pwmfile = os.path.join(outdir, "filtered_motifs.pwm")
    with open(pwmfile, "w") as fout:
        for motif in motifs:
            if motif.factors:
                fout.write("{}\n".format(motif.to_pwm()))
    if gcbins:
        seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=genome).ids]    
        s = Scanner()
        s.set_motifs(pwmfile)
        s.set_threshold(threshold=0.0)
        s.set_genome(genome)
        with open(pwmfile) as f:
            motifs = read_motifs(f)
        
        chunksize = 10000
        with open('{}/pwmScore.txt'.format(outdir), "w") as fout:
            for chunk in range(0, len(seqs), chunksize):
                chunk_seqs = seqs[chunk:chunk+chunksize]
                print(chunk, "-", chunk + chunksize)
                pwm_score = []
                it = s.best_score(chunk_seqs,zscore=True,gc=True)
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
                pwm_score[cols].to_csv(fout, sep='\t', header=write_header)
    else:
        print("Motif background distribution")
        motif_stats = os.path.join(outdir, "motif_distribution.txt")
        comput_peak_background(fin_regions_fa, motif_stats, genome, pwmfile=pwmfile)
        motif_bg = pd.read_table(motif_stats, index_col=0)
        
        motiflen=[]
        for m in motif_bg.index:
            motiflen.append(motifsdic[m])
        motif_bg["len"]=np.log2(motiflen)
        motif_bg.to_csv(motif_stats,sep="\t")    
        print("Motif scan")
        seqs = [s.split(" ")[0] for s in as_fasta(fin_regions_fa, genome=genome).ids]
        s = Scanner()
        s.set_motifs(pwmfile)
        s.set_threshold(threshold=0.0)
        s.set_genome(genome)
        with open(pwmfile) as f:
            motifs = read_motifs(f)
        
        chunksize = 10000
        with open('{}/pwmScore.txt'.format(outdir), "w") as fout:
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
                pwm_score[cols].to_csv(fout, sep='\t',
                        header=write_header)

def get_expression(fin_expression, features, outdir, min_tpm=1e-10, column="tpm"):
    df = pd.read_hdf(features)
    df=df[["source_target", "factor", "gene"]]
    df.source_target=[i.upper() for i in list(df.source_target)]
    df.gene=[i.upper() for i in list(df.gene)]
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
            pd.concat([
                pd.read_table(f, index_col=0)[[column]] for f in fin_expression],
                axis=1).mean(1), 
            columns=[column])
    expression.index=[i.upper() for i in list(expression.index)]
   # print(expression)
    expression[column] = np.log2(expression[column] + 1e-5)
    df = df.join(expression, on="factor")
    df = df.rename(columns={column:"factor_expression"})
    df = df.join(expression, on="gene")
    df = df.rename(columns={column:"target_expression"})

    df = df.dropna()

    for col in ["factor_expression", "target_expression"]:
        df[col + ".scale"] = minmax_scale(df[col])
        df[col + ".rank.scale"] = minmax_scale(rankdata(df[col]))

    outfile = os.path.join(outdir, "expression.txt")
    df.to_csv(outfile, sep="\t")

def get_factorExpression(fin_expression, motifs2factors, outdir):
    import numpy as np
    import pandas as pd
    from scipy.stats import rankdata
    from sklearn import preprocessing
    import warnings
    warnings.filterwarnings('ignore')
    factorsExpression = {}
    #for line in open('/home/george/data/cis-bp.vertebrate.clusters.v3.0.motif2factors.txt'):
    for line in open(motifs2factors):
        motif = line.split('\t')[0].upper()
        if not line.split('\t')[1].strip().split(',') == ['']:
            for factor in line.split('\t')[1].strip().split(','):
                factorsExpression[factor.upper()] = []

    for f in fin_expression: 
        with open(f) as fa:
            for line in fa:
                if not line.startswith('target_id'):
                    gene = line.split('\t')[0].upper()
                    expression = float(line.split('\t')[1])
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
    with open(os.path.join(outdir, "factorExpression0.txt"), 'w') as fout:
        fout.write('#factor\tfactorExpression\n')
        for factor in factorsExpression:
            if len(factorsExpression[factor]) == 0:
                fout.write('{}\t{}\n'.format(factor, np.log10(1e-10)))
            else:
                fout.write('{}\t{}\n'.format(factor, np.mean(factorsExpression[factor])))

    foutn=os.path.join(outdir, "factorExpression.txt")    
    scores_df = pd.read_table(os.path.join(outdir, "factorExpression0.txt"), sep="\t",index_col=0)
    #scores_df['factorExpressionRank'] = preprocessing.MinMaxScaler().fit_transform(rankdata(scores_df['factorExpression'], method='average'))
    scores_df['factorExpressionRank'] = preprocessing.MinMaxScaler().fit_transform(rankdata(scores_df['factorExpression'], method='average').reshape(-1,1))
    scores_df.to_csv(foutn, sep='\t')

def get_correlation(corrfiles, features, outdir):
    df = pd.read_hdf(features)
    df=df[["source_target"]]
    df.source_target=[i.upper() for i in list(df.source_target)]
    df = df.set_index("source_target")

    for i, corrfile in enumerate(corrfiles):
        corr = pd.read_table(corrfile, sep="\t", index_col=0)
        corr = corr.rename(columns={corr.columns[0]:"corr_file{}".format(i + 1)})
        df = df.join(corr)
    
    outfile = os.path.join(outdir, "correlation.txt")
    df.to_csv(outfile, sep="\t")

def get_peakRPKM(fin_rpkm, outdir):
    peaks = pd.read_table(fin_rpkm, names=["chrom", "start", "end", "peakRPKM"])
    peaks["peak"] = peaks["chrom"] + ":" + peaks["start"].astype(str) +  "-" + peaks["end"].astype(str)
    add = peaks["peakRPKM"][peaks["peakRPKM"] > 0].min()
    peaks["log10_peakRPKM"] = np.log10(peaks["peakRPKM"] + add)
    peaks["peakRPKMScale"] = minmax_scale(peaks["log10_peakRPKM"])
    peaks["peakRPKMRank"] = minmax_scale(rankdata(peaks["log10_peakRPKM"]))
    
    cols = ["peak", "peakRPKM", "log10_peakRPKM", "peakRPKMScale", "peakRPKMRank"]
    outname = os.path.join(outdir, "peakRPKM.txt")
    peaks[cols].to_csv(outname, sep="\t", index=False)

def get_binding_score(pwm_weight, peak_weight, pwmfile, factortable, model, outdir):
    
    # Load model
    with open(model, "rb") as f:
        clf = pickle.load(f)
    pwm = dd.read_csv(pwm_weight, sep="\t")
    peak = dd.read_csv(peak_weight, sep="\t")
    ft = dd.read_csv(factortable, sep="\t") 

    r = pwm.merge(peak, left_on="enhancer", right_on="peak")[["motif", "enhancer",  "zscore", "peakRPKMScale"]]
    # print(r.head())
    # print(ft.head())
    r = r.merge(ft, left_on="motif", right_on="motif")
    r = r.groupby(["factor", "enhancer"])[["zscore", "peakRPKMScale"]].max()
    r = r.dropna().reset_index()

    cache = Chest(available_memory=20e9)
    print("combining tables")
    table = r.compute()
    print("predict")
    table["binding"] = clf.predict_proba(table[["zscore", "peakRPKMScale"]])[:,1]
    print("save results")
    outfile = os.path.join(outdir, "binding.predicted.h5")
    table.to_hdf(outfile, key="/binding", format="table")
    table.to_csv(outfile.replace("h5","txt"), sep="\t", index=False)

def calculate_binding(fin_rpkm, gene_bed, outdir, genome="hg19", pwmfile=None, fpomoter=True, detail=True):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    config = MotifConfig()
    if pwmfile is None:
        params = config.get_default_params() 
        pwmfile = os.path.join(config.get_motif_dir(), params["motif_db"])
    
    motifs2factors = pwmfile.replace(".pfm", ".motif2factors.txt")
    factortable= pwmfile.replace(".pfm", ".factortable.txt")

    pwm_weight = os.path.join(outdir, "pwmScore.txt")
    peak_weight = os.path.join(outdir, "peakRPKM.txt")
    binding = os.path.join(outdir, "binding.predicted.h5")
    nfin_rpkm=os.path.join(outdir, "filtered_enahncers.txt")
    if not os.path.exists(nfin_rpkm):
        filter_promoter_peaks(gene_bed, fin_rpkm, fpomoter=fpomoter)

    if not os.path.exists(pwm_weight):
        get_PWMScore(nfin_rpkm, pwmfile, outdir, genome=genome)

    if not os.path.exists(peak_weight):
        get_peakRPKM(nfin_rpkm, outdir)
    
    pbar = ProgressBar()
    pbar.register()
    if not os.path.exists(binding):
        model = "/home/qxu/git/network2/db/dream_model.txt"
        get_binding_score(pwm_weight, peak_weight, pwmfile, factortable, model, outdir)

if __name__ == "__main__":
    description = ""
    usage = "%(prog)s [-h] [options]"
    parser = argparse.ArgumentParser(usage=usage,
                                    description=description,
                                    formatter_class=argparse.RawDescriptionHelpFormatter
                                    )

    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "-p", "--pwmfile",
        dest="pwmfile",
        help="PWM",
        metavar="FILE",
        default=None
    )
    
    parser.add_argument(
        "-r",
        required=True,
        dest="fin_rpkm",
        help="BED file with RPKM on the 4th column",
        metavar="FILE"
    )
    
    parser.add_argument(
        "-a",
        dest="annotation",
        help="Gene annotation in BED12 format",
        metavar="BED",
    )
    
    parser.add_argument(
        "-g",
        dest="genome",
        help="Genome",
        metavar="NAME",
        default="hg19",
    )
    
    parser.add_argument(
        "-o",
        required=True,
        dest="outdir",
        help="Output directory",
        metavar="DIR",
        default=None
    )

    parser.add_argument(
        "-f",
        dest="fpomoter",
        help="Filter promoters, True or False, input should be either 'True' or 'False'.",
        metavar="NAME",
        type=ast.literal_eval,
        default=True,
    )

    parser.add_argument(
        "-d",
        dest="detail",
        help="Keep detail files, True or False, input should be either 'True' or 'False'.",
        metavar="NAME",
        type=ast.literal_eval,
        default=True,
    )

    args = parser.parse_args()
    pwmfile = args.pwmfile
    fin_rpkm = args.fin_rpkm
    outdir = args.outdir
    genome = args.genome
    gene_bed = args.annotation
    fpomoter = args.fpomoter
    detail = args.detail

    calculate_binding(
            fin_rpkm,
            gene_bed, 
            outdir,
            genome=genome,
            pwmfile=pwmfile,
            fpomoter=fpomoter,
            detail=detail
            )           
