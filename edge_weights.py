#!/usr/bin/env python
import argparse
import os
import pickle
import subprocess
import sys
import math

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

def filter_promoter_peaks(fin_regions_fa,outdir, genome="hg19"):
    pass

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

def filter_enhancer(gene_bed, peak_bed, up=2000, down=2000):
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
            if line.split()[0]+":"+line.split()[1]+"-"+line.split()[2] not in vals:
                fl2.write(line)

# 1 Motif PWM score
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
        # print (seqs)

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
                # pwm_score = pwm_score.join(pd.read_table(motif_stats, index_col=0))
                # # scale motif score
                # pwm_score["zscore"] = (pwm_score["score"] - pwm_score["mean"]) / pwm_score["std"] 
                pwm_score["zscoreRank"] = minmax_scale(rankdata(pwm_score["zscore"]))
                cols = ["enhancer", "zscore", "zscoreRank"]
                write_header = False
                if chunk == 0:
                    write_header = True
                pwm_score[cols].to_csv(fout, sep='\t',
                        header=write_header)
    else:
        print("Motif background distribution")
        motif_stats = os.path.join(outdir, "motif_distribution.txt")
        # get_motif_distribution(motif_stats, genome, pwmfile=pwmfile)
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
    
# 10 expression
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

# 10 Factor expression
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
    peaks = pd.read_table(fin_rpkm, 
            names=["chrom", "start", "end", "peakRPKM"])
    peaks["peak"] = peaks["chrom"] + ":" + peaks["start"].astype(str) +  "-" + peaks["end"].astype(str)
    add = peaks["peakRPKM"][peaks["peakRPKM"] > 0].min()
    peaks["log10_peakRPKM"] = np.log10(peaks["peakRPKM"] + add)
    peaks["peakRPKMScale"] = minmax_scale(peaks["log10_peakRPKM"])
    peaks["peakRPKMRank"] = minmax_scale(rankdata(peaks["log10_peakRPKM"]))
    
    cols = ["peak", "peakRPKM", "log10_peakRPKM", "peakRPKMScale", "peakRPKMRank"]

    outname = os.path.join(outdir, "peakRPKM.txt")
    peaks[cols].to_csv(outname, sep="\t", index=False)

#clff=sklearn.linear_model.LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
#              intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
#              penalty='l2', random_state=None, solver='liblinear', tol=0.0001,
#              verbose=0, warm_start=False)

def get_binding_score(pwm_weight, peak_weight, pwmfile, factortable, model, outdir):
    
    # Load model
    with open(model, "rb") as f:
        clf = pickle.load(f)

    # factor_table = os.path.join(outdir, "factor_table.txt")
    # with open(pwmfile) as f:
    #     motifs = read_motifs(f)
    # with open(factor_table, "w") as f:
    #     f.write("motif\tfactor\n")
    #     for motif in motifs:
    #         for factor in motif.factors:
    #             f.write("{}\t{}\n".format(motif.id, factor))
    
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
    table.to_csv(outfile.replace("h5","txt"), sep="\t")

def get_promoter_dataframe(gene_bed, peak_bed, genome, up=2000, down=2000):
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
        vals.append([chrom, gene, peak_start, peak_end])
    prom = pd.DataFrame(vals, columns=["chrom", "gene", "peak_start", "peak_end"])
    prom["loc"] = prom["chrom"] + ":" + prom["peak_start"].astype(str) + "-" + prom["peak_end"].astype(str)
    return prom
   
def get_gene_dataframe(gene_bed, peak_bed, genome, up=100000, down=100000):
    #all overlap Enh-TSS(100000-tss-100000) pair distance
    g = Genome(genome)
    gsize = g.props["sizes"]["sizes"]

    peaks = BedTool(peak_bed)
    b = BedTool(gene_bed)
    b = b.flank(l=1, r=0, s=True, g=gsize).slop(l=up, r=down, g=gsize, s=True)
    #bedtools flank  -r 0 -l 1 -i b.bed -g 
    ##all gene upstream 1bp position (TSS), Chr01 12800   12801   in Chr01    4170    12800   Xetrov90000001m.g   0   -   
    
    #bedtools slop  -r down -l up -i b.bed -g 
    ## |100000--TSS--100000|

    vals = []
    for f in b.intersect(peaks, wo=True, nonamecheck=True):
        #bedtools intersect -wo -nonamecheck -b peaks.bed -a b.bed
        ##
        chrom = f[0]
        strand = f[5]
        if strand == "+":
            tss = f.start + up
        else:
            tss = f.start + down
        gene = f[3]
        peak_start, peak_end = int(f[13]), int(f[14])
        vals.append([chrom, tss, gene, peak_start, peak_end])
    p = pd.DataFrame(vals, columns=["chrom", "tss", "gene", "peak_start", "peak_end"])
    p["peak"] = [int(i) for i in (p["peak_start"] + p["peak_end"]) / 2]
    #peak with int function, let distance int
    p["dist"] = np.abs(p["tss"] - p["peak"])
    p["loc"] = p["chrom"] + ":" + p["peak_start"].astype(str) + "-" + p["peak_end"].astype(str)
    p = p.sort_values("dist").drop_duplicates(["loc", "gene"], keep="first")[["gene", "loc", "dist"]] 
    return p

# def distance_weight(binding,distance):
#     if distance<=5000:
#         wbinding=np.log(binding+1)
#     elif 5000<distance<100000:
#         dist5k=distance-5000
#         wbinding=np.log(1/(dist5k/1000+1)*binding+1)
#     else:
#         wbinding=0
#     return(wbinding)

def distance_weight(binding,distance,alpha=1e5,padding=int(2e5)):
    print (binding,distance)

    # alpha is the effect of distance on regulatory potential. Distance at which regulatory potential is 1/2, (default=10kb)' 
    # u = -math.log(1.0/3.0)*1e5/alpha
    # weight  = np.array( [ 2.0*math.exp(-u*math.fabs(z)/1e5)/(1.0+math.exp(-u*math.fabs(z)/1e5))  for z in range( -padding,padding+1) ] )
    # wbinding=binding*weight[distance]

    # u = -math.log(1.0/3.0)*1e5/alpha
    # weight  = np.array( [ 2.0*math.exp(-u*math.fabs(z)/1e5)/(1.0+math.exp(-u*math.fabs(z)/1e5))  for z in range( 0,padding+1) ] )
    # wbinding=binding*weight[distance]
    # print (wbinding)
    # return(wbinding)

    u = -math.log(1.0/3.0)*1e5/alpha
    weight  = np.array( [ 2.0*math.exp(-u*math.fabs(z)/1e5)/(1.0+math.exp(-u*math.fabs(z)/1e5))  for z in range( 1,padding+1) ] )
    wbinding= np.dot( binding, weight[distance])
    return(wbinding)

# distance_weight([1,2,3],[100,200,1])

def aggregate_binding(binding, gene_bed, peak_bed, genome, outdir, window_up=100000, window_down=100000, alpha=1e4, padding=int(1e5), keep1=5000, remove=2000):
    # Overlaps
    g = Genome(genome)
    gsize = g.props["sizes"]["sizes"]
    
    print("promoter_overlap")
    prom = get_promoter_dataframe(gene_bed, peak_bed, genome)
    prom.gene=[i.upper() for i in list(prom.gene)]

    print("gene overlap")
    p = get_gene_dataframe(gene_bed, peak_bed, genome, window_up, window_down)
    p=p[p["dist"]<99999]
    # remove distance more than 100k interaction, for weight calculate
    p.gene=[i.upper() for i in list(p.gene)]

    ddf = dd.read_hdf(binding, key="/binding")[["factor", "enhancer", "binding"]]

    prom_table = ddf.merge(prom, left_on="enhancer", right_on="loc")
    prom_table = prom_table.groupby(["factor", "gene"])[["binding"]].max()
    prom_table = prom_table.rename(columns={"binding":"max_binding_in_promoter"})
    prom_table = prom_table.reset_index()
    prom_table["source_target"] = prom_table["factor"].map(str) + "_" + prom_table["gene"].map(str)

    f_table = ddf.merge(p, left_on="enhancer", right_on="loc")
    sum_enh = f_table.groupby(["factor", "gene"])[["binding"]].count()
    f_table["sum_weighted_logodds"] = f_table["binding"].div(f_table["binding"].mean()).apply(np.log, meta=('binding', np.float64)).rmul(50000).div(f_table["dist"])
    f_table["sum_logodds"] = f_table["binding"].div(f_table["binding"].mean()).apply(np.log, meta=('binding', np.float64))
    
    #f_table["tmp"] = f_table["binding"].div(f_table["binding"].mean())
    #f_table['sun_dist_weight'] = f_table.apply(lambda row: distance_weight(row['binding'], row['dist']), axis=1)
    
    # f_table=f_table.compute()
    # sum_enh=sum_enh.compute()
    # prom_table=prom_table.compute()

    u = -math.log(1.0/3.0)*1e5/alpha
    weight1  = pd.DataFrame({"weight": [0 for z in range(1, remove+1)], "dist" : range(1, remove+1)})
    weight2  = pd.DataFrame({"weight": [1 for z in range(remove+1,keep1+1)], "dist" : range(remove+1,keep1+1)})
    weight3  = pd.DataFrame({"weight": [2.0*math.exp(-u*math.fabs(z)/1e5)/(1.0+math.exp(-u*math.fabs(z)/1e5)) for z in range( 1,padding-keep1+1)], 
                            "dist" : range(keep1+1, padding+1)})

    weight = pd.concat([weight1, weight2, weight3])
    
    weight.to_csv(os.path.join(outdir, "weight.csv"))

    weight=dd.read_csv(os.path.join(outdir, "weight.csv"))
    f_table=f_table.merge(weight,how='left',on='dist')
    f_table['sum_dist_weight'] = f_table['binding'] * f_table['weight']

    # f_table['sun_dist_weight'] = distance_weight(f_table['binding'], f_table['dist'])

    #f_table = f_table.drop(['temp'], axis=1)

    f_table_sum = f_table.groupby(["factor", "gene"]).sum()[["sum_weighted_logodds", "sum_logodds", "binding", "sum_dist_weight"]]
    # f_table_max = f_table.groupby(["factor", "gene"]).max()[["binding", "sum_dist_weight"]]

    f_table_max = f_table.groupby(["factor", "gene"])[["binding","sum_dist_weight"]].max()

    f_table_sum = f_table_sum.rename(columns={"binding":"sum_binding"})
    f_table_max = f_table_max.rename(columns={"binding":"max_binding"})
    f_table_max = f_table_max.rename(columns={"sum_dist_weight":"max_sum_dist_weight"})

    sum_enh = sum_enh.rename(columns={"binding":"enhancers"})
    f_table_sum = f_table_sum.reset_index()
    f_table_max = f_table_max.reset_index()

    f_table = f_table.reset_index()
    sum_enh = sum_enh.reset_index()

    f_table_sum["source_target"] = f_table_sum["factor"] + "_" + f_table_sum["gene"]
    f_table_max["source_target"] = f_table_max["factor"] + "_" + f_table_max["gene"]

    f_table["source_target"] = f_table["factor"] + "_" + f_table["gene"]
    sum_enh["source_target"] = sum_enh["factor"] + "_" + sum_enh["gene"]
    f_table_max = f_table_max.rename(columns={"factor":"factor2"})    
    f_table_max = f_table_max.rename(columns={"gene":"gene2"})

    f_table = f_table_sum.merge(f_table_max, left_on="source_target", right_on="source_target", how="outer")
    f_table = f_table.merge(sum_enh, left_on="source_target", right_on="source_target", how="outer")
    f_table = f_table.merge(prom_table, left_on="source_target", right_on="source_target", how="outer")
    f_table = f_table[["source_target", "factor", "gene", "sum_weighted_logodds", "sum_dist_weight", "sum_logodds", "sum_binding", "enhancers", 
                       "max_binding_in_promoter","max_binding","max_sum_dist_weight"]]
    f_table["log_sum_binding"] = f_table["sum_binding"].add(1e-5).apply(np.log, meta=('sum_binding', np.float64))
    f_table["log_enhancers"] = f_table["enhancers"].add(1).apply(np.log, meta=('enhancers', np.float64))
    f_table["factor"] = f_table["source_target"].str.replace("_.*", '')
    f_table["gene"] = f_table["source_target"].str.replace(".*_", '')
    f_table["max_binding_in_promoter"] = f_table["max_binding_in_promoter"].fillna(0)

    outfile = os.path.join(outdir, "features.h5")
    print("computing, output file {}".format(outfile))
    f_table.to_hdf(outfile, key="/features")
    # f_table.to_csv(outfile.replace("h5","txt"), sep="\t", index=False)

def join_features(features, other, outfile):
    network = pd.read_hdf(features)

    for fname in other:
        df = pd.read_table(fname, sep="\t")
        for col in ["factor", "gene"]:
            if col in df.columns:
                df = df.drop(col, 1)
        network = network.merge(df, 
                left_on="source_target", 
                right_on="source_target",)
    
    # Compute before saving, will result in an error otherwise
    # network = network.compute()   
    network.to_csv(outfile.replace("h5","txt"), sep="\t", index=False)
    network.to_hdf(outfile, key="/features")
 
def calculate_features(fin_rpkm, gene_bed, outdir, genome="hg19", pwmfile=None, fin_expression=None, corrfiles=None):
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
        filter_enhancer(gene_bed, fin_rpkm)
    
    if not os.path.exists(pwm_weight):
        get_PWMScore(nfin_rpkm, pwmfile, outdir, genome=genome)
    

    if not os.path.exists(peak_weight):
        get_peakRPKM(nfin_rpkm, outdir)
    
    pbar = ProgressBar()
    pbar.register()
    
    if not os.path.exists(binding):
        model = "/home/qxu/projects/regulatoryNetwork/run_dream3/dream_model.txt"
        get_binding_score(pwm_weight, peak_weight, pwmfile, factortable, model, outdir)
    
    features = os.path.join(outdir, "features.h5")
    if not os.path.exists(features):
        aggregate_binding(binding, gene_bed, nfin_rpkm, genome, outdir, window_up=50000, window_down=50000)
     
    if fin_expression is not None:
        get_expression(fin_expression, features, outdir)

    if not os.path.exists(os.path.join(outdir, "factorExpression.txt")):
        get_factorExpression(fin_expression, motifs2factors, outdir)


    if not os.path.exists(os.path.join(outdir, "correlation.txt")):
        if corrfiles is not None: 
            get_correlation(corrfiles, features, outdir)

    other = [
        os.path.join(outdir, "expression.txt"),
        os.path.join(outdir, "correlation.txt"),
    ]
    outfile = os.path.join(outdir, 'full_features.h5')
    join_features(features, other, outfile) 
 
if __name__ == "__main__":
    description = ""
    usage = "%(prog)s [-h] [options]"
    parser = argparse.ArgumentParser(usage=usage,
    description=description,
    #epilog=epilog,
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
        "-e",
        dest="fin_expression",
        help="Expression scores",
        metavar="FILE",
        nargs='*'
    )
    
    parser.add_argument(
        "-c",
        dest="corrfiles",
        help="Files with correlation",
        metavar="FILE",
        nargs='*'
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
    
    args = parser.parse_args()
    pwmfile = args.pwmfile
    fin_expression = args.fin_expression
    fin_rpkm = args.fin_rpkm
    outdir = args.outdir
    genome = args.genome
    gene_bed = args.annotation
    corrfiles = args.corrfiles
    
    calculate_features(
            fin_rpkm,
            gene_bed, 
            outdir,
            genome=genome,
            pwmfile=pwmfile,
            fin_expression=fin_expression,
            corrfiles=corrfiles,
            )           
    # -r /home/qxu/projects/regulatoryNetwork/run20180806/data/p300/p300stage9.bed 
    # -e /home/qxu/projects/regulatoryNetwork/run20180806/data/RNAseq/ClutchA_polyA_4_5_hpf_chr1.tsv 
    # -o /home/qxu/projects/regulatoryNetwork/run20180806/results -a data/xtchr1n.bed 
    # -g xtchr1 
    # -c /home/qxu/projects/regulatoryNetwork/run20180806/results/expressioncorrelation.txt 
    # -p data/gimme.vertebrate.v3.3.1.xt.pwm
