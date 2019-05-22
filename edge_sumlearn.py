#!/usr/bin/env python
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, roc_curve
import numpy as np
from scipy.stats import rankdata
import dask.dataframe as dd
import sys
import os
import argparse
from sklearn.preprocessing import minmax_scale
import warnings
warnings.filterwarnings('ignore')


def create_network(featurefile, outdir, impute=False): 

    network = pd.read_hdf(featurefile, key="/features")
    
    exclude_cols = ["sum_weighted_logodds","enhancers", "log_enhancers", 
                    "sum_binding", 
                    "sum_logodds", 
                    "log_sum_binding", "factorExpression", "targetExpression", "factor", "gene", "factor_expression", 
                    "target_expression", 
                    "factor_expression.scale", "target_expression.scale",
        #                 "factor_expression.rank.scale", "target_expression.rank.scale",
                    "corr_file1", "correlation",
        #                 "correlationRank",
                        "max_binding_in_promoter",
                    "max_binding","max_sum_dist_weight",
        #                "sum_dist_weight"
                   ]
    network = network[[c for c in network.columns if c not in exclude_cols]]
    network = network.set_index("source_target")
    network['binding'] = minmax_scale(rankdata(network['sum_dist_weight'], method='dense'))
    network.drop(['sum_dist_weight'],axis=1,inplace=True)
    
    bp=network.mean(axis=1)
    bpd=pd.DataFrame(bp)
    bpd=bpd.rename(columns={0:"binding"})
    bpd['prob'] = minmax_scale(rankdata(bpd['binding'], method='dense'))
    
    bpd.to_csv(os.path.join(outdir, "SUMlearn_network.txt"), sep="\t")  
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--features",
        dest="features",
        help="HDF5 file with features",
        metavar="FILE",
        default=None,
        required=True,
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
        "-i", "--impute",
        dest="impute",
        help="Impute missing values",
        default=False,
        action="store_true",
    )

    args = parser.parse_args()

    featurefile = args.features
    outdir = args.outdir
    impute = args.impute

    create_network(
        featurefile,
        outdir,
        impute,
        )


# featurefile="../results/full_features.h5"
# outfile="../results/logistic_regression_network.txt"
# network = pd.read_hdf(featurefile, key="/features")


# ##python edge_sumlearn.py -f ../results/full_features.h5 \
#                                          -o ../results


