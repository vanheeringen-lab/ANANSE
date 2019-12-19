#!/usr/bin/env python
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
import warnings

warnings.filterwarnings("ignore")


class Network(object):
    def __init__(self):
        pass

    def create_network(self, featurefile, outfile, impute=False):

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
            #                 "factor_expression.rank.scale", "target_expression.rank.scale",
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

        bpd.to_csv(outfile, sep="\t")

    def run_network(self, featurefile, outfile):
        self.create_network(featurefile, outfile)
