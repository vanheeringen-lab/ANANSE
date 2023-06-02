import os
import re
import sys
import urllib.request
from glob import glob

import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from sklearn.metrics import average_precision_score, roc_auc_score

import ananse

logger.remove()
logger.add(
    sys.stderr, format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | {level} | {message}"
)

TFP_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=TF_Perturbations_Followed_by_Expression"
TRRUST_URL = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
MSIGDB_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/c3.all.v7.4.symbols.gmt"


def download_trrust_reference(outfile):
    edges = []
    with urllib.request.urlopen(
        TRRUST_URL,
    ) as f:
        for line in f.readlines():
            tf, target, regtype, pmid = line.decode().strip().split("\t")
            # Just skip repression for now
            if regtype in ["Activation", "Unknown"]:
                edges.append([tf, target, 1])
    edges = pd.DataFrame(edges, columns=["tf", "target", "interaction"])
    edges.to_csv(outfile, sep="\t", index=False)


def download_msigdb_reference(outfile):
    with urllib.request.urlopen(MSIGDB_URL) as gmt, open(outfile, "w") as fl1:
        for line in gmt:
            a = line.decode("utf-8").split()
            tf = a[0].split("_")[0]
            targets = a[2:]
            for target in targets:
                fl1.write(f"{tf}\t{target}\n")


def fix_columns(df):
    """Make sure network has a tf and a target column."""

    df.columns = df.columns.str.lower()
    df = df.rename(
        columns={
            "source": "tf",
            "source_target": "tf_target",
            "target_gene": "target",
        }
    )

    if "tf_target" in df.columns:
        for sep in [ananse.SEPARATOR, "_", "-"]:
            if df["tf_target"].head(100).str.contains(sep).sum() == 100:
                break

        df[["tf", "target"]] = df["tf_target"].str.split(sep, expand=True).iloc[:, :2]
        df = df.drop(columns=["tf_target"])

    if "tf" not in df.columns:
        raise ValueError("Expect a column named 'source' or 'tf'")
    if "target" not in df.columns:
        raise ValueError("Expect a column named 'target' or 'target_gene'")

    return df


def prepare_reference_network(network, filter_tfs=True):
    """Generate reference network.

    This network contains all possible edges, based on the TFs
    and the target genes in the input. TFs are optionally filtered
    to contain only validated TFs.

    Returns
    -------
        DataFrame with column `"interaction"` having 1 for a validated
        edge and 0 otherwise.
    """
    if isinstance(network, pd.DataFrame):
        df = network.reset_index()
    elif isinstance(network, str):
        if network.endswith("feather"):
            df = pd.read_feather(network)
        else:
            df = pd.read_table(network)
    else:
        raise ValueError("Unknown network type, need DataFrame or filename.")

    df = fix_columns(df)

    interaction_column = None
    for col in df.columns:
        if col in ["tf", "target"]:
            continue
        vals = df[col].unique()

        if len(vals) in [1, 2] and 1 in vals:
            interaction_column = col
            break

    tfs = set(df["tf"].unique())
    if filter_tfs:
        valid_tfs = set(get_tfs())
        tfs = list(tfs.intersection(valid_tfs))

    targets = df["target"].unique()

    # logger.info(
    #     f"{os.path.split(network)[-1]} reference - {len(tfs)} TFs, {len(targets)} targets, {df.shape[0]} edges."
    # )

    total = []
    for tf in tfs:
        for target in targets:
            total.append([tf, target])
    total = pd.DataFrame(total, columns=["tf", "target"]).set_index(["tf", "target"])

    if interaction_column is not None:
        logger.info(f"Using '{interaction_column}' as interaction column.")
        df = df.set_index(["tf", "target"])[[interaction_column]].rename(
            columns={interaction_column: "interaction"}
        )
    else:
        logger.info("No column with 1 found, assuming all lines are positive edges.")
        df = df.set_index(["tf", "target"])
        df["interaction"] = 1

    return total.join(df[["interaction"]]).fillna(0)


def _read_dorothea_reference(fname):
    dorothea = pd.read_table(fname)
    cols = [
        "is_evidence_chip_seq",
        "is_evidence_curated",
        "is_evidence_inferred",
        "is_evidence_tfbs",
    ]
    dorothea = dorothea.set_index(["tf", "target"])[cols]
    for col in cols:
        dorothea[col] = dorothea[col].astype(int)
    dorothea["dorothea"] = np.any(dorothea[cols] == 1, 1).astype(int)

    dorothea = dorothea.reset_index()

    tfs = set(dorothea["tf"].unique())
    valid_tfs = set(get_tfs())
    tfs = list(tfs.intersection(valid_tfs))
    targets = dorothea["target"].unique()
    logger.info(
        f"Dorothea reference - {len(tfs)} TFs, {len(targets)} targets, {dorothea.shape[0]} edges."
    )

    total = []
    for tf in tfs:
        for target in targets:
            total.append([tf, target])
    total = pd.DataFrame(total, columns=["tf", "target"]).set_index(["tf", "target"])

    dorothea = dorothea.set_index(["tf", "target"])
    dorothea = total.join(dorothea)
    dorothea = dorothea.fillna(0)
    return dorothea


def _read_enrichr_perturbation_reference(fname=None):
    """Uses the TF perturbations from Enrichr[1,2] to create reference edges.

    Targets are defined by up- or down-regulated gened from the following sets:
    Up: INDUCTION, ACTIVATION, OE.
    Down: KD, KO, INACTIVATION, DEPLETION, SIRNA, SHRNA, KNOCKOUT, DELETION INHIBITION.

    The TF and targets in the DataFrame consists of the Cartesian product of
    all TFs and target genes that occur in the set.

    Returns
    -------
        DataFrame with tf-target edges.

    References
    ----------
    .. [1] Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NfR, Ma'ayan A.
       "Enrichr: interactive and collaborative HTML5 gene list enrichment analysis
       tool." BMC Bioinformatics. 2013;128(14)

    .. [2] Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z,
       Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD,
       Gundersen GW, Ma'ayan A. "Enrichr: a comprehensive gene set enrichment
       analysis web server 2016 update." Nucleic Acids Research. 2016; gkw377.
    """
    use_online = False
    if fname:
        fopen = open(fname)
    else:
        logger.info(
            "No filename provided for TF perturbations, downloading from Enrichr"
        )
        fopen = urllib.request.urlopen(TFP_URL)
        use_online = True
    p = re.compile(r"(\w+)\s+(\w+)\s+(.+)\s+(\w+)")
    all_info = []
    edges = []
    with fopen as f:
        for line in f:
            if use_online:
                line = line.decode("utf-8")
            vals = line.strip().split("\t")
            m = re.search(p, vals[0])
            all_info.append(m.groups(0))
            if (
                m.group(2) in ["INDUCTION", "ACTIVATION", "OE"] and m.group(4) == "UP"
            ) or (
                m.group(2)
                in [
                    "KD",
                    "KO",
                    "INACTIVATION",
                    "DEPLETION",
                    "SIRNA",
                    "SHRNA",
                    "KNOCKOUT",
                    "DELETION",
                    "INHIBITION",
                ]
                and m.group(4) == "DOWN"
            ):
                tf = m.group(1)
                for target in vals[2:]:
                    edges.append([tf, target])
    all_info = pd.DataFrame(all_info, columns=["tf", "exp", "info", "up_down"])

    perturb_df = pd.DataFrame(edges, columns=["tf", "target"])
    tfs = set(perturb_df["tf"].unique())
    targets = perturb_df["target"].unique()

    logger.info(
        f"TF perturbation reference - {len(tfs)} TFs, {len(targets)} targets, {perturb_df.shape[0]} edges."
    )

    perturb_df["experiments"] = 1
    perturb_df = perturb_df.groupby(["tf", "target"]).count()
    perturb_df["interaction"] = 1
    perturb_df.columns = ["perturb_experiments", "perturb_interaction"]

    valid_tfs = set(get_tfs())
    tfs = list(tfs.intersection(valid_tfs))

    total = []
    for tf in tfs:
        for target in targets:
            total.append([tf, target])
    total = pd.DataFrame(total, columns=["tf", "target"]).set_index(["tf", "target"])

    perturb_df = total.join(perturb_df).fillna(0)
    return perturb_df


def get_tfs():
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


def read_network(fname, name=None):
    network = fname
    if fname.endswith("feather"):
        df = pd.read_feather(network)
    else:
        df = pd.read_table(network)
    df = fix_columns(df)

    # Make sure there are no duplicate interactions
    df.drop_duplicates(subset=["tf", "target"], inplace=True)
    df = df.set_index(["tf", "target"])

    # Assuming last column is the edge weight
    df = df.iloc[:, [-1]]
    if name is not None:
        df.columns = [name]

    return df


def _read_correlation_reference(network, corCutoff=0.6):
    tfs_name = f"{os.path.dirname(ananse.__file__)}/db/tfs.txt"
    tfs = pd.read_csv(tfs_name, header=None)[0].tolist()
    edb = pd.read_csv(network, sep="\t")

    edb["iscorrelation"] = [1 if i > corCutoff else 0 for i in edb["correlationRank"]]

    edb[["tf", "target"]] = edb["source_target"].str.split("_", expand=True).iloc[:, :2]
    edb = edb.drop(
        columns=["source_target", "ocorrelation", "correlation", "correlationRank"]
    )
    edb = edb[edb.tf.isin(tfs)]
    edb = edb.set_index(["tf", "target"])
    return edb


def _read_goterm_reference(network, goCutoff=0):
    tfs_name = f"{os.path.dirname(ananse.__file__)}/db/tfs.txt"
    tfs = pd.read_csv(tfs_name, header=None)[0].tolist()
    gdb = pd.read_csv(network, sep="\t", header=None)

    gdb["isgo"] = [1 if i > goCutoff else 0 for i in gdb[2]]

    gdb = gdb.rename(columns={3: "tf", 1: "target"})
    gdb = gdb[gdb.tf.isin(tfs)]

    gdb = gdb.drop(columns=[0, 2])
    gdb = gdb.set_index(["tf", "target"])
    return gdb


def _read_msigdb_reference(network):
    msidb = pd.read_csv(network, sep="\t", header=None)
    msidb = msidb.rename(columns={0: "tf", 1: "target"})
    msidb = msidb.set_index(["tf", "target"])
    msidb["interaction"] = 1
    return msidb


def _read_regnet_reference(network):
    regnet = pd.read_csv(network)
    regnet = regnet.rename(
        columns={"regulator_symbol": "tf", "target_symbol": "target"}
    )
    regnet = regnet.set_index(["tf", "target"])
    regnet["interaction"] = 1
    return regnet[["interaction"]]


def read_reference(name, fname=None):
    """
    Valid reference networks (name):
    - dorothea
    - perturbation
    - correlation
    - goterm
    - msigdb
    - regnet
    - trrust
    """
    if name.lower() == "dorothea":
        return _read_dorothea_reference(fname)
    if name.lower() == "perturbation":
        return prepare_reference_network(_read_enrichr_perturbation_reference(fname))
    if name.lower() == "correlation":
        return prepare_reference_network(_read_correlation_reference(fname, 0.6))
    if name.lower() == "goterm":
        return prepare_reference_network(_read_goterm_reference(fname, 0))
    if name.lower() == "msigdb":
        return prepare_reference_network(_read_msigdb_reference(fname))
    if name.lower() == "regnet":
        return prepare_reference_network(_read_regnet_reference(fname))
    if name.lower() == "trrust":
        return prepare_reference_network(fname)


def validate_files(fnames, ignore_missing=False):
    file_error = False
    for fname in fnames:
        if not os.path.exists(fname):
            logger.error(f"file {fname} does not exist")
            file_error = True
    if not ignore_missing and file_error:
        raise ValueError("One or more files not found!")


def read_networks(network_dict, ignore_missing=False):
    """Read predicted networks.

    Input is a dictionary with name as key and filename as value.
    """
    # Validate files first
    validate_files(network_dict.values(), ignore_missing=ignore_missing)

    df = pd.DataFrame({"tf": [], "target": []}).set_index(["tf", "target"])

    for name, fname in network_dict.items():
        if os.path.exists(fname):
            logger.info(f"Reading {name}")
            tmp = read_network(fname, name=name)
            logger.info(f"Merging {name}")
            df = df.join(tmp, how="outer")

    return df


class NetworkBenchmark:
    def __init__(self, ref_dir, outdir):
        logger.info("Reading reference networks")
        self.dorothea = read_reference(
            "dorothea", fname=f"{ref_dir}/dorothea/entire_database.txt"
        )
        self.perturb = read_reference("perturbation")
        self.trrust = read_reference("trrust", fname=f"{ref_dir}/trrust.txt")

        self.references = [
            [self.dorothea, "dorothea", "is_evidence_curated"],
            [self.perturb, "TF perturbations", "interaction"],
            [self.trrust, "TRRUST", "interaction"],
        ]
        self.outdir = outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.benchmark_results = []

    def _merge_files(self):
        for network in self.networks:
            outfile = f"{self.outdir}/{network}.feather"

            if os.path.exists(outfile):
                continue

            logger.info(f"Creating merged file for {network} networks.")
            network_dfs = []
            fnames = glob(self.network_files[network])
            for fname in fnames:
                name = fname
                for part in self.network_files[network].split("*"):
                    name = name.replace(part, "")
                logger.info(fname)
                df = read_network(fname)
                df.columns = [f"{network}.{name}"]
                network_dfs.append(df)

            logger.info("Concatenating...")
            df = pd.concat(network_dfs, axis=1)
            df.reset_index().to_feather(outfile)
            logger.info("Done.")

    def benchmark(self, network_files):
        self.network_files = network_files
        self.networks = list(self.network_files.keys())

        # Create a feather file with all individual GRNs merged
        # Saves time if (part of) the benchmark needs to be rerun
        self._merge_files()

        for network in self.networks:  # networks:
            fname = f"{self.outdir}/{network}.feather"
            logger.info(f"reading {fname}")
            df = pd.read_feather(fname)
            logger.info("done")
            df.set_index(df.columns[:2].tolist(), inplace=True)
            cols = df.columns
            min_value = df.min().min() - 1

            for ref_network, name, ref_col in self.references:
                logger.info(f"joining {name}")
                test_network = ref_network.join(df)
                logger.info("done")
                test_network = test_network.fillna(min_value)

                for col in cols:
                    pr_auc = average_precision_score(
                        test_network[ref_col], test_network[col]
                    )
                    roc_auc = roc_auc_score(test_network[ref_col], test_network[col])
                    exp_name = os.path.splitext(os.path.basename(col))[0]
                    logger.info(
                        f"{name}\t{network}\t{exp_name}\t{pr_auc:0.5f}\t{roc_auc:0.2f}"
                    )
                    self.benchmark_results.append(
                        [name, network, exp_name, pr_auc, roc_auc]
                    )

    def plot(self, outname=None):
        bench_df = pd.DataFrame(
            [row[:2] + [row[2][-1]] + row[3:] for row in self.benchmark_results],
            columns=["reference", "network", "tissue", "pr_auc", "roc_auc"],
        )

        for metric in ["pr_auc", "roc_auc"]:
            sns.set_context(
                "paper",
                rc={
                    "font.size": 16,
                    "axes.titlesize": 16,
                    "axes.labelsize": 16,
                    "xtick.labelsize": 16,
                    "ytick.labelsize": 16,
                },
            )
            order = (
                bench_df[bench_df["network"] != "baseline"]
                .groupby("network")
                .mean()[metric]
                .sort_values()
                .index[::-1]
            )

            g = sns.catplot(
                data=bench_df[bench_df["network"] != "baseline"],
                y="network",
                x=metric,
                col="reference",
                kind="box",
                order=order,
                fliersize=0,
                boxprops={"alpha": 0.5},
                sharex=False,
                aspect=1.3,
            )

            for subplot in g.axes:
                for ax in subplot:
                    for line in ax.lines:
                        line.set_alpha(0.3)

            g.map(
                sns.stripplot,
                metric,
                "network",
                color="black",
                # data=bench_df[bench_df["network"] != "baseline"],
                s=8,
                jitter=0.25,
                alpha=0.7,
                order=order,
            )
            # g.set_xticklabels(rotation=30)

            if outname:
                extension = os.path.splitext(outname)[1]
                g.savefig(outname.replace(extension, f".{metric}{extension}"))
