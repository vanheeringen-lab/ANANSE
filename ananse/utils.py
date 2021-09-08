import atexit
import getpass
import os
import pwd
import shutil
import subprocess as sp
import tempfile
import warnings

import genomepy.utils
from pybedtools import BedTool
import pysam
import pandas as pd


def check_path(arg, error_missing=True):
    """Expand all paths. Can check for existence."""
    if arg is None:
        return arg

    args = [arg] if isinstance(arg, str) else arg
    paths = [cleanpath(arg) for arg in args]
    if error_missing:
        for path in paths:
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"'{os.path.basename(path)}' not found in '{os.path.dirname(path)}'."
                )
    return paths[0] if isinstance(arg, str) else paths


def cleanpath(path):
    """Expand any path input to a literal path output"""
    return os.path.abspath(  # expand relative paths (inc './' and '../')
        os.path.expanduser(  # expand '~'
            os.path.expandvars(path)  # expand '$VARIABLES'
        )
    )


def shhh_bedtool(func):
    """
    Decorator that silences pybedtools RuntimeWarnings such as
    `line buffering (buffering=1) isn't supported in binary mode`
    """

    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            func(*args, **kwargs)

    return wrapper


@shhh_bedtool
def bed_sort(bed):
    """
    Sort a bed file.
    """
    tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
    try:
        tmpfile = os.path.join(tmpdir, os.path.basename(bed))
        BedTool(bed).sort(output=tmpfile)
        shutil.copy2(tmpfile, bed)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@shhh_bedtool
def bed_merge(list_of_beds, merged_bed):
    """
    merge any number of bed files (merges overlapping regions)
    """
    bed = BedTool(list_of_beds[0])
    if list_of_beds[1:]:
        bed = bed.cat(*list_of_beds[1:])
    bed.saveas(merged_bed)


@shhh_bedtool
def count_reads(bams, peakfile, bed_output):
    """
    Count bam reads in putative enhancer regions
    """
    # replace with gimmemotifs.preprocessing.coverage_table()
    bed = BedTool(peakfile)
    bam_list = bams if isinstance(bams, list) else [bams]
    bed.multi_bam_coverage(bams=bam_list, output=bed_output)


def samc(ncore):
    """set decent samtools range for samtools functions (1-5 total threads)"""
    return max(0, min(ncore - 1, 4))


def bam_index(bam, force=True, ncore=1):
    if force or not os.path.exists(f"{bam}.bai"):
        index_parameters = [f"-@ {samc(ncore)}", bam]
        pysam.index(*index_parameters)  # noqa


def bam_sort(bam, ncore=1):
    tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
    try:
        sorted_bam = os.path.join(tmpdir, os.path.basename(bam))
        sort_parameters = [f"-@ {samc(ncore)}", "-o", sorted_bam, bam]
        pysam.sort(*sort_parameters)  # noqa: pysam bug

        shutil.copy2(sorted_bam, bam)
        bam_index(bam, force=True, ncore=ncore)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def bam_merge(list_of_bams, merged_bam, ncore=1):
    """
    merge any number of (sorted) bam files
    """
    [bam_index(bam, force=False, ncore=ncore) for bam in list_of_bams]
    if len(list_of_bams) > 1:
        merge_parameters = ["-f", f"-@ {samc(ncore)}", merged_bam] + list_of_bams
        pysam.merge(*merge_parameters)  # noqa: pysam bug
        bam_index(merged_bam)
    else:
        # os.symlink() doesn't work with multi_bam_coverage()
        bam = list_of_bams[0]
        shutil.copy2(bam, merged_bam)
        shutil.copy2(f"{bam}.bai", f"{merged_bam}.bai")


def mosdepth(bed, bam, bed_output, ncore=1):
    """
    Count (median per base overlap of) bam reads in putative enhancer regions
    """
    ncore = min(4, ncore)
    tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
    try:
        prefix = os.path.join(tmpdir, "bam_coverage")
        cmd = f"mosdepth -nxm -t {ncore} -b {bed} {prefix} {bam}"
        sp.check_call(cmd, shell=True)

        tmp_bed_output = f"{prefix}.regions.bed"
        cmd = f"gunzip -f {tmp_bed_output}.gz"
        sp.check_call(cmd, shell=True)

        shutil.copy2(tmp_bed_output, bed_output)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# def bed_sum_coverages(multi_bam_coverage, sum_bam_coverage):
#     """
#     MultiBamCov returns a BED3+n with one column per bam.
#     This function sums up all bam counts and returns a BED3+1.
#     """
#     bed = pd.read_csv(multi_bam_coverage, header=None, sep="\t")
#     columns = bed.shape[1]
#     if columns != 4:
#         bed3 = bed.iloc[:, :3]
#         scores = bed.iloc[:, 3:].sum(axis=1)
#         bed = pd.concat([bed3, scores], axis=1)
#     bed.to_csv(sum_bam_coverage, sep="\t", header=False, index=False)


# def non_empty_files(files, error_msg, size_threshold=10, verbose=True):
#     """Check list of files for content
#
#     Args:
#         files: list of filepaths
#         error_msg: message for all empty files
#         size_threshold: minimum size to be considered non-empty
#         verbose: return warnings?
#
#     Returns:
#         list of non-empty files
#     """
#     ret_files = []
#     for file in files:
#         if os.path.getsize(file) > size_threshold:
#             ret_files.append(file)
#         elif verbose:
#             logger.warning(f"Empty file: '{os.path.basename(file)}'")
#
#     if not ret_files:
#         logger.exception(error_msg)
#         exit(1)
#     return ret_files


def mytmpdir():
    """
    returns a temp directory that is removed when the process is completed
    the directory is not removed if the process is killed by the user
    """
    if not hasattr(mytmpdir, "dir") or not mytmpdir.dir:
        # can also be removed with clean_tmp()
        mytmpdir.dir = tempfile.mkdtemp(prefix=f"ANANSE_{os.getpid()}.")
        atexit.register(shutil.rmtree, mytmpdir.dir, ignore_errors=True)
    return mytmpdir.dir


def clean_tmp():
    """
    remove leftover temp dirs
    temp dirs are left if ANANSE was killed by the user
    """
    user = getpass.getuser()
    tempdir = tempfile.gettempdir()

    # all tmp files/directories starting with "ANANSE_" & owner by the user
    tmp_files = os.listdir(tempdir)
    ananse_files = [
        os.path.join(tempdir, f) for f in tmp_files if f.startswith("ANANSE_")
    ]
    user_files = [
        f
        for f in ananse_files
        if os.path.exists(f) and pwd.getpwuid(os.stat(f).st_uid).pw_name == user
    ]

    # delete
    _ = [genomepy.utils.rm_rf(f) for f in user_files]


def get_motif_factors(motif, indirect=True):
    """Return all TFs that are associated with a motif."""
    motif_factors = []
    for factor_type, factors in motif.factors.items():
        if factor_type == "direct" or indirect:
            motif_factors += factors
    return motif_factors


def check_input_factors(factors):
    """Check factors.

    Factors can either be a list of transcription factors, or a filename of a
    file that contains TFs. Returns a list of factors.

    Returns
    -------
    list
        List of TF names.
    """
    # Load factors
    if factors is None:
        return

    # if factors is a string, assume it's a filename
    if isinstance(factors, str):
        fname = factors

    # if factors is a list of 1, and it exists, assume it's a filename
    elif isinstance(factors, list) and len(factors) == 1:
        fname = factors[0]

    # It's a list with more than one value, assuming it's a list of TF names.
    else:
        return factors

    if not os.path.exists(fname):
        raise ValueError(f"Factors file '{factors}' does not exist")

    factors = [line.strip() for line in open(fname)]
    return factors


def view_h5(fname, tfs=None, fmt="wide"):
    """Extract information from an ANANSE binding.h5 file.

    Parameters
    ----------
    fname : str
        File name (binding.h5).

    tfs : list, optional
        List of transcription factor names to extract. All TFs are used
        by default.

    fmt : str, optional
        Return output in 'wide' or in 'long' format. Default is 'wide'.

    Returns
    -------
    pandas.DataFrame
    """
    if fmt not in ["wide", "long"]:
        raise ValueError("fmt should be either 'wide' or 'long'")

    with pd.HDFStore(fname) as hdf:
        if tfs is None:
            tfs = [x for x in dir(hdf.root) if not x.startswith("_")]

        idx = hdf.get("_index")

        df = pd.DataFrame(index=idx.index)
        for tf in tfs:
            df[tf] = hdf.get(tf).values

    if fmt == "long":
        df.index.rename("loc", inplace=True)
        df = df.reset_index().melt(
            id_vars=["loc"], value_name="prob", var_name="factor"
        )

    return df
