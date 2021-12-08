import atexit
import getpass
import os
import pwd
import shutil
import re
import pandas as pd
import tempfile
from typing import Union

import genomepy.utils
from ananse.bed import CombineBedFiles


def check_cores(ncore=None):
    if ncore is None:
        ncore = min(os.cpu_count(), 4)
    return int(ncore)


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


def parse_input(
    anything: Union[str, list] = None,
    none_ok=True,
    parse_list=lambda x: x,
    parse_files=None,
    *args,
    **kwargs,
) -> Union[None, list]:
    if anything is None:
        if not none_ok:
            raise ValueError("Input was None")
        return None

    if isinstance(anything, str):
        anything = [anything]
    elif not isinstance(anything, list):
        raise ValueError(
            f"Only None, string or list input accepted, received type {type(anything)}"
        )

    as_files = [cleanpath(a) for a in anything]
    if not all(os.path.exists(f) for f in as_files):
        # input is in list format
        out = parse_list(anything)
    else:
        # input is in file format
        out = parse_files(as_files, *args, **kwargs)
    return sorted(set(out))


def parse_tf_files(tf_files: list) -> list:
    """convert one or more files with one element per line into a list of elements"""
    out = []
    for infile in tf_files:
        for line in open(infile):
            out.append(line.strip().strip(","))
    return out


def test_tfs(tfs: list) -> list:
    # if it looks like it could be a list of files, check for file separators
    # (those should never be in TFs right?)
    if len(tfs) < 50:
        for tf in tfs:
            if os.sep in tf:
                raise ValueError(
                    f"Mixed input detected: some look like a file (e.g. {tf}), others do not."
                )
    return tfs


def merge_regionfiles(regionfiles, genome, outdir):
    cbed = CombineBedFiles(genome=genome, peakfiles=regionfiles)
    combined_bed = os.path.join(outdir, "regions_combined.bed")
    cbed.run(outfile=combined_bed, width=200, force=True)
    return combined_bed


def parse_regions(regionfiles, *args, **kwargs):
    # combine bed files into one
    infile = regionfiles[0]
    if len(regionfiles) > 1:
        infile = merge_regionfiles(regionfiles, *args, **kwargs)

    df = pd.read_table(infile, header=None, sep="\t", comment="#", dtype=str)
    # skip potential header line
    test = str(df.at[0, 0] if len(df) == 1 else df.at[1, 0])
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
    if not bool(re.match(r"^.*:\d+-\d+$", str(regions[0]))) and len(regions) > 1:
        regions = regions[1:]

    return regions


def test_regions(regions: list) -> list:
    if not bool(re.match(r"^.*:\d+-\d+$", str(regions[0]))):
        raise ValueError(
            f"A region list should be formatted like 'chr1:100-200'. Received {regions[0]}"
        )
    return regions


def load_tfs(tfs: Union[str, list] = None) -> list:
    """
    Unpacks a string or list of TFs or files with TFs
    """
    return parse_input(
        tfs,
        parse_list=test_tfs,
        parse_files=parse_tf_files,
    )


def load_regions(
    regions: Union[str, list] = None, genome: str = None, outdir: str = None
) -> list:
    """
    Unpacks a string or list of regions or files with regions
    """
    return parse_input(
        regions,
        parse_list=test_regions,
        parse_files=parse_regions,
        genome=genome,
        outdir=outdir,
    )
