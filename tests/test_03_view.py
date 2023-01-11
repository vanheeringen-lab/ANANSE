import os
from collections import namedtuple

import pandas as pd

import ananse.view
from ananse.commands import view


def test_get_binding_tfs():
    cmd = ananse.view.get_binding_tfs
    binding = "tests/data/network/binding.h5"
    ret = cmd(binding)
    assert len(ret) == 33
    assert "NFKB2" in ret


def test_view_h5():
    cmd = ananse.view.view_h5
    binding = "tests/data/network/binding.h5"

    # regions and TFs
    regions = cmd(binding, list_regions=True)
    assert regions.shape == (61786, 1)
    tfs = cmd(binding, list_tfs=True)
    assert tfs.shape == (33, 1)

    # all = regions x TFs
    complete_view = cmd(binding)
    assert complete_view.shape == (regions.shape[0], tfs.shape[0])

    # various subsets in either format
    head_wide = cmd(binding, n=10)
    assert head_wide.shape == (10, 10)
    head_long = cmd(binding, n=10, fmt="long")
    assert head_long.shape == (10 * 10, 3)

    subset_wide = cmd(
        binding,
        tfs=["NFKB2", "KLF6"],
        regions=[
            "chr10:100000171-100000371",
            "chr10:100001843-100002043",
        ],
    )
    assert subset_wide.shape == (2, 2)
    subset_long = cmd(
        binding,
        tfs=["NFKB2", "KLF6"],
        regions=[
            "chr10:100000171-100000371",
            "chr10:100001843-100002043",
        ],
        fmt="long",
    )
    assert subset_long.shape == (2 * 2, 3)


def test_command_view(outdir):
    Args = namedtuple(
        "args",
        "infile outfile tfs regions format n list_regions list_tfs activity",
    )
    outfile = os.path.join(outdir, "view.tsv")
    tfs = ["NFKB2", "KLF6"]
    regions = [
        "chr10:100000171-100000371",
        "chr10:100001843-100002043",
    ]
    args = Args(
        infile="tests/data/network/binding.h5",
        outfile=outfile,
        tfs=tfs,
        regions=regions,
        format="wide",
        n=None,
        list_regions=False,
        list_tfs=False,
        activity=False,
    )
    view(args)

    df = pd.read_csv(outfile, sep="\t").set_index("loc")
    assert set(df.columns) == set(tfs)
    assert set(df.index) == set(regions)
