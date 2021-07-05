import os
from collections import namedtuple
from tempfile import NamedTemporaryFile

import numpy as np
import pytest
import pandas as pd

from ananse.network import Network
from ananse.commands import network
from .test_02_utils import write_file


@pytest.fixture
def binding_fname():
    return "tests/example_data/binding2.tsv"


@pytest.fixture
def network_obj():
    genome = "tests/data/genome.fa"
    if not os.path.exists(genome):
        write_file(genome, [">chr1", "N"])

    return Network(genome=genome, gene_bed="ananse/db/hg38.genes.bed")


def test_unique_enhancer(network_obj, binding_fname):
    regions = network_obj.unique_enhancers(binding_fname)
    regions = regions.as_df()
    assert regions.shape[0] == 6
    assert sorted(list(regions["Chromosome"].unique())) == ["chr1", "chr10", "chr17"]
    assert sorted(list(regions["Start"].unique())) == [7677184, 7687827]


def test_distance_weight(network_obj):
    dw = network_obj.distance_weight(
        include_promoter=True,
        promoter_region=20,
        full_weight_region=50,
        maximum_distance=100,
        alpha=5,
    )

    assert list(dw.columns) == ["weight", "dist"]
    dw = dw.set_index("dist")

    assert dw.loc[0, "weight"] == 1
    assert dw.loc[25, "weight"] == 1
    assert dw.loc[50, "weight"] == 1
    assert dw.loc[51, "weight"] < 1
    assert np.isclose(dw.loc[100, "weight"], 0, atol=1e-4)
    assert dw.shape[0] == 101

    dw = network_obj.distance_weight(
        include_promoter=False,
        promoter_region=20,
        full_weight_region=50,
        maximum_distance=100,
        alpha=5,
    )

    assert list(dw.columns) == ["weight", "dist"]
    dw = dw.set_index("dist")
    assert dw.loc[0, "weight"] == 0
    assert dw.loc[20, "weight"] == 0
    assert dw.loc[21, "weight"] == 1
    assert dw.shape[0] == 101


def test_command():
    with NamedTemporaryFile() as tmp:
        fname = tmp.name
        Args = namedtuple(
            "args",
            "genome annotation include_promoter include_enhancer binding fin_expression outfile ncore",
        )
        args = Args(
            genome="hg38",
            annotation=None,
            include_promoter=True,
            include_enhancer=True,
            binding="tests/data/network/binding.tsv.gz",
            fin_expression="tests/data/network/heart_expression.txt",
            outfile=fname,
            ncore=2,
        )
        network(args)

        df = pd.read_table(fname)
        assert df.shape[0] == 30690
        assert df.shape[1] == 2
        assert df.columns[0] == "tf_target"
        assert df.columns[1] == "prob"
