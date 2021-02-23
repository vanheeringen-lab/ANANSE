import os

import numpy as np
import pytest

from ananse.network import Network
from .test_02_utils import write_file


@pytest.fixture
def binding_fname():
    return "tests/example_data/binding2.tsv"


@pytest.fixture
def network():
    genome = "tests/data/genome.fa"
    if not os.path.exists(genome):
        write_file(genome, [">chr1", "N"])

    return Network(genome=genome, gene_bed="ananse/db/hg38.genes.bed")


def test_unique_enhancer(network, binding_fname):
    regions = network.unique_enhancers(binding_fname)
    regions = regions.as_df()
    assert regions.shape[0] == 6
    assert sorted(list(regions["Chromosome"].unique())) == ["chr1", "chr10", "chr17"]
    assert sorted(list(regions["Start"].unique())) == [7677184, 7687827]


def test_distance_weight(network):
    dw = network.distance_weight(
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

    dw = network.distance_weight(
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
