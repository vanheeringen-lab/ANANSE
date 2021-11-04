from copy import deepcopy
import os
import pytest
from string import ascii_lowercase

import networkx as nx

import ananse.peakpredictor
from . import write_file


def test__check_input_regions(outdir, genome, bed1, bed2):
    regions = ananse.peakpredictor._check_input_regions(
        None, None, verbose=True, force=True
    )
    assert regions is None

    # >1 file (all BED)
    regionfiles = [bed1, bed2]
    regions = ananse.peakpredictor._check_input_regions(
        regionfiles, genome, outdir, verbose=True, force=True
    )
    assert regions == ["chr1:400-600", "chr1:2400-2600", "chr1:4400-4600"]

    # 1 file (BED)
    bed3 = os.path.join(outdir, "bed3.bed")
    write_file(bed3, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n", "chr1\t4000\t5000\n"])

    regions = ananse.peakpredictor._check_input_regions(
        [bed3], None, verbose=True, force=True
    )
    assert regions == ["chr1:0-1000", "chr1:2000-3000", "chr1:4000-5000"]

    # 1 file (regions /w header)
    bed4 = os.path.join(outdir, "bed4.bed")
    write_file(
        bed4, ["header", "chr1:0-1000\n", "chr1:2000-3000\n", "chr1:4000-5000\n"]
    )

    regions = ananse.peakpredictor._check_input_regions(
        [bed4], None, verbose=True, force=True
    )
    assert regions == ["chr1:0-1000", "chr1:2000-3000", "chr1:4000-5000"]


def test__check_input_files():
    missing_files = ["this_file_does_not_exist"]
    with pytest.raises(ReferenceError):
        try:
            _ = ananse.peakpredictor._check_input_files(missing_files)
        except SystemExit:
            raise ReferenceError  # arbitrary error

    present_files = [
        "tests/data/binding/hg38_testgenome.fa",
        "tests/data/binding/test.pfm",
    ]
    _ = ananse.peakpredictor._check_input_files(present_files)


def test__get_species():
    s = ananse.peakpredictor._get_species("hg38")
    assert s == "human"
    s = ananse.peakpredictor._get_species("mm9")
    assert s == "mouse"
    s = ananse.peakpredictor._get_species("GRCz11")
    assert s is None


def test__load_human_factors():
    valid_factors = ananse.peakpredictor._load_human_factors()
    # assert isinstance(valid_factors, list)
    assert "TP53" in valid_factors
    assert all([tf not in ascii_lowercase for tf in valid_factors])


def test_peakpredictor_init(peakpredictor):
    with pytest.raises(ValueError):
        ananse.peakpredictor.PeakPredictor(reference="ananse/db/default_reference")

    p = peakpredictor  # initialized in conftest.py

    # boring stuff
    assert p.data_dir == "ananse/db/default_reference"
    assert p.genome == "tests/data/binding/hg38_testgenome.fa"
    assert p.pfmfile == "tests/data/binding/test.pfm"
    assert p.ncore == 1

    # set in test.motif2factors.txt
    assert p.factors() == ["SOX12"]
    assert p.f2m == {"SOX12": ["GM.5.0.Sox.0001"]}
    assert "GM.5.0.Sox.0001" in p.motifs

    # parsed from input
    assert p.regions == ["chr1:1010-1020"]
    assert p.region_type == "custom"
    assert p.species == "human"  # because hg38 is in the genome name

    assert p._histone_data is None
    assert "chr1:1010-1020" in p._atac_data.index
    assert p._X_columns == ["ATAC", "motif"]
    assert p._model_type == "ATAC_motif"

    # assert len(p.motif_graph) == 707


def test__scan_motifs(peakpredictor):
    region = "chr1:1010-1020"
    tf = "SOX12"
    peakpredictor._scan_motifs([region], zscore=False, gc=False)
    score = peakpredictor._motifs.at[region, tf]
    # rounding to reduce flakiness
    assert round(float(score), 2) == -0.99


def test__load_prescanned_motifs(peakpredictor):
    region = "chr1:1010-1020"
    tf = "SOX12"
    score = peakpredictor._motifs.at[region, tf]
    # rounding to reduce flakiness
    assert round(float(score), 2) == -0.99


@pytest.mark.skip("reference.factor.feather missing")
def test__load_reference_data(peakpredictor):
    p = deepcopy(peakpredictor)
    p._load_reference_data()

    # TODO: nice tests based on the output of these properties
    print(p._motifs)
    print(p._avg)
    print(p._dist)
    print(p.regions)
    exit(1)


def test_factors(peakpredictor):
    p = deepcopy(peakpredictor)
    assert p.factors() == ["SOX12"]

    # dynamically updated
    p.species = None
    p.f2m = {"whatever_I_want": []}
    assert p.factors() == ["whatever_I_want"]

    # filtered
    p.species = "mouse"
    p.f2m = {"mouse_gene": [], "HUMAN_GENE": []}
    assert p.factors() == ["mouse_gene"]


def test__jaccard_motif_graph(peakpredictor):
    factor = "SOX12"
    tfs = nx.single_source_dijkstra_path_length(peakpredictor.motif_graph, factor, 1)
    assert tfs[factor] == 0
    assert "SOX7" in tfs
    assert tfs["SOX7"] == 0.75
