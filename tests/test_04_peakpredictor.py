from collections import namedtuple
from copy import deepcopy
import os
import pytest
from string import ascii_lowercase

import networkx as nx

import ananse.peakpredictor
from ananse.view import view_h5
from ananse.commands import binding
from tests import write_file


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


def test__load_reference_data(peakpredictor):
    p = deepcopy(peakpredictor)
    with pytest.raises(FileNotFoundError):
        p._load_reference_data()

    # # TODO: nice tests based on the output of these properties
    # print(p._motifs)
    # print(p._avg)
    # print(p._dist)
    # print(p.regions)
    # exit(1)


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


def test_command_binding(outdir, genome):
    Args = namedtuple(
        "args",
        "outdir atac_bams histone_bams regionfiles reference factors genome pfmfile pfmscorefile jaccard_cutoff ncore",
    )
    out_dir = os.path.join(outdir, "binding")
    bed = os.path.join(outdir, "bed3.bed")

    # # hg38
    # write_file(
    #     bed,
    #     [
    #         "chr1\t50\t250\n",  # regions encompasses the bam reads
    #         "chr1\t70\t300\n",
    #         "chr1\t450\t500\n",
    #     ],
    # )
    # scorefile = os.path.join(outdir, "scorefile.tsv")
    # write_file(
    #     scorefile,
    #     [
    #         "\tGM.5.0.Sox.0001\n",
    #         "chr1:50-250\t0.2\n",
    #         "chr1:70-300\t0.0\n",
    #         "chr1:450-500\t-0.2\n",
    #     ],
    # )
    # args = Args(
    #     outdir=out_dir,
    #     atac_bams=["tests/data/binding/bam1.bam"],  # chr 1, 3 reads
    #     histone_bams=None,
    #     regionfiles=[bed],
    #     reference=None,
    #     factors=None,
    #     genome="tests/data/binding/hg38_testgenome.fa",  # chr 1 (fake)
    #     pfmfile="tests/data/binding/test.pfm",  # 1 motif (GM.5.0.Sox.0001), 1 factor (SOX12)
    #     pfmscorefile=scorefile,
    #     jaccard_cutoff=0.1,
    #     ncore=1,
    # )
    # binding(args)
    #
    # bindingfile = os.path.join(out_dir, "binding.h5")
    # assert os.path.exists(bindingfile)
    # df = view_h5(bindingfile)
    # # regions overlapping between bam files, regionfiles and pfmscorefile
    # assert set(df.index) == {"chr1:50-250", "chr1:70-300", "chr1:450-500"}
    # # TFs overlapping between factors and pfmfile's motif2factors.txt
    # assert set(df.columns) == {"SOX12"}

    # GRCz11
    write_file(
        bed,
        [
            "9:2802-3002\n",
            "9:3612-3812\n",
        ],
    )
    args = Args(
        outdir=out_dir,
        atac_bams=["tests/data/GRCz11_chr9/chr9.bam"],
        histone_bams=None,
        regionfiles=[bed],  # ["tests/data/GRCz11_chr9/GRCz11_chr9_regions.bed"],
        reference=None,
        factors=["pou2f1b", "pou1f1", "pou3f3a"],
        genome="tests/data/GRCz11_chr9/GRCz11/GRCz11.fa",
        pfmfile="tests/data/GRCz11_chr9/GRCz11_chr9.pfm",
        pfmscorefile="tests/data/GRCz11_chr9/GRCz11_chr9_scan.bed",
        jaccard_cutoff=0.1,
        ncore=1,
    )
    binding(args)

    bindingfile = os.path.join(out_dir, "binding.h5")
    assert os.path.exists(bindingfile)
    df = view_h5(bindingfile)
    # regions overlapping between bam files, regionfiles and pfmscorefile
    assert len(df.index) == 2
    # TFs overlapping between factors and pfmfile's motif2factors.txt
    assert len(df.columns) == 3
