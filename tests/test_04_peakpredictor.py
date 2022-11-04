from types import SimpleNamespace
from copy import deepcopy
import os
import pytest
from string import ascii_lowercase

import pandas as pd
import networkx as nx

import ananse.peakpredictor
from ananse.view import view_h5
from ananse.commands import binding
from tests import write_file


def test_read_factor2motifs():
    f2m = ananse.peakpredictor.read_factor2motifs(
        "tests/data/GRCz11_chr9/GRCz11_chr9.pfm"
    )
    assert len(f2m) == 4
    assert len(f2m["pou3f3a"]) == 18

    with pytest.raises(ValueError):
        ananse.peakpredictor.read_factor2motifs(
            "tests/data/GRCz11_chr9/GRCz11_chr9.pfm", indirect=False
        )

    f2m = ananse.peakpredictor.read_factor2motifs(
        "tests/data/GRCz11_chr9/GRCz11_chr9.pfm", factors=["pou3f3a"]
    )
    assert len(f2m) == 1


def test__prune_f2m():
    f2m = ananse.peakpredictor.read_factor2motifs()

    mouse_myc = ["GM.5.0.bHLH.0008"]
    assert f2m["Myc"] == mouse_myc

    human_myx = [
        "GM.5.0.bZIP.0004",
        "GM.5.0.bHLH.0008",
        "GM.5.0.bHLH.0008",
        "GM.5.0.C2H2_ZF.0024",
        "GM.5.0.bHLH.0026",
        "GM.5.0.bHLH.0030",
        "GM.5.0.bHLH.0037",
        "GM.5.0.bHLH.0052",
        "GM.5.0.bHLH.0058",
        "GM.5.0.bHLH.0079",
        "GM.5.0.bHLH.0101",
        "GM.5.0.bHLH.0114",
        "GM.5.0.bHLH.0117",
        "GM.5.0.bHLH.0122",
        "GM.5.0.bHLH.0131",
        "GM.5.0.bHLH.0138",
    ]
    assert f2m["MYC"] == human_myx

    pruned_f2m = ananse.peakpredictor._prune_f2m(f2m, "mm10")
    assert "MYC" not in pruned_f2m
    assert pruned_f2m["Myc"] == mouse_myc

    pruned_f2m = ananse.peakpredictor._prune_f2m(f2m, "hg38")
    assert "Myc" not in pruned_f2m
    assert set(pruned_f2m["MYC"]).__eq__(set(mouse_myc + human_myx))


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


def test__remove_invalid_regions():
    r_in = ["1", "2", "3"]
    other = ["1", "2", "3"]
    regions = ananse.peakpredictor._remove_invalid_regions(r_in, other)
    assert regions == r_in

    r_in = ["1", "2"]
    other = ["1", "2", "3"]
    regions = ananse.peakpredictor._remove_invalid_regions(r_in, other)
    assert regions == r_in

    r_in = ["1", "2", "3"]
    other = ["1", "2"]
    regions = ananse.peakpredictor._remove_invalid_regions(r_in, other)
    assert regions == other


def test__load_cage_tpm(outdir):
    cage_tpms = os.path.join(outdir, "cage.tsv")
    pd.DataFrame(
        {
            "regions": [
                "chr1:600-930",
                "chr1:610-920",
                "chr1:1010-1020",
                "chr1:1200-1400",
            ],
            "TPM": [10, 20, 0, 0],
        }
    ).to_csv(cage_tpms, sep="\t", index=False)

    df = ananse.peakpredictor._load_cage_tpm(cage_tpms)
    assert df.shape == (3, 1)
    assert df.index.to_list() == ["chr1:1200-1400", "chr1:665-865", "chr1:915-1115"]
    assert df.at["chr1:665-865", "CAGE"] == 15


def test__load_cage_remap_data():
    pass  # TODO


def test__istable():
    assert ananse.peakpredictor._istable("table.tsv")
    assert ananse.peakpredictor._istable(["table.tsv"])

    assert not ananse.peakpredictor._istable(["table.tsv", "and something else"])
    assert not ananse.peakpredictor._istable("not a table")
    assert not ananse.peakpredictor._istable(["not a table"])


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


def test_peakpredictor_init(peakpredictor):
    with pytest.raises(ValueError):
        ananse.peakpredictor.PeakPredictor(reference="ananse/db/default_reference")

    p = peakpredictor  # initialized in conftest.py

    # boring stuff
    assert p.data_dir == os.path.abspath("ananse/db/default_reference")
    assert p.genome == "tests/data/binding/hg38_testgenome.fa"
    assert p.pfmfile == "tests/data/binding/test.pfm"
    assert p.ncore == 1

    # set in test.motif2factors.txt
    assert p.f2m == {"SOX12": ["GM.5.0.Sox.0001"]}

    # parsed from input
    assert p.regions == ["chr1:1010-1020"]
    assert p.factor_model_db == "default"
    assert "chr1:1010-1020" in p.atac_data.index
    assert p.histone_data is None
    assert p.p300_data is None
    assert p.cage_data is None
    assert p._avg is None
    assert p._dist is None
    assert p.all_data is None
    assert p.all_data_columns == ["ATAC", "motif"]

    assert len(p.factor_models) > 0, "factor_models"
    assert len(p.motif_graph) > 0, "motif_graph"


def test__jaccard_motif_graph(peakpredictor):
    factor = "SOX12"
    tfs = nx.single_source_dijkstra_path_length(peakpredictor.motif_graph, factor, 1)
    assert tfs[factor] == 0
    assert "SOX7" in tfs
    assert tfs["SOX7"] == 0.75


def test__load_counts(peakpredictor):
    p = deepcopy(peakpredictor)

    regionsfile = "tests/data/GRCz11_chr9/GRCz11_chr9_regions.bed"
    regions = ananse.peakpredictor.load_regions(regionsfile)
    p.regions = regions
    countsfile = "tests/data/GRCz11_chr9/GRCz11_chr9_raw.tsv"
    df = p._load_counts(countsfile)

    assert len(df) == len(regions)

    regions = ["9:2802-3002"]
    p.regions = regions
    countsfile = "tests/data/GRCz11_chr9/GRCz11_chr9_raw.tsv"
    df = p._load_counts(countsfile)

    assert len(df) == len(regions)


def test__scan_motifs(peakpredictor):
    region = "chr1:1010-1020"
    tf = "SOX12"
    peakpredictor._scan_motifs(zscore=False, gc=False)
    score = peakpredictor.region_factor_df.at[region, tf]
    # rounding to reduce flakiness
    assert round(float(score), 2) == -0.99  # -0.86 new scores with gimme 0.17.0


def test__load_prescanned_motifs(peakpredictor):
    region = "chr1:1010-1020"
    tf = "SOX12"
    score = peakpredictor.region_factor_df.at[region, tf]
    # rounding to reduce flakiness
    assert round(float(score), 2) == -0.99  # -0.86 new scores with gimme 0.17.0


def test__load_cage(outdir, peakpredictor):
    # create some cage data
    cage_tpms = os.path.join(outdir, "cage.tsv")
    pd.DataFrame(
        {
            "regions": [
                "chr1:10-20",  # not matching regions
                "chr1:610-920",  # not matching regions, raw window size
                "chr1:1010-1020",
                "chr1:1200-1400",  # not matching pfmscorefile
            ],
            "TPM": [12.3, 4.56, 7.89, 0],
        }
    ).to_csv(cage_tpms, sep="\t", index=False)

    pfmscorefile = os.path.join(outdir, "pfmscorefile.tsv")
    pd.DataFrame(
        {
            "": ["chr1:10-20", "chr1:760-770", "chr1:1010-1020"],
            "GM.5.0.Sox.0001": [12.3, 4.56, 7.89],
        }
    ).to_csv(pfmscorefile, sep="\t", index=False)

    p = deepcopy(peakpredictor)
    p.genome = "asd"  # we dont want to download 750 MB

    # only pfmfile (with 3 normalized regions)
    p._load_cage(
        cage_tpms=cage_tpms, regions=None, pfmscorefile=pfmscorefile, window=10
    )
    # all regions in the pfmscore file match (after normalizing the tpm regions)
    assert sorted(p.regions) == ["chr1:10-20", "chr1:1010-1020", "chr1:760-770"]
    # CAGE TPM values scaled and ordered correctly
    assert round(p.cage_data.at["chr1:760-770", "CAGE"], 1) == 0.0
    assert round(p.cage_data.at["chr1:1010-1020", "CAGE"], 1) == 0.5
    assert round(p.cage_data.at["chr1:10-20", "CAGE"], 1) == 1.0

    # pfmfile and regions (with 1 matching normalized region)
    p._load_cage(
        cage_tpms=cage_tpms,
        regions=[
            "chr1:1010-1020",
            "chr1:1234-1245",  # invalid
        ],
        pfmscorefile=pfmscorefile,
        window=10,
    )
    # only overlapping regions kept
    assert p.regions == ["chr1:1010-1020"]


def test__load_reference_data(peakpredictor):
    p = deepcopy(peakpredictor)
    with pytest.raises(FileNotFoundError):
        p._load_custom_data()

    # # TODO: nice tests based on the output of these properties
    # print(p.region_factor_df)
    # print(p._avg)
    # print(p._dist)
    # print(p.regions)


def test__model_input(peakpredictor):
    mi = peakpredictor._model_input()
    assert mi == "ATAC_motif"

    p = deepcopy(peakpredictor)
    p.all_data_columns = ["p300", "dist", "CAGE"]
    mi = p._model_input()
    assert mi == "CAGE_H3K27ac_dist_motif"
    assert p.all_data_columns == ["CAGE", "p300", "dist", "motif"]

    p.all_data_columns = ["CAGE", "motif"]
    mi = p._model_input()
    assert mi == "CAGE_average_motif"


def test_predict_binding_probability(peakpredictor):
    pred = peakpredictor.predict_binding_probability("SOX12")
    assert round(pred.at["chr1:1010-1020", 0], 0) == round(1.0, 0)
    assert peakpredictor.all_data.shape == (1, 1)


def test_predict_factor_activity(peakpredictor):
    act = peakpredictor.predict_factor_activity(1)
    assert act.columns.tolist() == ["factor", "activity"]
    assert act.shape == (1, 2)


def test_predict_peaks():
    pass  # TODO


def test_command_binding(outdir, genome):
    out_dir = os.path.join(outdir, "binding")
    bed = os.path.join(outdir, "bed3.bed")
    regions = ["9:2802-3002", "9:3612-3812"]
    tfs = ["pou2f1b", "pou1f1", "pou3f3a"]
    write_file(bed, regions)

    args = SimpleNamespace()  # mimic argparse layer
    args.outdir = out_dir
    args.atac_bams = ["tests/data/GRCz11_chr9/chr9.bam"]
    args.histone_bams = None
    args.cage_tpms = None
    args.columns = None
    args.regions = [bed]  # ["tests/data/GRCz11_chr9/GRCz11_chr9_regions.bed"],
    args.reference = None
    args.tfs = tfs
    args.genome = "tests/data/GRCz11_chr9/GRCz11/GRCz11.fa"
    args.pfmfile = "tests/data/GRCz11_chr9/GRCz11_chr9.pfm"
    args.pfmscorefile = "tests/data/GRCz11_chr9/GRCz11_chr9_scan.bed"
    args.jaccard_cutoff = 0.1
    args.ncore = 1
    binding(args)

    bindingfile = os.path.join(out_dir, "binding.h5")
    assert os.path.exists(bindingfile)
    df = view_h5(bindingfile)
    # regions overlapping between bam files, regionfiles and pfmscorefile
    assert len(df.index) == len(regions)
    # TFs overlapping between factors and pfmfile's motif2factors.txt
    assert len(df.columns) == len(tfs)
