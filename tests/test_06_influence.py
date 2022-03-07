import os

import pandas as pd
import pytest

import ananse.influence


@pytest.fixture
def influence_obj(outdir):
    i = ananse.influence.Influence(
        outfile=os.path.join(outdir, "influence.tsv"),
        degenes="tests/data/influence/degenes.tsv",
        grn_target_file="tests/data/influence/network.tsv",
        edges=10,
        gene_gtf="tests/data/GRCz11_chr9/GRCz11/GRCz11.annotation.gtf",
    )
    return i


def test_read_network():
    grn = ananse.influence.read_network("tests/data/influence/network.tsv", edges=10)
    assert len(grn.nodes) == 10 + 1  # including a self-edge

    top_int = ["FOXK2—AL935186.11", "FOXK2—ABCA7"]
    grn = ananse.influence.read_network(
        "tests/data/influence/network.tsv", interactions=top_int
    )
    assert len(grn.nodes) == 3


def test_read_top_interactions():
    fname = "tests/data/influence/network.tsv"
    top = ananse.influence.read_top_interactions(fname, edges=10)
    assert len(top) == 10

    top = ananse.influence.read_top_interactions(fname, edges=1)
    assert top == {"FOXK2—AL935186.11"}

    top = ananse.influence.read_top_interactions(fname, edges=1, sort_by="activity")
    assert top == {"FOXK2—ABCA7"}


def test_difference():
    pass  # TODO


def test_read_expression(influence_obj):
    scores = influence_obj.expression_change["foxg1b"]
    assert round(scores.score, 2) == 1.83
    assert round(scores.absfc, 2) == 1.83
    assert round(scores.realfc, 2) == 1.83


def test_influence_scores(influence_obj):
    de_genes = set(
        g
        for g in influence_obj.expression_change
        if influence_obj.expression_change[g].score > 0
    )
    de_genes = de_genes & influence_obj.grn.nodes
    line = ananse.influence.influence_scores(
        "FOXK2", influence_obj.grn, influence_obj.expression_change, de_genes
    )
    assert line[0] == "FOXK2"
    assert line[1] == 10  # all test edges
    assert len(line) == 8


def test_filter_tf():
    pass  # TODO


def test_save_reg_network():
    pass  # TODO


def test_run_target_score(outdir):
    i = ananse.influence.Influence(
        outfile=os.path.join(outdir, "influence.tsv"),
        degenes="tests/data/influence/degenes.tsv",
        grn_target_file="tests/data/influence/network.tsv",
        edges=10,
        gene_gtf="tests/data/GRCz11_chr9/GRCz11/GRCz11.annotation.gtf",
        padj_cutoff=0.05,
    )
    i.run_target_score()
    assert os.path.exists(i.outfile)
    df = pd.read_table(i.outfile, index_col=0)
    assert df.loc["FOXK2"]["directTargets"] == 10
    assert round(df.loc["FOXK2"]["pval"], 2) == 0.18

    # same output with multiprocessing
    os.remove(i.outfile)
    i.ncore = 2
    i.run_target_score()
    assert os.path.exists(i.outfile)
    df2 = pd.read_table(i.outfile, index_col=0)
    assert all(df.eq(df2))


def test_run_influence_score():
    pass  # TODO


def test_run_influence():
    pass  # TODO


def test_command_influence():
    pass  # TODO
