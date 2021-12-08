from collections import namedtuple
import math
from copy import deepcopy
import os
import numpy as np
import pandas as pd
import pyranges as pr
import pytest

import ananse.network
from ananse.utils import cleanpath
from tests import write_file

from ananse.commands import network


@pytest.fixture
def network_obj():
    return ananse.network.Network(genome="hg38")


def test_get_bed(outdir):
    cmd = ananse.network.get_bed

    # accepts existing files only
    test_bed = os.path.join(outdir, "test.annotation.bed")
    assert not os.path.exists(test_bed)
    with pytest.raises((ValueError, FileNotFoundError)):
        cmd(test_bed, None)

    test_fa = os.path.join(outdir, "test.fa")
    assert not os.path.exists(test_fa)
    with pytest.raises(FileNotFoundError):
        cmd(None, test_fa)

    # accepts a genomepy genome
    write_file(test_fa, [">chr1\n", "AATTCCGG\n"])
    assert os.path.exists(test_fa)
    assert not os.path.exists(test_bed)
    with pytest.raises(TypeError):
        cmd(None, test_fa)
    write_file(test_bed, "")
    assert os.path.exists(test_bed)
    assert cmd(None, test_fa) == cleanpath(test_bed)

    # accepts a bed file
    assert os.path.exists(test_bed)
    assert cmd(test_bed, None) == test_bed

    # default BED files
    for genome in ["hg38", "hg19"]:
        bed = f"ananse/db/{genome}.genes.bed"
        assert cmd(None, genome) == cleanpath(bed)


def test_get_factors(outdir):
    cmd = ananse.network.get_factors
    with pytest.raises(ValueError):
        cmd()
    with pytest.raises(ValueError):
        cmd(None, None)

    # len(out_tfs) == 0
    with pytest.raises(ValueError):
        cmd(tfs=[])

    # specify tfs
    tfs = ["my_favourite_TF"]
    out_tfs = cmd(tfs=tfs)
    assert out_tfs == tfs

    # binding tfs
    tfs = ["my_favourite_TF", "some_dumb_TF"]
    binding = os.path.join(outdir, "binding.h5")
    with pd.HDFStore(binding, complib="lzo", complevel=9) as hdf:
        for tf in tfs:
            hdf.put(
                key=tf,
                value=pd.Series([0, 0, 0]),
                format="table",
            )
    out_tfs = cmd(binding=binding)
    assert len(set(out_tfs) & set(tfs)) == 2

    # filter binding tfs with specified tfs
    tfs = ["my_favourite_TF"]
    out_tfs = cmd(binding, tfs)
    assert out_tfs == tfs


def test_region_gene_overlap(outdir):
    cmd = ananse.network.region_gene_overlap

    test_bed = os.path.join(outdir, "test.annotation.bed")
    write_file(test_bed, ["chr1\t100\t200\tmy_favourite_TF\t0\t+\t0\t0\t0\t0\t0\t0\n"])
    df = pd.DataFrame(
        {
            "Chromosome": ["chr1"],
            "Start": [300],
            "End": [400],
        }
    )
    region_pr = pr.PyRanges(df)

    # no overlap
    df = cmd(region_pr, test_bed, 0, 0)
    assert len(df) == 0

    # overlap after extension
    df = cmd(region_pr, test_bed, 0, 201)
    assert len(df) == 1


def test_combine_expression_files(outdir):
    cmd = ananse.network.combine_expression_files

    # 1 file
    df1 = pd.DataFrame(
        {
            "gene": ["gene1", "gene2", "gene3"],
            "my_column": [0.0, 100.0, 200.0],  # <- non-default expression column name
        }
    ).set_index("gene")
    f1 = os.path.join(outdir, "expression1.tsv")
    df1.to_csv(f1, sep="\t")
    e1 = cmd(f1, column="my_column")
    assert df1["my_column"].equals(e1["expression"])

    df2 = pd.DataFrame(
        {
            "gene": ["gene1", "gene4"],
            "my_column": [100.0, 100.0],
        }
    ).set_index("gene")
    f2 = os.path.join(outdir, "expression2.tsv")
    df2.to_csv(f2, sep="\t")
    e2 = cmd(f2, column="my_column")
    assert df2["my_column"].equals(e2["expression"])

    # 1 file as list
    e1 = cmd([f1], column="my_column")
    assert df1["my_column"].equals(e1["expression"])

    # 2 files (with mean of duplicate genes)
    exp_files = [f1, f2]
    e3 = cmd(exp_files, column="my_column")
    df3 = pd.DataFrame(
        {
            "gene": ["gene1", "gene2", "gene3", "gene4"],
            "my_column": [50.0, 100.0, 200.0, 100.0],
        }
    ).set_index("gene")
    assert df3["my_column"].equals(e3["expression"])


def test_network_init():
    n = ananse.network.Network(genome="hg38")
    assert cleanpath(n.gene_bed) == cleanpath("ananse/db/hg38.genes.bed")

    n = ananse.network.Network(gene_bed="ananse/db/hg38.genes.bed")
    assert cleanpath(n.gene_bed) == cleanpath("ananse/db/hg38.genes.bed")


# def test_unique_enhancer(network_obj, binding_fname):
#     regions = network_obj.unique_enhancers(binding_fname)
#     regions = regions.as_df()
#     assert regions.shape[0] == 6
#     assert sorted(list(regions["Chromosome"].unique())) == ["chr1", "chr10", "chr17"]
#     assert sorted(list(regions["Start"].unique())) == [7677184, 7687827]


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


def test_enhancer2gene(network_obj):
    # TODO: expand
    df = pd.DataFrame(
        {
            "Chromosome": ["chr1"],
            "Start": [300],
            "End": [400],
        }
    )
    peak_pr = pr.PyRanges(df)
    genes = network_obj.enhancer2gene(peak_pr)
    assert genes.at["chr1:300-400", "gene"] == "OR4F5"


def test_aggregate_binding(network_obj):
    with pytest.raises(ValueError):
        network_obj.aggregate_binding(
            binding="tests/data/network/binding.h5",
            tfs=["GATA3"],
            regions=["chr10:128538766-128538966"],
        )

    df = network_obj.aggregate_binding(
        binding="tests/data/network/binding.h5",
        tfs=["GATA3"],
    ).head()
    assert "weighted_binding" in df.columns
    assert ananse.SEPARATOR in df.index.to_list()[0]
    genes = pd.DataFrame(df.index.str.split(ananse.SEPARATOR).to_list())
    assert set(genes[0]) == {"GATA3"}
    assert len(set(genes[1])) == len(genes[1])


def test_gene_overlap(network_obj):
    # TFs in expression match TFs in gene_bed and TFs in variable tfs
    expression_in = pd.DataFrame({None: ["OR4F5"], "tpm": [10]}).set_index(None)
    tfs = list(expression_in.index)
    expression_out = network_obj.gene_overlap(expression_in, tfs)
    assert expression_out.equals(expression_in)

    # TFs dont overlap, no genomepy annotation
    expression_in = pd.DataFrame({None: ["not_in_there"], "tpm": [10]}).set_index(None)
    tfs = list(expression_in.index)
    with pytest.raises(SystemExit):
        network_obj.gene_overlap(expression_in, tfs)

    # TFs dont overlap, genomepy annotation
    n = deepcopy(network_obj)
    # 2 transcript IDs for 1 gene name (pou3f3a)
    expression_in = pd.DataFrame(
        {None: ["ENSDART00000113914", "ENSDART00000133178"], "tpm": [10, 5]}
    ).set_index(None)
    tfs = ["pou3f3a"]
    n.gene_bed = "tests/data/GRCz11_chr9/GRCz11/GRCz11.annotation.bed"
    expression_out = n.gene_overlap(expression_in, tfs)
    assert expression_in.shape == (2, 1)
    assert "pou3f3a" not in expression_in.index
    assert expression_out.shape == (1, 1)  # duplicate gene names merged
    assert "pou3f3a" in expression_out.index
    assert expression_out.at["pou3f3a", "tpm"] == 15  # duplicate gene names merged


def test_create_expression_network(network_obj):
    column = "my_fav_column"
    tfs = ["tf1", "tf2", "tf3"]
    genes = tfs + ["gene1", "gene2", "gene3"]
    expression = pd.DataFrame(
        {None: genes, column: [math.pow(10, n) for n in range(len(genes))]}
    ).set_index(None)
    df = network_obj.create_expression_network(expression, tfs, column).compute()
    assert len(df) == len(tfs) * len(genes)
    assert len(df.tf_target) == len(set(df.tf_target))
    assert int(df.tf_expression.min()) == 0
    assert int(df.tf_expression.max()) == 1
    assert int(df.target_expression.min()) == 0
    assert int(df.target_expression.max()) == 1


def test_run_network(network_obj, outdir):
    n = deepcopy(network_obj)

    # expression network (minimal)
    n.full_output = False
    df = n.run_network(
        binding=None,
        fin_expression="tests/data/network/heart_expression.tsv",
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
    )
    assert df.index.name is None
    assert sorted(df.columns) == sorted(["tf_target", "prob"])

    # expression network (full)
    n.full_output = True
    df = n.run_network(
        binding=None,
        fin_expression="tests/data/network/heart_expression.tsv",
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
    )
    assert df.index.name is None
    assert sorted(df.columns) == sorted(
        ["tf_target", "prob", "tf_expression", "target_expression"]
    )

    # binding network (minimal)
    n.full_output = False
    df = n.run_network(
        binding="tests/data/network/binding.h5",
        fin_expression=None,
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
    )
    assert df.index.name is None
    assert sorted(df.columns) == sorted(["tf_target", "prob"])

    # binding network (full)
    n.full_output = True
    df = n.run_network(
        binding="tests/data/network/binding.h5",
        fin_expression=None,
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
    )
    assert df.index.name is None
    assert sorted(df.columns) == sorted(["tf_target", "prob", "weighted_binding"])

    # expression-binding network (minimal)
    n.full_output = False
    df = n.run_network(
        binding="tests/data/network/binding.h5",
        fin_expression="tests/data/network/heart_expression.tsv",
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
    )
    assert df.index.name is None
    assert sorted(df.columns) == sorted(["tf_target", "prob"])

    # expression-binding network (full)
    n.full_output = True
    df = n.run_network(
        binding="tests/data/network/binding.h5",
        fin_expression="tests/data/network/heart_expression.tsv",
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
    )
    assert df.index.name is None
    assert sorted(df.columns) == sorted(
        [
            "tf_target",
            "prob",
            "tf_expression",
            "target_expression",
            "weighted_binding",
            "activity",
        ]
    )

    # expression-binding network without enhancers/promoters + outfile
    n.full_output = True
    n.include_promoter = False
    n.include_enhancer = False
    outfile = os.path.join(outdir, "network.tsv")
    n.run_network(
        binding="tests/data/network/binding.h5",
        fin_expression="tests/data/network/heart_expression.tsv",
        tfs=["GATA3"],
        regions=["chr10:11716757-11716957"],
        outfile=outfile,
    )
    df = pd.read_csv(outfile, sep="\t")
    assert df.index.name is None
    assert sorted(df.columns) == sorted(
        ["tf_target", "prob", "tf_expression", "target_expression", "activity"]
    )


def test_command_network(outdir):
    outfile = os.path.join(outdir, "network.tsv")
    expression = "tests/data/network/heart_expression.tsv"
    tfs = ["GATA3"]
    regionsfile = os.path.join(outdir, "regions.txt")
    write_file(regionsfile, ["chr10:11716757-11716957\n", "chr10:10-100\n"])
    Args = namedtuple(
        "args",
        "genome annotation include_promoter include_enhancer binding "
        "fin_expression tfs regions column outfile full_output ncore",
    )
    args = Args(
        genome="hg38",
        annotation=None,
        include_promoter=False,
        include_enhancer=False,
        binding="tests/data/network/binding.h5",
        fin_expression=expression,
        tfs=tfs,
        regions=regionsfile,
        column="tpm",
        full_output=False,
        outfile=outfile,
        ncore=1,
    )
    network(args)

    df = pd.read_table(outfile, sep="\t")
    assert sorted(df.columns) == sorted(["tf_target", "prob"])
    tf_target = df["tf_target"].str.split(ananse.SEPARATOR, expand=True)
    genes = pd.read_csv(expression, sep="\t", usecols=[0])["target_id"]
    assert set(tf_target[0]) == set(tfs)
    assert len(tf_target[1]) == len(set(tf_target[1]))
    assert len(tf_target[1]) <= len(set(genes))
    assert all(tf_target[1].isin(genes))
