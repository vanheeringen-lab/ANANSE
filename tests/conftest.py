"""
Global fixtures and functions for pytest
pytest can only share fixtures between modules if they are declared here.
"""
import os
import pytest

import ananse.peakpredictor
from tests import write_file


@pytest.fixture(scope="session")
def outdir(tmp_path_factory):
    """tmp directory"""
    outdir = tmp_path_factory.mktemp("pytest_output")
    return outdir


# for utils


@pytest.fixture(scope="session")
def unsorted_bed(outdir):
    unsorted_bed = os.path.join(outdir, "unsorted.bed")
    write_file(unsorted_bed, ["chr1\t817046\t817246\n", "chr1\t778558\t778758\n"])
    return unsorted_bed


@pytest.fixture(scope="session")
def sorted_bed(outdir):
    sorted_bed = os.path.join(outdir, "sorted.bed")
    write_file(sorted_bed, ["chr1\t778558\t778758\n", "chr1\t817046\t817246\n"])
    return sorted_bed


# for enhancer binding & peakpredictor


@pytest.fixture(scope="session")
def genome(outdir):
    genome = os.path.join(outdir, "genome.fa")
    write_file(genome, [">chr1", "N" * 50000])
    return genome


@pytest.fixture(scope="session")
def bed1(outdir):
    bed1 = os.path.join(outdir, "bed1.bed")
    write_file(bed1, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n"])
    return bed1


@pytest.fixture(scope="session")
def bed2(outdir):
    bed2 = os.path.join(outdir, "bed2.bed")
    write_file(bed2, ["chr1\t4000\t5000\n", "chr1\t2000\t3000\n"])
    return bed2


@pytest.fixture(scope="package")
def peakpredictor():
    pp = ananse.peakpredictor.PeakPredictor(
        reference="ananse/db/default_reference",
        atac_bams=["tests/data/binding/bam1.bam"],
        regions=["chr1:1010-1020"],  # 1 (tiny) region
        genome="tests/data/binding/hg38_testgenome.fa",  # 1 chr
        pfmfile="tests/data/binding/test.pfm",  # 1 motif, 1 factor
        pfmscorefile="tests/data/binding/scan.tsv",  # precomputed
        ncore=1,
    )
    return pp
