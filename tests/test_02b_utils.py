import os
import tempfile

import pytest
from gimmemotifs.motif import read_motifs

import ananse.utils
from tests import write_file


def test_check_path():
    cmd = ananse.utils.check_path
    assert cmd(None) is None

    # accepts a string, returns a string
    abspath = cmd("tests", error_missing=False)
    assert os.path.exists(abspath) and os.path.isabs(abspath)

    # accepts a list, returns a list
    relpaths = ["tests/data", "tests"]
    abspaths = cmd(relpaths, error_missing=False)
    assert all(os.path.exists(p) for p in abspaths)
    assert all(os.path.isabs(p) for p in abspaths)

    # error on missing file/dir
    with pytest.raises(FileNotFoundError):
        cmd("missing_dir")


def test_cleanpath():
    path = "./tests/test_02b_utils.py"
    expected = __file__
    res = ananse.utils.cleanpath(path)
    assert res == expected

    path = "~/../.."
    expected = "/"
    res = ananse.utils.cleanpath(path)
    assert res == expected


def test_mytmpdir():
    tmpdir = ananse.utils.mytmpdir()
    assert os.path.exists(tmpdir)
    assert tempfile.gettempdir() in tmpdir


def test_clean_tmp():
    tmpdir = ananse.utils.mytmpdir()
    assert os.path.exists(tmpdir)
    ananse.utils.clean_tmp()
    assert not os.path.exists(tmpdir)


def test_get_motif_factors():
    pfmfile = "tests/data/debug.pfm"
    motifs = read_motifs(pfmfile, as_dict=True)
    motif = motifs["GM.5.0.Homeodomain.0001"]
    cmd = ananse.utils.get_motif_factors

    indirect_tfs = cmd(motif, indirect=True)
    assert indirect_tfs == ["TGIF1"]

    direct_tfs = cmd(motif, indirect=False)
    assert direct_tfs == []


def test_load_tfs(outdir):
    cmd = ananse.utils.load_tfs
    assert cmd(None) is None

    # accepts a TF
    in_tfs = "TF1"
    output = cmd(in_tfs)
    assert [in_tfs] == output

    # accepts a list of TFs
    in_tfs = ["TF1", "TF2"]
    output = cmd(in_tfs)
    assert in_tfs == output

    # accepts a file
    in_file = os.path.join(outdir, "tfs.txt")
    write_file(in_file, ["TF1\n", "TF2\n"])
    output = cmd(in_file)
    assert in_tfs == output

    # accepts a list of files
    output = cmd([in_file])
    assert in_tfs == output
    output = cmd([in_file, in_file])
    assert in_tfs == output

    # error on bad file
    in_tfs = "/file"
    with pytest.raises(ValueError):
        cmd(in_tfs)


def test_load_regions(outdir, genome, bed1, bed2):
    cmd = ananse.utils.load_regions

    regions = cmd(None, None, None)
    assert regions is None

    # >1 file (all BED)
    regionfiles = [bed1, bed2]
    regions = cmd(regionfiles, genome, outdir)
    assert regions == sorted(["chr1:400-600", "chr1:2400-2600", "chr1:4400-4600"])

    # 1 file (BED)
    bed3 = os.path.join(outdir, "bed3.bed")
    write_file(bed3, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n", "chr1\t4000\t5000\n"])

    regions = cmd(bed3, None, None)
    assert regions == sorted(["chr1:0-1000", "chr1:2000-3000", "chr1:4000-5000"])

    # 1 file (regions /w header)
    bed4 = os.path.join(outdir, "bed4.bed")
    write_file(
        bed4, ["header", "chr1:0-1000\n", "chr1:2000-3000\n", "chr1:4000-5000\n"]
    )

    regions = cmd(bed4, None, None)
    assert regions == sorted(["chr1:0-1000", "chr1:2000-3000", "chr1:4000-5000"])

    # error on bad regions
    regions = "not:a-region"
    with pytest.raises(ValueError):
        cmd(regions, None, None)
