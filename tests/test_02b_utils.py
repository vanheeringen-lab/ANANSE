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


# # test BAM functions
#
#
# h0 = "@HD	VN:1.6	SO:coordinate"
# h1 = "@SQ	SN:chr1	LN:50000"
# line1 = (
#     "read1	147	chr1	10003	40	11S90M	=	10048	-46	"
#     + "CCCTACCCTCTCCCTATCCCTAACCCTAACCCCAACCCTAACCCTATCCCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA	"
#     + "A77--7-7---7-7---A77---AA7----<7-AAJAA-7JJFF<--F-A-AFFFF<FJJJJF-AFJF7F-JJJFJFFFJFF<FJJJJFJJFJJFFFFFAA	"
# )
# line2 = (
#     "read2	83	chr1	10004	30	2S45M1D54M	=	10074	-30	"
#     + "ATCCCTAACCCTAACCCTAACCCTAACCCTACCCCTACCCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT	"
#     + "--JAA7F-FAFA-7JJFA--F<7-FF<<FAF7<7F7A-FFAF7-FJJJFJJ----J<JFA-JAF7JFJFJF<<JFJF<JJJFFJJJAAAA-JFFFA-FAA-	"
# )
# line3 = (
#     "read3	163	chr1	10027	40	100M	=	10032	105	"
#     + "ACCCGAACCCTAACCCTAACCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCGAACCCA	"
#     + "AAFFFJJJJJJJJJJJFJJJFJJJFJFJJFJJJJ<-FJJFJAFFJA7AFAJJJJFJFJ-<F-AAJJ<FF7-J-AAJ--<JJJ--AAJ-77-AA-7A<-A-	"
# )
#
# unsorted_bam = os.path.join(outdir, "unsorted.bam")
# write_bam(unsorted_bam, [h0, h1, line2, line1])
#
# sorted_bam = os.path.join(outdir, "sorted.bam")
# write_bam(sorted_bam, [h0, h1, line1, line2])
#
# second_bam = os.path.join(outdir, "second.bam")
# write_bam(second_bam, [h0, h1, line3])
#
#
# def test_bam_index():
#     ncores = os.cpu_count()  # test max cores
#
#     genomepy.utils.rm_rf(f"{sorted_bam}.bai")
#     assert not os.path.exists(f"{sorted_bam}.bai")
#
#     ananse.utils.bam_index(sorted_bam, ncore=ncores)
#     assert os.path.exists(f"{sorted_bam}.bai")
#
#     # test force
#     t0 = os.path.getmtime(f"{sorted_bam}.bai")
#     time.sleep(1)
#
#     ananse.utils.bam_index(sorted_bam, force=False, ncore=ncores)
#     t1 = os.path.getmtime(f"{sorted_bam}.bai")
#     assert t0 == t1
#
#     ananse.utils.bam_index(sorted_bam, force=True, ncore=ncores)
#     t1 = os.path.getmtime(f"{sorted_bam}.bai")
#     assert t0 != t1
#
#
# def test_bam_sort():
#     ncores = -999  # test min cores
#
#     assert not compare_contents(sorted_bam, unsorted_bam, ftype="bam")
#     ananse.utils.bam_sort(unsorted_bam, ncore=ncores)
#     assert compare_contents(sorted_bam, unsorted_bam, ftype="bam")
#     assert os.path.exists(f"{unsorted_bam}.bai")  # bam is indexed
#
#     # bam is identical to the already sorted bam
#     ananse.utils.bam_index(sorted_bam, force=False)
#     assert os.path.getsize(f"{unsorted_bam}.bai") == os.path.getsize(
#         f"{sorted_bam}.bai"
#     )
#
#
# def test_bam_merge():
#     ncores = min(2, os.cpu_count())  # test average cores
#     merged_bam = os.path.join(outdir, "merged.bam")
#
#     # 1 bam: copy
#     ananse.utils.bam_merge([sorted_bam], merged_bam, ncore=ncores)
#     assert compare_contents(sorted_bam, merged_bam, ftype="bam")
#     assert os.path.getsize(f"{sorted_bam}.bai") == os.path.getsize(f"{merged_bam}.bai")
#
#     # >1 bam: merge
#     ananse.utils.bam_merge([sorted_bam, second_bam], merged_bam, ncore=ncores)
#     l1 = pysam.view(sorted_bam).strip().split("\n")
#     l2 = pysam.view(second_bam).strip().split("\n")
#     l3 = pysam.view(merged_bam).strip().split("\n")
#     assert len(l1) + len(l2) == len(l3) == 3
#
#
# def test_mosdepth():
#     bed_input = os.path.join(outdir, "mosdepth_input.bed")
#     write_file(bed_input, ["chr1\t10003\t10203\n", "chr1\t10203\t10403\n"])
#
#     # bam = sorted & indexed (required)
#     bam_input = os.path.join(outdir, "mosdepth_input.bam")
#     write_bam(bam_input, [h0, h1, line1, line2, line3])
#     ananse.utils.bam_index(bam_input, ncore=os.cpu_count())
#
#     bed_output = os.path.join(outdir, "mosdepth_output.bed")
#     ananse.utils.mosdepth(bed_input, bam_input, bed_output, ncore=1)
#
#     with open(bed_output) as f:
#         score = f.readlines()[0].strip().split("\t")[3]
#     assert score == "1.00"


# test other functions


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


def test_check_input_factors(outdir):
    cmd = ananse.utils.check_input_factors
    assert cmd(None) is None

    # accepts a list of TFs (len>1)
    in_tfs = ["TF1", "TF2"]
    output = cmd(in_tfs)
    assert in_tfs == output

    # accepts a file
    in_file = os.path.join(outdir, "tfs.txt")
    write_file(in_file, ["TF1\n", "TF2\n"])
    output = cmd(in_file)
    assert in_tfs == output

    # accepts a list of files len==1
    output = cmd([in_file])
    assert in_tfs == output


def test__check_input_regions(outdir, genome, bed1, bed2):
    regions = ananse.utils.check_input_regions(None, None, verbose=True, force=True)
    assert regions is None

    # >1 file (all BED)
    regionfiles = [bed1, bed2]
    regions = ananse.utils.check_input_regions(
        regionfiles, genome, outdir, verbose=True, force=True
    )
    assert sorted(regions) == sorted(
        ["chr1:400-600", "chr1:2400-2600", "chr1:4400-4600"]
    )

    # 1 file (BED)
    bed3 = os.path.join(outdir, "bed3.bed")
    write_file(bed3, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n", "chr1\t4000\t5000\n"])

    regions = ananse.utils.check_input_regions([bed3], None, verbose=True, force=True)
    assert sorted(regions) == sorted(
        ["chr1:0-1000", "chr1:2000-3000", "chr1:4000-5000"]
    )

    # 1 file (regions /w header)
    bed4 = os.path.join(outdir, "bed4.bed")
    write_file(
        bed4, ["header", "chr1:0-1000\n", "chr1:2000-3000\n", "chr1:4000-5000\n"]
    )

    regions = ananse.utils.check_input_regions([bed4], None, verbose=True, force=True)
    assert sorted(regions) == sorted(
        ["chr1:0-1000", "chr1:2000-3000", "chr1:4000-5000"]
    )
