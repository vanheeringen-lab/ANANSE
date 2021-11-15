import os
import tempfile

import pytest
from gimmemotifs.motif import read_motifs

import ananse.utils
from . import compare_contents, write_file


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


def test_bed_sort(unsorted_bed, sorted_bed):
    assert not compare_contents(unsorted_bed, sorted_bed, ftype="bed")
    ananse.utils.bed_sort(unsorted_bed)
    assert compare_contents(unsorted_bed, sorted_bed, ftype="bed")


def test_bed_merge(outdir, sorted_bed, unsorted_bed):
    merged_bed = os.path.join(outdir, "merged.bed")

    # 1 bed = nothing changes
    ananse.utils.bed_merge([sorted_bed], merged_bed)
    assert compare_contents(sorted_bed, merged_bed, ftype="bed")

    # >1 bed, same content
    ananse.utils.bed_merge([unsorted_bed, sorted_bed], merged_bed)
    assert compare_contents(sorted_bed, merged_bed, ftype="bed")
    with open(merged_bed) as mb:
        assert len(mb.readlines()) == 2

    # >1 beds, different content
    second_bed = os.path.join(outdir, "second.bed")
    write_file(second_bed, ["chr1\t827457\t827657\n"])

    ananse.utils.bed_merge([unsorted_bed, second_bed], merged_bed)
    with open(merged_bed) as mb:
        assert len(mb.readlines()) == 3


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
    path = "./tests/continuous_integration/test_02_utils.py"
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


def test_get_binding_tfs():
    cmd = ananse.utils.get_binding_tfs
    binding = "tests/data/network/binding.h5"
    ret = cmd(binding)
    assert len(ret) == 30
    assert "NFKB2" in ret


def test_view_h5():
    cmd = ananse.utils.view_h5
    binding = "tests/data/network/binding.h5"

    # regions and TFs
    regions = cmd(binding, list_regions=True)
    assert regions.shape == (61786, 1)
    tfs = cmd(binding, list_tfs=True)
    assert tfs.shape == (30, 1)

    # all = regions x TFs
    complete_view = cmd(binding)
    assert complete_view.shape == (regions.shape[0], tfs.shape[0])

    # various subsets in either format
    head_wide = cmd(binding, n=10)
    assert head_wide.shape == (10, 10)
    head_long = cmd(binding, n=10, fmt="long")
    assert head_long.shape == (10 * 10, 3)

    subset_wide = cmd(
        binding,
        tfs=["NFKB2", "KLF6"],
        regions=[
            "chr10:100000171-100000371",
            "chr10:100001843-100002043",
        ],
    )
    assert subset_wide.shape == (2, 2)
    subset_long = cmd(
        binding,
        tfs=["NFKB2", "KLF6"],
        regions=[
            "chr10:100000171-100000371",
            "chr10:100001843-100002043",
        ],
        fmt="long",
    )
    assert subset_long.shape == (2 * 2, 3)
