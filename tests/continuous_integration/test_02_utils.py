import os
import time

import genomepy.utils
import pysam

import ananse.utils

# prep

test_dir = os.path.dirname(os.path.dirname(__file__))
outdir = os.path.join(test_dir, "output")
genomepy.utils.mkdir_p(outdir)


def write_file(filename, lines):
    with open(filename, "w") as f:
        for line in lines:
            if not line.endswith("\n"):
                line = line + "\n"
            f.write(line)


def compare_contents(file1, file2, ftype="bed"):
    if ftype == "bed":
        with open(file1) as f:
            contents1 = f.readlines()
        with open(file2) as f:
            contents2 = f.readlines()
    else:
        contents1 = pysam.view(file1)
        contents2 = pysam.view(file2)
    return contents1 == contents2


# test BED functions


unsorted_bed = os.path.join(outdir, "unsorted.bed")
write_file(unsorted_bed, ["chr1\t817046\t817246\n", "chr1\t778558\t778758\n"])

sorted_bed = os.path.join(outdir, "sorted.bed")
write_file(sorted_bed, ["chr1\t778558\t778758\n", "chr1\t817046\t817246\n"])

second_bed = os.path.join(outdir, "second.bed")
write_file(second_bed, ["chr1\t827457\t827657\n"])


def test_bed_sort():
    assert not compare_contents(unsorted_bed, sorted_bed, ftype="bed")
    ananse.utils.bed_sort(unsorted_bed)
    assert compare_contents(unsorted_bed, sorted_bed, ftype="bed")


def test_bed_merge():
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
    ananse.utils.bed_merge([unsorted_bed, second_bed], merged_bed)
    with open(merged_bed) as mb:
        assert len(mb.readlines()) == 3


# test BAM functions


h0 = "@HD	VN:1.6	SO:coordinate"
h1 = "@SQ	SN:chr1	LN:248956422"
line1 = (
    "GADGET:53:HF5JFBBXX:2:1228:15189:39594	147	chr1	10003	40	11S90M	=	10048	-46	"
    + "CCCTACCCTCTCCCTATCCCTAACCCTAACCCCAACCCTAACCCTATCCCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA	"
    + "A77--7-7---7-7---A77---AA7----<7-AAJAA-7JJFF<--F-A-AFFFF<FJJJJF-AFJF7F-JJJFJFFFJFF<FJJJJFJJFJJFFFFFAA	"
    + "MC:Z:58M1I42M	MD:Z:5A15T13A3T50	PG:Z:MarkDuplicates	NM:i:4	AS:i:70	XS:i:70"
)
line2 = (
    "GADGET:53:HF5JFBBXX:2:2224:31010:36622	83	chr1	10004	30	2S45M1D54M	=	10074	-30	"
    + "ATCCCTAACCCTAACCCTAACCCTAACCCTACCCCTACCCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT	"
    + "--JAA7F-FAFA-7JJFA--F<7-FF<<FAF7<7F7A-FFAF7-FJJJFJJ----J<JFA-JAF7JFJFJF<<JFJF<JJJFFJJJAAAA-JFFFA-FAA-	"
    + "XA:Z:chr3,-10471,31M1I69M,4;chr3,-10535,101M,5;chr21,+46699903,82M19S,2;	MC:Z:73M1D27M	"
    + "MD:Z:29A5A9^T54	PG:Z:MarkDuplicates	NM:i:3	AS:i:82	XS:i:81"
)
line3 = (
    "GADGET:53:HF5JFBBXX:3:1110:5051:22168	163	chr1	10027	40	100M	=	10032	105	"
    + "ACCCGAACCCTAACCCTAACCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCGAACCCA	"
    + "AAFFFJJJJJJJJJJJFJJJFJJJFJFJJFJJJJ<-FJJFJAFFJA7AFAJJJJFJFJ-<F-AAJJ<FF7-J-AAJ--<JJJ--AAJ-77-AA-7A<-A-	"
    + "MC:Z:8M1D91M	MD:Z:4T29T58T5T0	PG:Z:MarkDuplicates	NM:i:4	AS:i:84	XS:i:82"
)

unsorted_bam = os.path.join(outdir, "unsorted.bam")
tmp_sam = os.path.join(outdir, "tmp.sam")
write_file(tmp_sam, [h0, h1, line2, line1])
pysam.view(tmp_sam, "-b", "-o", unsorted_bam, catch_stdout=False)

sorted_bam = os.path.join(outdir, "sorted.bam")
write_file(tmp_sam, [h0, h1, line1, line2])
pysam.view(tmp_sam, "-b", "-o", sorted_bam, catch_stdout=False)

second_bam = os.path.join(outdir, "second.bam")
write_file(tmp_sam, [h0, h1, line3])
pysam.view(tmp_sam, "-b", "-o", second_bam, catch_stdout=False)

genomepy.utils.rm_rf(tmp_sam)


def test_bam_index():
    ncores = os.cpu_count()  # test max cores

    genomepy.utils.rm_rf(f"{sorted_bam}.bai")
    assert not os.path.exists(f"{sorted_bam}.bai")

    ananse.utils.bam_index(sorted_bam, ncore=ncores)
    assert os.path.exists(f"{sorted_bam}.bai")

    # test force
    t0 = os.path.getmtime(f"{sorted_bam}.bai")
    time.sleep(1)

    ananse.utils.bam_index(sorted_bam, force=False, ncore=ncores)
    t1 = os.path.getmtime(f"{sorted_bam}.bai")
    assert t0 == t1

    ananse.utils.bam_index(sorted_bam, force=True, ncore=ncores)
    t1 = os.path.getmtime(f"{sorted_bam}.bai")
    assert t0 != t1


def test_bam_sort():
    ncores = -999  # test min cores

    assert not compare_contents(sorted_bam, unsorted_bam, ftype="bam")
    ananse.utils.bam_sort(unsorted_bam, ncore=ncores)
    assert compare_contents(sorted_bam, unsorted_bam, ftype="bam")
    assert os.path.exists(f"{unsorted_bam}.bai")  # bam is indexed

    # bam is identical to the already sorted bam
    ananse.utils.bam_index(sorted_bam, force=False)
    assert os.path.getsize(f"{unsorted_bam}.bai") == os.path.getsize(
        f"{sorted_bam}.bai"
    )


def test_bam_merge():
    ncores = min(2, os.cpu_count())  # test average cores
    merged_bam = os.path.join(outdir, "merged.bam")

    # 1 bam: copy
    ananse.utils.bam_merge([sorted_bam], merged_bam, ncore=ncores)
    assert compare_contents(sorted_bam, merged_bam, ftype="bam")
    assert os.path.getsize(f"{sorted_bam}.bai") == os.path.getsize(f"{merged_bam}.bai")

    # >1 bam: merge
    ananse.utils.bam_merge([sorted_bam, second_bam], merged_bam, ncore=ncores)
    l1 = pysam.view(sorted_bam).strip().split("\n")
    l2 = pysam.view(second_bam).strip().split("\n")
    l3 = pysam.view(merged_bam).strip().split("\n")
    assert len(l1) + len(l2) == len(l3) == 3


def test_cleanpath():
    path = "./tests/continuous_integration/test_02_utils.py"
    expected = __file__
    res = ananse.utils.cleanpath(path)
    assert res == expected

    path = "~/../.."
    expected = "/"
    res = ananse.utils.cleanpath(path)
    assert res == expected


# test other functions


def test_cleanup():
    genomepy.utils.rm_rf(outdir)
