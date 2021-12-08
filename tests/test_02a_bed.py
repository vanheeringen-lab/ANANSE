import os

import ananse.bed
from tests import compare_contents, write_file


def test_bed_sort(unsorted_bed, sorted_bed):
    assert not compare_contents(unsorted_bed, sorted_bed, ftype="bed")
    ananse.bed.bed_sort(unsorted_bed)
    assert compare_contents(unsorted_bed, sorted_bed, ftype="bed")


def test_bed_merge(outdir, sorted_bed, unsorted_bed):
    merged_bed = os.path.join(outdir, "merged.bed")

    # 1 bed = nothing changes
    ananse.bed.bed_merge([sorted_bed], merged_bed)
    assert compare_contents(sorted_bed, merged_bed, ftype="bed")

    # >1 bed, same content
    ananse.bed.bed_merge([unsorted_bed, sorted_bed], merged_bed)
    assert compare_contents(sorted_bed, merged_bed, ftype="bed")
    with open(merged_bed) as mb:
        assert len(mb.readlines()) == 2

    # >1 beds, different content
    second_bed = os.path.join(outdir, "second.bed")
    write_file(second_bed, ["chr1\t827457\t827657\n"])

    ananse.bed.bed_merge([unsorted_bed, second_bed], merged_bed)
    with open(merged_bed) as mb:
        assert len(mb.readlines()) == 3


def test_is_narrowpeak(outdir, genome):
    np = os.path.join(outdir, "f.narrowPeak")
    write_file(np, ["chr1\t629812\t630105\tnarrowPeak1\t6047\t.\t0\t0\t0\t122"])
    bp = os.path.join(outdir, "f.broadPeak")
    write_file(
        bp,
        [
            "chr1\t778061\t779255\tbroadRegion1\t660\t.\t778061\t"
            + "779255\t0\t3\t1,638,1\t0,17,1193\t0\t0\t0"
        ],
    )

    cbed = ananse.bed.CombineBedFiles(genome=genome, peakfiles=[])
    assert cbed.is_narrowpeak(np) is True
    assert cbed.is_narrowpeak(bp) is False


def test_bed_resize(outdir, genome, bed1):
    cbed = ananse.bed.CombineBedFiles(genome=genome, peakfiles=[])
    bed_out = os.path.join(outdir, "bed_out.bed")

    # default width, extended width with outlier, extended width with fixed outlier
    for n in range(3):
        width = [200, 2000, 2000][n]
        fix_outliers = [False, False, True][n]
        nlines = [2, 1, 2][n]
        estart = [400, 1500, 0][n]
        estop = [600, 3500, 2000][n]

        cbed.bed_resize(genome, bed1, bed_out, width=width, fix_outliers=fix_outliers)
        with open(bed_out) as f:
            lines = f.readlines()
        assert len(lines) == nlines
        chrom, start, stop = lines[0].split()[0:3]
        assert int(start) == estart
        assert int(stop) == estop
        assert int(stop) - int(start) == width


def test_cbedf(outdir, genome, bed1, bed2):
    cbed = ananse.bed.CombineBedFiles(genome=genome, peakfiles=[bed1, bed2])
    combined_bed = os.path.join(outdir, "combined.bed")
    width = 200
    cbed.run(outfile=combined_bed, width=width, force=True)

    with open(combined_bed) as f:
        lines = f.readlines()

    # 3 unique regions over the 2 bed files
    assert len(lines) == 3

    # width is set correctly
    for line in lines:
        chrom, start, stop = line.split()[0:3]
        assert int(stop) - int(start) == width
