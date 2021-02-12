import os
import filecmp

import genomepy.utils

import ananse.utils

# prep
test_dir = os.path.dirname(__file__)
data_dir = os.path.join(test_dir, "data")
genomepy.utils.mkdir_p(data_dir)
outdir = os.path.join(test_dir, "output")

unsorted_bed = os.path.join(outdir, "unsorted.bed")
with open(unsorted_bed, "w") as bed:
    bed.write("chr1\t817046\t817246\n")
    bed.write("chr1\t778558\t778758\n")

sorted_bed = os.path.join(outdir, "sorted.bed")
with open(unsorted_bed, "w") as bed:
    bed.write("chr1\t778558\t778758\n")
    bed.write("chr1\t817046\t817246\n")

second_bed = os.path.join(outdir, "second.bed")
with open(unsorted_bed, "w") as bed:
    bed.write("chr1\t827457\t827657\n")


def test_bed_sort():
    assert not filecmp.cmp(unsorted_bed, sorted_bed)
    ananse.utils.bed_sort(unsorted_bed)

    assert filecmp.cmp(unsorted_bed, sorted_bed)


def test_bed_merge():
    merged_bed = os.path.join(outdir, "merged.bed")

    # only 1 bed = nothing changes
    ananse.utils.bed_merge([unsorted_bed], merged_bed)
    assert filecmp.cmp(merged_bed, sorted_bed)

    # 2 bed, same content
    ananse.utils.bed_merge([unsorted_bed, sorted_bed], merged_bed)
    assert filecmp.cmp(merged_bed, sorted_bed)

    # 2 beds, different content
    ananse.utils.bed_merge([unsorted_bed, second_bed], merged_bed)
    assert not filecmp.cmp(merged_bed, sorted_bed)
    assert not filecmp.cmp(merged_bed, unsorted_bed)
    with open(merged_bed) as mb:
        assert len(mb.readlines()) == 3


# def test_cleanpath():
#     path = "."
#     expected = __file__
#     res = ananse.utils.cleanpath(path)
#     assert res == expected
