import os

import genomepy.utils
import pysam

import ananse.enhancer_binding


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


genome = os.path.join(outdir, "genome.fa")
write_file(genome, [">chr1", "N"*50000])

bed1 = os.path.join(outdir, "bed1.bed")
write_file(bed1, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n"])

bed2 = os.path.join(outdir, "bed2.bed")
write_file(bed2, ["chr1\t4000\t5000\n", "chr1\t2000\t3000\n"])


def test_is_narrowpeak():
    np = os.path.join(outdir, "f.narrowPeak")
    write_file(np, ["1	629812	630105	narrowPeak1	6047	.	0	0	0	122"])

    cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=[])
    assert cbed.is_narrowpeak(np) is True
    assert cbed.is_narrowpeak(bed1) is False


def test_bed_resize():
    cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=[])

    bed_in = os.path.join(outdir, "bed_in.bed")
    write_file(bed_in, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n"])
    bed_out = os.path.join(outdir, "bed_out.bed")

    width = 200
    cbed.bed_resize(genome, bed_in, bed_out, width=width)
    with open(bed_out) as f:
        lines = f.readlines()

    chrom, start, stop = lines[0].split()[0:3]
    assert int(stop) - int(start) == width


def test_cbf():
    cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=[bed1, bed2])
    combined_bed = os.path.join(outdir, "combined.bed")
    width = 200
    cbed.run(outfile=combined_bed, width=width, force=True)

    with open(combined_bed) as f:
        lines = f.readlines()

    # 3 unique regions over the 2 bed files
    assert len(lines) == 3

    for line in lines:
        chrom, start, stop = line.split()[0:3]
        assert int(stop) - int(start) == width




# test_dir = os.path.dirname(os.path.dirname(__file__))
# data_dir = os.path.join(test_dir, "data")
# genomepy.utils.mkdir_p(data_dir)
# outdir = os.path.join(test_dir, "output")
# genomepy.utils.mkdir_p(outdir)
# examples_dir = os.path.join(test_dir, "example_data")


# # H3K27Ac data
# genome = os.path.join(data_dir, "hg38.fa")
# peakfiles = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_peaks.broadPeak")
# bams = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_rep1.samtools-coordinate.bam")
#
# assert "" == bams

# # download test data locally
# for file in [genome, peakfiles, bams]:
#     if not os.path.exists(file):
#         url = "https://mbdata.science.ru.nl/ANANSE/tests/data/" + file
#         genomepy.utils.download_file(url, file)
#
#
# # group 1 (can run simultaneously)
# cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=peakfiles)
# combined_bed = os.path.join(outdir, "combined.bed")
# cbed.run(outfile=combined_bed, width=200, force=True)
#
# cbam = ananse.enhancer_binding.CombineBamFiles(bams=bams)
# combined_bam = os.path.join(outdir, "combined.bam")
# cbam.run(outfile=combined_bam, force=True)
#
# # group 2 (can run when input is ready)
# sp = ananse.enhancer_binding.ScorePeaks(bed=combined_bed, bam=combined_bam)
# scored_peaks = os.path.join(outdir, "scoredpeaks.bed")
# sp.run(
#     outfile=scored_peaks,
#     dist_func="peak_rank_file_dist",
#     **{"dist": "loglaplace"},
#     force=True
# )
# # distplot(scored_peaks)
#
# sm = ananse.enhancer_binding.ScoreMotifs(
#     genome=genome, bed=combined_bed, ncore=max(os.cpu_count() - 2, 1)
# )
# scored_motifs = os.path.join(outdir, "scoredmotifs.bed")
# sm.run(outfile=scored_motifs, force=True)
#
# # group 3 (end result)
# b = ananse.enhancer_binding.Binding(
#     peak_weights=scored_peaks,
#     motif_weights=scored_motifs,
#     ncore=max(os.cpu_count() - 2, 1),
# )
# outfile = os.path.join(outdir, "binding.tsv")
# b.run(outfile=outfile, force=True)
