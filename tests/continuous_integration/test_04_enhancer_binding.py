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


genome = os.path.join(outdir, "genome.fa")
write_file(genome, [">chr1", "N" * 50000])

bed1 = os.path.join(outdir, "bed1.bed")
write_file(bed1, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n"])

bed2 = os.path.join(outdir, "bed2.bed")
write_file(bed2, ["chr1\t4000\t5000\n", "chr1\t2000\t3000\n"])


h0 = "@HD	VN:1.6	SO:coordinate"
h1 = "@SQ	SN:chr1	LN:50000"
line1 = (
    "read1	147	chr1	10003	40	11S90M	=	10048	-46	"
    + "CCCTACCCTCTCCCTATCCCTAACCCTAACCCCAACCCTAACCCTATCCCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA	"
    + "A77--7-7---7-7---A77---AA7----<7-AAJAA-7JJFF<--F-A-AFFFF<FJJJJF-AFJF7F-JJJFJFFFJFF<FJJJJFJJFJJFFFFFAA	"
)
line2 = (
    "read2	83	chr1	10004	30	2S45M1D54M	=	10074	-30	"
    + "ATCCCTAACCCTAACCCTAACCCTAACCCTACCCCTACCCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT	"
    + "--JAA7F-FAFA-7JJFA--F<7-FF<<FAF7<7F7A-FFAF7-FJJJFJJ----J<JFA-JAF7JFJFJF<<JFJF<JJJFFJJJAAAA-JFFFA-FAA-	"
)
line3 = (
    "read3	163	chr1	10027	40	100M	=	10032	105	"
    + "ACCCGAACCCTAACCCTAACCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCGAACCCA	"
    + "AAFFFJJJJJJJJJJJFJJJFJJJFJFJJFJJJJ<-FJJFJAFFJA7AFAJJJJFJFJ-<F-AAJJ<FF7-J-AAJ--<JJJ--AAJ-77-AA-7A<-A-	"
)

unsorted_bam = os.path.join(outdir, "unsorted.bam")
tmp_sam = os.path.join(outdir, "tmp.sam")
write_file(tmp_sam, [h0, h1, line2, line1])
pysam.view(tmp_sam, "-b", "-o", unsorted_bam, catch_stdout=False)

sorted_bam = os.path.join(outdir, "sorted.bam")
write_file(tmp_sam, [h0, h1, line1, line2])
pysam.view(tmp_sam, "-b", "-o", sorted_bam, catch_stdout=False)

genomepy.utils.rm_rf(tmp_sam)


# all expected outputs
combined_bed = os.path.join(outdir, "combined.bed")
combined_bam = os.path.join(outdir, "combined.bam")
# coverage_file = os.path.join(outdir, "raw_scoredpeaks.bed")
scored_peaks = os.path.join(outdir, "scoredpeaks.bed")
raw_motif_scores = os.path.join(outdir, "raw_scoredmotifs.bed")
scored_motifs = os.path.join(outdir, "scoredmotifs.bed")
outfile = os.path.join(outdir, "binding.tsv")


def test_is_narrowpeak():
    np = os.path.join(outdir, "f.narrowPeak")
    write_file(np, ["1	629812	630105	narrowPeak1	6047	.	0	0	0	122"])

    cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=[])
    assert cbed.is_narrowpeak(np) is True
    assert cbed.is_narrowpeak(bed1) is False


def test_bed_resize():
    cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=[])
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


def test_cbedf():
    cbed = ananse.enhancer_binding.CombineBedFiles(
        genome=genome, peakfiles=[bed1, bed2]
    )
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


def test_cbamf():
    cbam = ananse.enhancer_binding.CombineBamFiles(bams=[unsorted_bam, sorted_bam])
    cbam.run(outfile=combined_bam, force=True)
    assert os.path.exists(combined_bam)
    assert os.path.exists(f"{combined_bam}.bai")


def test_normalize_peaks():
    sp = ananse.enhancer_binding.ScorePeaks(bed=None, bam=None)
    bam_cov = os.path.join(outdir, "bam_cov.bed")
    write_file(bam_cov, ["chr1\t0\t200\t10", "chr1\t0\t200\t20", "chr1\t0\t200\t30"])
    norm_cov = os.path.join(outdir, "norm_cov.bed")
    sp.normalize_peaks(bam_cov, norm_cov, dist_func="scale_dist")

    scores = []
    norm_scores = []
    with open(norm_cov) as bed:
        for line in bed:
            line = line.strip().split()
            scores.append(line[1])
            norm_scores.append(line[2])

    assert len(scores) == 3 + 1  # lines + header
    assert scores[1:] == ["10", "20", "30"]
    assert norm_scores[1:] == ["0.0", "0.5", "1.0"]


def test_sp():
    sp = ananse.enhancer_binding.ScorePeaks(bed=combined_bed, bam=combined_bam)
    sp.run(outfile=scored_peaks, force=True)

    with open(scored_peaks) as f:
        lines = f.readlines()
    assert len(lines[0].split()) == 6


def test_get_motif_scores():
    # TODO: cant get gimmemotifs to skip making a GC background index
    # fake_cg_index = "~/.cache/gimmemotifs/genome.fa.gcfreq.100.feather"
    # try:
    #     import pandas as pd
    #     import numpy as np
    #     df = pd.DataFrame({
    #         "chrom": ["chr1"],
    #         "start": ["0"],
    #         "end": ["100"],
    #         "w100": ["0.0"],
    #         "n100": ["0.0"],
    #         "w200": [np.NaN],
    #         "n200": [np.NaN],
    #         "w500": [np.NaN],
    #         "n500": [np.NaN],
    #     })
    #     df.to_feather(fake_cg_index)
    #
    #     pfmfile = os.path.join(test_dir, "example_data", "debug.pfm")
    #     sm = ananse.enhancer_binding.ScoreMotifs(genome, combined_bed, pfmfile=pfmfile)
    #     sm.get_motif_scores(combined_bed, raw_motif_scores)
    # finally:
    #     genomepy.utils.rm_rf(fake_cg_index)
    #
    #     write_file(
    #         raw_motif_scores,
    #         [
    #             "motif	region	zscore",
    #             "GM.5.0.Sox.0001	chr1:400-600	-0.5444524936254616",
    #             "GM.5.0.Homeodomain.0001	chr1:2400-2600	-0.3774763844954927",
    #             "GM.5.0.Sox.0001	chr1:4400-4600	-0.5444524936254616",
    #         ],
    #     )

    write_file(
        raw_motif_scores,
        [
            "motif	region	zscore",
            "GM.5.0.Sox.0001	chr1:400-600	-0.5444524936254616",
            "GM.5.0.Homeodomain.0001	chr1:2400-2600	-0.3774763844954927",
            "GM.5.0.Sox.0001	chr1:4400-4600	-0.5444524936254616",
        ],
    )


def test_normalize_motifs():
    sm = ananse.enhancer_binding.ScoreMotifs(genome, combined_bed)
    sm.normalize_motifs(raw_motif_scores, scored_motifs)

    with open(raw_motif_scores) as f:
        lines1 = f.readlines()
    with open(scored_motifs) as f:
        lines2 = f.readlines()

    assert len(lines2[0].split()) == len(lines1[0].split()) + 1


def test_sm():
    pfmfile = os.path.join(test_dir, "data", "debug.pfm")
    sm = ananse.enhancer_binding.ScoreMotifs(genome, combined_bed, pfmfile=pfmfile)
    sm.run(outfile=scored_motifs, force=False)


def test_filter_transcription_factors():
    pfmfile = os.path.join(test_dir, "data", "debug.pfm")
    b = ananse.enhancer_binding.Binding(None, None, pfmfile=pfmfile)

    # curation filter
    m2f = b.filter_transcription_factors(curation_filter=None)
    assert m2f.shape[0] == 9  # all TFs in the file
    m2f = b.filter_transcription_factors(curation_filter=True)
    assert m2f.shape[0] == 8  # all curated TFs
    m2f = b.filter_transcription_factors(curation_filter=False)
    assert m2f.shape[0] == 1  # all non-curated TFs

    # tf filter
    tf_list = os.path.join(outdir, "tf_list.txt")
    write_file(tf_list, ["SOX12"])
    m2f = b.filter_transcription_factors(tf_list=tf_list, whitelist=True)
    assert m2f.shape[0] == 1
    m2f = b.filter_transcription_factors(tf_list=tf_list, whitelist=False)
    assert m2f.shape[0] == 8


def test_get_binding_score():
    pfmfile = os.path.join(test_dir, "data", "debug.pfm")
    b = ananse.enhancer_binding.Binding(None, None, pfmfile=pfmfile)
    b.get_binding_score(scored_motifs, scored_peaks, outfile)

    assert os.path.exists(outfile)
