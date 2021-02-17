import os

import genomepy.utils
import pytest

import ananse.enhancer_binding
from ananse.commands.enhancer_binding import run_binding
import ananse.utils
from .test_02_utils import write_file, write_bam, h0, h1, line1, line2, line3

# prep

test_dir = os.path.dirname(os.path.dirname(__file__))
outdir = os.path.join(test_dir, "output")
genomepy.utils.mkdir_p(outdir)

# beds

genome = os.path.join(outdir, "genome.fa")
write_file(genome, [">chr1", "N" * 50000])

bed1 = os.path.join(outdir, "bed1.bed")
write_file(bed1, ["chr1\t0\t1000\n", "chr1\t2000\t3000\n"])

bed2 = os.path.join(outdir, "bed2.bed")
write_file(bed2, ["chr1\t4000\t5000\n", "chr1\t2000\t3000\n"])

sp_bed_input = os.path.join(outdir, "sp_input.bed")
write_file(sp_bed_input, ["chr1\t10003\t10203\n", "chr1\t10203\t10403\n"])

# bams

bam1 = os.path.join(outdir, "bam1.bam")
write_bam(bam1, [h0, h1, line1, line2, line2])
ananse.utils.bam_index(bam1)

bam2 = os.path.join(outdir, "bam2.bam")
write_bam(bam2, [h0, h1, line1, line3, line3])
ananse.utils.bam_index(bam2)

# shared in/outputs

combined_bed = os.path.join(outdir, "combined.bed")
raw_peak_scores = os.path.join(outdir, "raw_scoredpeaks.bed")
scored_peaks = os.path.join(outdir, "scoredpeaks.bed")
raw_motif_scores = os.path.join(outdir, "raw_scoredmotifs.bed")
scored_motifs = os.path.join(outdir, "scoredmotifs.bed")
outfile = os.path.join(outdir, "binding.tsv")


def test_is_narrowpeak():
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

    cbed = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=[])
    assert cbed.is_narrowpeak(np) is True
    assert cbed.is_narrowpeak(bp) is False


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


def test_peaks_count():
    sp = ananse.enhancer_binding.ScorePeaks(
        bams=[bam1, bam2], bed=sp_bed_input, ncore=1, verbose=True
    )
    sp.peaks_count(outdir)

    files = os.listdir(outdir)
    counts = [f for f in files if f.endswith(".regions.bed")]
    assert len(counts) == 2


def test_peaks_merge():
    sp = ananse.enhancer_binding.ScorePeaks(bams=[], bed=None, ncore=1, verbose=True)
    sp.peaks_merge(outdir, raw_peak_scores, sp.ncore)

    with open(raw_peak_scores) as f:
        content = f.readlines()[0]

    assert len(content.strip().split("\t")) == 4


def test_normalize_peaks():
    sp = ananse.enhancer_binding.ScorePeaks(bams=[], bed=None, ncore=1, verbose=True)
    raw_cov = os.path.join(outdir, "raw_cov.bed")
    write_file(raw_cov, ["chr1\t0\t200\t10", "chr1\t0\t200\t20", "chr1\t0\t200\t30"])
    norm_cov = os.path.join(outdir, "norm_cov.bed")
    sp.peaks_fit(raw_cov, norm_cov, dist_func="scale_dist")

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
    sp = ananse.enhancer_binding.ScorePeaks(
        bams=[bam1, bam2], bed=sp_bed_input, ncore=1, verbose=True
    )
    sp.run(outfile=scored_peaks, dist_func="scale_dist", force=True)

    with open(scored_peaks) as f:
        lines = f.readlines()
    peak1 = lines[1].split()
    assert len(peak1) == 6
    assert peak1[1] == "0.375"  # raw score
    assert peak1[2] == "1.0"  # norm score (scaled)


def test_motifs_get_scores():
    # TODO: cant get gimmemotifs to skip making a GC background index
    # fake_cg_index = "~/.cache/gimmemotifs/genome.fa.gcfreq.100.feather"
    # try:
    #     import pandas as pd
    #     import numpy as np
    #     df = pd.DataFrame({
    #         "chrom": ["chr1"], "start": ["0"], "end": ["100"],
    #         "w100": ["0.0"], "n100": ["0.0"], "w200": [np.NaN],
    #         "n200": [np.NaN], "w500": [np.NaN], "n500": [np.NaN],
    #     })
    #     df.to_feather(fake_cg_index)
    #
    #     pfmfile = os.path.join(test_dir, "example_data", "debug.pfm")
    #     sm = ananse.enhancer_binding.ScoreMotifs(genome, combined_bed, pfmfile=pfmfile)
    #     sm.get_motif_scores(combined_bed, raw_motif_scores)
    # finally:
    #     genomepy.utils.rm_rf(fake_cg_index)

    write_file(
        raw_motif_scores,
        [
            "motif	region	zscore",
            "GM.5.0.Sox.0001	chr1:400-600	-0.5444524936254616",
            "GM.5.0.Homeodomain.0001	chr1:2400-2600	-0.3774763844954927",
            "GM.5.0.Sox.0001	chr1:10003-10203	-0.5444524936254616",
        ],
    )
    pass


def test_normalize_motifs():
    sm = ananse.enhancer_binding.ScoreMotifs(genome, combined_bed)
    sm.motifs_normalize(raw_motif_scores, scored_motifs)

    with open(raw_motif_scores) as f:
        lines1 = f.readlines()
    with open(scored_motifs) as f:
        lines2 = f.readlines()

    assert len(lines2[0].split()) == len(lines1[0].split()) + 1


def test_sm():
    # TODO: fix when get_motif_score can be run with toy data
    # pfmfile = os.path.join(test_dir, "data", "debug.pfm")
    # sm = ananse.enhancer_binding.ScoreMotifs(genome, combined_bed, pfmfile=pfmfile)
    # sm.run(outfile=scored_motifs, force=True)
    pass


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


def test_run_binding():
    # test API wrapper
    run_binding(genome=genome, bams=[bam1], peakfiles=bed1, outdir=outdir, force=False)

    with pytest.raises(FileNotFoundError):
        run_binding(
            genome="/not/a/real/genome.fa",
            bams=[bam1],
            peakfiles=bed1,
            outdir=outdir,
            force=False,
        )
