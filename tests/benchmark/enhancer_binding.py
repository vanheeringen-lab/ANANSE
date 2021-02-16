import os

import genomepy.utils

from ananse.enhancer_binding import (
    CombineBedFiles,
    CombineBamFiles,
    ScorePeaks,
    ScoreMotifs,
    Binding,
)
from tests.benchmark.utils import distplot


# prep

# test_dir = os.path.join(os.getcwd(), "tests")
test_dir = os.path.dirname(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
genomepy.utils.mkdir_p(data_dir)
outdir = os.path.join(test_dir, "output")
genomepy.utils.mkdir_p(outdir)
intermediate_dir = os.path.join(outdir, "intermediate_results")
genomepy.utils.mkdir_p(intermediate_dir)
ncore = max(1, os.cpu_count() - 2)


# H3K27Ac data
genome = os.path.join(data_dir, "hg38.fa")
peakfiles = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_peaks.broadPeak")
bams = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_rep1.samtools-coordinate.bam")


# download test data locally
for file in [genome, peakfiles, bams]:
    if not os.path.exists(file):
        url = "https://mbdata.science.ru.nl/ANANSE/tests/data/" + file
        genomepy.utils.download_file(url, file)


# group 1 (can run simultaneously)
cbed = CombineBedFiles(genome=genome, peakfiles=peakfiles, verbose=True)
combined_bed = os.path.join(intermediate_dir, "combined.bed")
cbed.run(outfile=combined_bed, width=200, force=False)

cbam = CombineBamFiles(bams=bams, ncore=ncore, verbose=True)
combined_bam = os.path.join(intermediate_dir, "combined.bam")
cbam.run(outfile=combined_bam, force=False)

# group 2 (can run when input is ready)
sp = ScorePeaks(bed=combined_bed, bam=combined_bam, ncore=ncore, verbose=True)

# benchmark peak normalization
for func, kwargs in zip(
    [
        "peak_rank_file_dist",  # Quan's file
        "peak_rank_dist",  # scalable version of Quan's file
        "scale_dist",  # no normalization
        "log_scale_dist",  # no normalization, but with a log
        "scipy_dist",  # see kwargs
        "scipy_dist",  # see kwargs
    ],
    [{}, {}, {}, {}, {"dist": "loglaplace"}, {"dist": "lognorm"}],  # fits best
):
    suffix = kwargs.get("dist")
    if not suffix:
        suffix = func
    scored_peaks = os.path.join(intermediate_dir, f"scoredpeaks_{suffix}.bed")
    sp.run(outfile=scored_peaks, dist_func=func, force=False, **kwargs)
    distplot(scored_peaks, score_col=5)

run_gimme = False
if run_gimme:
    scored_peaks = os.path.join(intermediate_dir, "scoredpeaks.bed")
    sp = ScorePeaks(bed=combined_bed, bam=combined_bam, ncore=ncore, verbose=True)
    sp.run(outfile=scored_peaks, dist_func="peak_rank_dist", force=False)

    sm = ScoreMotifs(
        genome=genome, bed=combined_bed, pfmfile=None, ncore=ncore, verbose=True
    )
    scored_motifs = os.path.join(intermediate_dir, "scoredmotifs.bed")
    sm.run(outfile=scored_motifs, force=True)

    # group 3 (end result)
    b = Binding(
        peak_weights=scored_peaks,
        motif_weights=scored_motifs,
        pfmfile=None,
        model=None,
        curation_filter=None,
        tf_list=None,
        whitelist=True,
        ncore=ncore,
        verbose=True,
    )
    outfile = os.path.join(outdir, "binding.tsv")
    b.run(outfile=outfile, force=True)
