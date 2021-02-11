import os

import genomepy.utils

import ananse.enhancer_binding
from tests.benchmark import distplot


# prep
test_dir = os.path.dirname(__file__)
data_dir = os.path.join(test_dir, "data")
genomepy.utils.mkdir_p(data_dir)
out_dir = os.path.join(test_dir, "output")
genomepy.utils.mkdir_p(out_dir)


# H3K27Ac data
genome = os.path.join(data_dir, "hg38.fa")
peaks_file = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_peaks.broadPeak")
bam_file = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_rep1.samtools-coordinate.bam")

# # ATAC-seq data
# genome = ""
# peaks_file = ""
# bam_file = ""
#
# # p300 data
# genome = ""
# peaks_file = ""
# bam_file = ""


# download test data locally
for file in [genome, peaks_file, bam_file]:
    if not os.path.exists(file):
        url = "https://mbdata.science.ru.nl/ANANSE/tests/data/" + file
        genomepy.utils.download_file(url, file)


# group 1 (can run simultaneously)
print("Combining bed files")
cbedf = ananse.enhancer_binding.CombineBedFiles(genome=genome, peakfiles=peaks_file)
combined_bed = os.path.join(out_dir, "combined.bed")
cbedf.run(outfile=combined_bed, width=200)

print("Combining bam files")
cbamf = ananse.enhancer_binding.CombineBamFiles(bams=bam_file)
combined_bam = os.path.join(out_dir, "combined.bam")
cbamf.run(outfile=combined_bam)

# group 2 (can run when input is ready)
print("Scoring peaks")
sp = ananse.enhancer_binding.ScorePeaks(bed=combined_bed, bam=combined_bam)
scored_peaks = os.path.join(out_dir, "scoredpeaks_peakrankfile.bed")
sp.run(outfile=scored_peaks, dist_func="peak_rank_file_dist", **{"dist": "loglaplace"})
distplot(scored_peaks)

print("Scoring motifs")
sm = ananse.enhancer_binding.ScoreMotifs(genome=genome, bed=combined_bed, ncore=min(os.cpu_count() - 2, 1))
scored_motifs = os.path.join(out_dir, "scoredmotifs.bed")
sm.run(outfile=scored_motifs)

# group 3 (end result)
print("Predict TF binding")
b = ananse.enhancer_binding.Binding(
    peak_weights=scored_peaks,
    motif_weights=scored_motifs,
    ncore=min(os.cpu_count() - 2, 1)
)
outfile = os.path.join(out_dir, "table.txt")
b.run(outfile=outfile)



# # run enhancer_binding.py
# cpf = ananse.enhancer_binding.CombinePeakFiles(genome, peaks_file)
# outfile = "tests/output/CombinePeaks.out.bed"
# cpf.run(outfile)
#
#
# peaks = "tests/output/CombinePeaks.out.bed"
# sp = ananse.enhancer_binding.ScorePeaks(bam_file, peaks)
# # sp.run(outfile="tests/output/ScorePeaks_scale.out.bed", dist_func="scale_dist")
# # sp.run(outfile="tests/output/ScorePeaks_logscale.out.bed", dist_func="log_scale_dist")
# # sp.run(outfile="tests/output/ScorePeaks_lognorm.out.bed", dist_func="scipy_dist", **{"dist": "lognorm"})
# # sp.run(outfile="tests/output/ScorePeaks_loglaplace.out.bed", dist_func="scipy_dist", **{"dist": "loglaplace"})
# # sp.run(outfile="tests/output/ScorePeaks_peakrank.out.bed", dist_func="peak_rank_dist")
# sp.run(outfile="tests/output/ScorePeaks_peakrankfile.out.bed", dist_func="peak_rank_file_dist")
#
# # distplot("tests/output/ScorePeaks_scale.out.bed")
# # distplot("tests/output/ScorePeaks_logscale.out.bed")
# # distplot("tests/output/ScorePeaks_lognorm.out.bed")
# # distplot("tests/output/ScorePeaks_loglaplace.out.bed")
# # distplot("tests/output/ScorePeaks_peakrank.out.bed")
# distplot("tests/output/ScorePeaks_peakrankfile.out.bed")
#
# scored_peaks = "tests/output/ScorePeaks_peakrankfile.out.bed"  # from sp
# b = ananse.enhancer_binding.Binding(genome, scored_peaks, ncore=min(os.cpu_count() - 2, 1))
#
# # gms_out = "tests/output/get_motif_scores.out"
# # b.get_motif_scores(scored_peaks, gms_out)
# outfile = "tests/output/table.txt"
# b.run(outfile)
