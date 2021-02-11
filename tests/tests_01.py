import os
import genomepy.utils
import ananse.enhancer_binding
from tests.benchmark import distplot


# H3K27Ac
genome = "tests/data/hg38.fa"
peaks_file = "tests/data/hg38-keratinocyte_H3K27ac_peaks.broadPeak"
bam_file = "tests/data/hg38-keratinocyte_H3K27ac_rep1.samtools-coordinate.bam"

# # ATAC-seq
# genome = "tests/data/GRCh38.p13.fa"
# peaks_file = "tests/data/GRCh38.p13-GSM2264802_peaks.narrowPeak"
# bam_file = "tests/data/GRCh38.p13-GSM2264802.samtools-coordinate.bam"
#
# # p300
# genome = ""
# peaks_file = ""
# bam_file = ""

# download test data locally
for file in [genome, peaks_file, bam_file]:
    if not os.path.exists(file):
        url = "https://mbdata.science.ru.nl/ANANSE/" + file  # TODO upload data there
        genomepy.utils.download_file(url, file)

# genomepy.utils.rm_rf("tests/output")
genomepy.utils.mkdir_p("tests/output")

# run enhancer_binding.py
cpf = ananse.enhancer_binding.CombinePeakFiles(genome, peaks_file)
outfile = "tests/output/CombinePeaks.out.bed"
cpf.run(outfile)


peaks = "tests/output/CombinePeaks.out.bed"
sp = ananse.enhancer_binding.ScorePeaks(bam_file, peaks)
# sp.run(outfile="tests/output/ScorePeaks_scale.out.bed", dist_func="scale_dist")
# sp.run(outfile="tests/output/ScorePeaks_logscale.out.bed", dist_func="log_scale_dist")
# sp.run(outfile="tests/output/ScorePeaks_lognorm.out.bed", dist_func="scipy_dist", **{"dist": "lognorm"})
# sp.run(outfile="tests/output/ScorePeaks_loglaplace.out.bed", dist_func="scipy_dist", **{"dist": "loglaplace"})
# sp.run(outfile="tests/output/ScorePeaks_peakrank.out.bed", dist_func="peak_rank_dist")
sp.run(outfile="tests/output/ScorePeaks_peakrankfile.out.bed", dist_func="peak_rank_file_dist")

# distplot("tests/output/ScorePeaks_scale.out.bed")
# distplot("tests/output/ScorePeaks_logscale.out.bed")
# distplot("tests/output/ScorePeaks_lognorm.out.bed")
# distplot("tests/output/ScorePeaks_loglaplace.out.bed")
# distplot("tests/output/ScorePeaks_peakrank.out.bed")
distplot("tests/output/ScorePeaks_peakrankfile.out.bed")

scored_peaks = "tests/output/ScorePeaks_peakrankfile.out.bed"  # from sp
b = ananse.enhancer_binding.Binding(genome, scored_peaks, ncore=min(os.cpu_count() - 2, 1))

# gms_out = "tests/output/get_motif_scores.out"
# b.get_motif_scores(scored_peaks, gms_out)
outfile = "tests/output/table.txt"
b.run(outfile)
