# import os
#
# import genomepy.utils
#
# import ananse.enhancer_binding
# # from tests.benchmark import distplot
#
#
# # prep
# test_dir = os.path.dirname(__file__)
# data_dir = os.path.join(test_dir, "data")
# genomepy.utils.mkdir_p(data_dir)
# outdir = os.path.join(test_dir, "output")
# genomepy.utils.mkdir_p(outdir)
# examples_dir = os.path.join(test_dir, "example_data")
#
#
# # H3K27Ac data
# genome = os.path.join(data_dir, "hg38.fa")
# peakfiles = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_peaks.broadPeak")
# bams = os.path.join(data_dir, "hg38-keratinocyte_H3K27ac_rep1.samtools-coordinate.bam")
#
#
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
