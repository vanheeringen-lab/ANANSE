import ananse.enhancer_binding
import os


genome = "/home/siebrenf/genomes/hg38/hg38.fa"
list_of_enhancer_region_peakfiles = ["testdata/hg38-fibroblast_H3K27ac_peaks.broadPeak",
                                     "testdata/hg38-keratinocyte_H3K27ac_peaks.broadPeak"]
outfile = "testdata/CombinePeaks.out.bed"
if os.path.exists(outfile):
    os.remove(outfile)

cpf = ananse.enhancer_binding.CombinePeakFiles(genome, list_of_enhancer_region_peakfiles, outfile)
cpf.run()


peaks = outfile
# bams = "testdata/hg38-fibroblast_H3K27ac_rep1.samtools-coordinate.bam"
bams = ["testdata/hg38-fibroblast_H3K27ac_rep1.samtools-coordinate.bam",
        "testdata/hg38-keratinocyte_H3K27ac_rep1.samtools-coordinate.bam"]
outfile = "testdata/ScorePeaks.out.bed"
if os.path.exists(outfile):
    os.remove(outfile)

sp = ananse.enhancer_binding.ScorePeaks(bams, peaks, outfile)
sp.run()


scored_peaks = outfile
outfile = "testdata/table.txt"
b = ananse.enhancer_binding.Binding(genome, scored_peaks, outfile)
b.run()
