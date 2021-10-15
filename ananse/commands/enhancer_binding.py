# import os
#
# import genomepy.utils
# from loguru import logger
#
# from ananse.enhancer_binding import (
#     CombineBedFiles,
#     ScorePeaks,
#     ScoreMotifs,
#     Binding,
# )
# from ananse.utils import clean_tmp
#
#
# @logger.catch
# def run_binding(
#     genome,
#     peakfiles,
#     bams,
#     outdir,
#     peak_width=200,
#     dist_func="peak_rank_file_dist",
#     pfmfile=None,
#     curation_filter=None,
#     tf_list=None,
#     whitelist=True,
#     model=None,
#     ncore=1,
#     force=False,
#     keep_intermediates=True,
#     verbose=True,
#     **kwargs,
# ):
#     """
#     Predict transcription factor binding in specified regions
#
#     Args:
#         genome: path to the genome fasta used to align the bams and peaks to
#         peakfiles: one or more BED format files with putative enhancer regions (e.g. narrowPeak, broadPeak)
#         bams: one or more BAM format files where reads mark enhancer activity (H3K27Ac/p300 ChIP-seq or ATAC-seq)
#         outdir: directory where you wish to store the output
#         peak_width: peakfiles are resized to this width (default 200 bp)
#         dist_func: bam reads are normalized to the selected distribution (default: an empirical distribution)
#         pfmfile: the pfm file of the transcription factors to search for (default gimme.vertebrate.v5.0)
#         curation_filter: True = curated TFs, False = no curated TFs, None = all TFs (default: None)
#         tf_list: optional file with single column TF names
#         whitelist: True = use tf_list as a whitelist. False = use tf_list as a blacklist
#         model: classification model to use (default: dream)
#         ncore: number of cores to use
#         force: overwrite earlier intermediate data? (default: False)
#         keep_intermediates: keep intermediate data after completion? (default: True)
#         verbose: keep you informed of the progress? (default: True)
#         **kwargs: passed to the selected dist_func
#
#     Returns:
#         binding.tsv: the strongest transcription factor and its binding score for each region in the peakfile(s)
#     """
#     # clean up previous ANANSE tmp files
#     clean_tmp()
#
#     # check input file paths
#     files = []
#     for arg in [genome, peakfiles, bams, pfmfile, tf_list, model]:
#         if arg:
#             if isinstance(arg, list):
#                 files.extend(arg)
#             else:
#                 files.append(arg)
#     for file in files:
#         if not os.path.exists(file):
#             logger.exception(f"Could not find {file}!")
#             exit(1)
#
#     outfile = os.path.join(outdir, "binding.tsv")
#     intermediate_dir = os.path.join(outdir, "intermediate_results")
#     if force or not os.path.exists(outfile):
#         genomepy.utils.mkdir_p(intermediate_dir)
#
#         cbed = CombineBedFiles(genome=genome, peakfiles=peakfiles, verbose=verbose)
#         combined_bed = os.path.join(intermediate_dir, "combined.bed")
#         cbed.run(outfile=combined_bed, width=peak_width, force=force)
#
#         sp = ScorePeaks(bams=bams, bed=combined_bed, ncore=ncore, verbose=verbose)
#         scored_peaks = os.path.join(intermediate_dir, "scoredpeaks.bed")
#         sp.run(outfile=scored_peaks, dist_func=dist_func, force=force, **kwargs)
#
#         sm = ScoreMotifs(
#             genome=genome,
#             bed=scored_peaks,
#             pfmfile=pfmfile,
#             ncore=ncore,
#             verbose=verbose,
#         )
#         scored_motifs = os.path.join(intermediate_dir, "scoredmotifs.bed")
#         sm.run(outfile=scored_motifs, force=force)
#
#         b = Binding(
#             peak_weights=scored_peaks,
#             motif_weights=scored_motifs,
#             pfmfile=pfmfile,
#             model=model,
#             curation_filter=curation_filter,
#             tf_list=tf_list,
#             whitelist=whitelist,
#             ncore=ncore,
#             verbose=verbose,
#         )
#         b.run(outfile=outfile, force=force)
#
#     if not keep_intermediates:
#         genomepy.utils.rm_rf(intermediate_dir)
#     if verbose:
#         logger.info("ANANSE binding finished successfully!")
