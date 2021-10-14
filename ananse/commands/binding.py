from loguru import logger

from ananse.peakpredictor import predict_peaks
from ananse.utils import check_path, check_input_factors


@logger.catch
def binding(args):
    predict_peaks(
        check_path(args.outdir, error_missing=False),
        atac_bams=check_path(args.atac_bams),
        histone_bams=check_path(args.histone_bams),
        regionfiles=check_path(args.regionfiles),
        reference=check_path(args.reference),
        factors=check_input_factors(args.factors),
        genome=args.genome,  # checked in CLI
        pfmfile=check_path(args.pfmfile),
        pfmscorefile=check_path(args.pfmscorefile),
        jaccard_cutoff=args.jaccard_cutoff,
        ncore=args.ncore,
    )
