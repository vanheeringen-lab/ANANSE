from loguru import logger

from ananse.peakpredictor import predict_peaks


@logger.catch
def binding(args):
    # all arguments are checked in the function
    predict_peaks(
        args.outdir,
        atac_bams=args.atac_bams,
        histone_bams=args.histone_bams,
        cage_tpms=args.cage_tpms,
        columns=args.columns,
        regions=args.regions,
        reference=args.reference,
        factors=args.tfs,
        genome=args.genome,
        pfmfile=args.pfmfile,
        pfmscorefile=args.pfmscorefile,
        jaccard_cutoff=args.jaccard_cutoff,
        ncore=args.ncore,
    )
# example
