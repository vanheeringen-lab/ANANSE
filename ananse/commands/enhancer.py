
from __future__ import print_function
import sys
import os

import ananse.enhancer


def enhancer(args):
    # if not os.path.exists(args.fin_rpkm):
    #     print("File %s does not exist!" % args.fin_rpkm)
    #     sys.exit(1)
    genome=args.genome
    etype=args.etype

    if genome == "hg38" and etype == "H3K27ac":
        b = ananse.enhancer.Enhancer(
            genome=args.genome, 
            bam_input=args.bam_input,
            epeak=args.epeak,
            bed_output=args.bed_output
        )
        b.run_enhancer(
            args.bam_input, args.epeak, args.bed_output
        )
    elif etype == "p300":
        b = ananse.enhancer.P300Enhancer(
            genome=args.genome, 
            bam_input=args.bam_input,
            epeak=args.epeak,
            bed_output=args.bed_output
        )
        b.run_enhancer(
            args.bam_input, args.epeak, args.bed_output
        )
    elif etype == "ATAC":
        b = ananse.enhancer.AtacEnhancer(
            genome=args.genome, 
            bam_input=args.bam_input,
            epeak=args.epeak,
            bed_output=args.bed_output
        )
        b.run_enhancer(
            args.bam_input, args.epeak, args.bed_output
        )

       
        