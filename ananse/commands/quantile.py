
from __future__ import print_function
import sys
import os

import ananse.quantile


def quantile(args):
    # if not os.path.exists(args.fin_rpkm):
    #     print("File %s does not exist!" % args.fin_rpkm)
    #     sys.exit(1)

    b = ananse.quantile.Quantile(
        bam_input=args.bam_input,
        broadPeak=args.broadPeak,
        bed_output=args.bed_output
    )
    b.run_quantile(
        args.bam_input, args.broadPeak, args.bed_output
    )
