
from __future__ import print_function
import sys
import os

import ananse.quantile


def quantile(args):
    # if not os.path.exists(args.fin_rpkm):
    #     print("File %s does not exist!" % args.fin_rpkm)
    #     sys.exit(1)

    b = ananse.quantile.Quantile(
        bed_input=args.bed_input,
        bed_output=args.bed_output
    )
    b.run_quantile(
        args.bed_input, args.bed_output
    )
