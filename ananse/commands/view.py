from ananse.utils import view_h5
import sys


def view(args):
    df = view_h5(args.infile, tfs=args.factors, fmt=args.format)
    index = True
    if args.format == "long":
        index = False

    if args.outfile is None:
        args.outfile = sys.stdout

    df.to_csv(args.outfile, sep="\t", index=index)
