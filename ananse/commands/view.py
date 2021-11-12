from ananse.utils import view_h5
import sys


def view(args):
    df = view_h5(args.infile, args.tfs, args.format, args.regions, args.factors)

    index = not (args.format == "long" or args.regions or args.factors)
    header = not (args.regions or args.factors)
    if args.outfile is None:
        args.outfile = sys.stdout

    df.to_csv(args.outfile, sep="\t", index=index, header=header)
