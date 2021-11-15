from ananse.view import view_h5
import sys


def view(args):
    df = view_h5(
        args.infile,
        args.tfs,
        args.regions,
        args.format,
        args.n,
        args.list_regions,
        args.list_tfs,
    )

    index = not (args.format == "long" or args.list_regions or args.list_tfs)
    header = not (args.list_regions or args.list_tfs)
    if args.outfile is None:
        args.outfile = sys.stdout

    df.to_csv(args.outfile, sep="\t", index=index, header=header)
