import atexit
import logging
import os
import shutil
from tempfile import mkdtemp

import pysam
from genomepy import Genome
from genomepy.utils import mkdir_p
from ananse import exceptions

# use logger settings defined in __init__
logger = logging.getLogger(__name__)


def mytmpdir():
    if not hasattr(mytmpdir, "dir") or not mytmpdir.dir:
        mytmpdir.dir = mkdtemp(prefix="ananse.{0}.".format(os.getpid()))
        atexit.register(shutil.rmtree, mytmpdir.dir)
    return mytmpdir.dir


def set_width(genome, bed_in, bed_out, summit_col=-1, seqlen=200):
    """
    Set BED region width to SEQLEN bp.
    If summit_col is set, centers on the summit, else on the region center.
    Shifts the region to fit the chromosome near edges.
    """
    half_seqlen = seqlen // 2
    chrom_sizes = Genome(genome).sizes

    if summit_col != -1:
        def get_summit(s, e, summit_offset):
            return s + int(summit_offset)
        logger.info(f"Using column {summit_col} as summit")
    else:
        def get_summit(s, e, sc):
            return (s + e) // 2

    missing_chrm = set()
    with open(bed_in) as old, open(bed_out, "w") as new:
        for line in old:
            line = line.split()
            chrm = str(line[0])
            if chrm not in chrom_sizes.keys():
                if chrm not in missing_chrm:
                    missing_chrm.add(chrm)
                    logger.info(f"{chrm} found in BED file that is missing from the genome. Skipping...")
                continue

            start = int(line[1])
            end = int(line[2])
            rest = line[3:]

            chrm_len = chrom_sizes[chrm]
            if seqlen == end - start:
                nstart = str(start)
                nend = str(end)
            elif chrm_len <= seqlen:
                nstart = str(0)
                nend = str(chrm_len)
            else:
                summit = get_summit(start, end, line[summit_col])

                # adjust the summit for the chromosome boundaries
                summit = max(summit, 0 + half_seqlen)
                summit = min(summit, chrm_len - half_seqlen)

                nstart = str(summit - half_seqlen)
                nend = str(summit + half_seqlen)

            new.write("\t".join([chrm, nstart, nend] + rest) + "\n")

    return bed_out


def check_genome(genome):
    check_arg(genome, "genome")
    if ".fa" not in genome.lower():
        exceptions.error("genome file must be of type FASTA")
    genome = check_file(genome, "genome")
    return genome


def check_type(etype):
    """
    Check if any accepted_type is a substring of etype.
    Case insensitive.
    e.g. atac-seq is recognized as ATAC.
    """
    check_arg(etype, "type")
    accepted_types = ["H3K27ac", "p300", "ATAC"]
    for accepted_type in accepted_types:
        if accepted_type.lower() in etype.lower():
            return accepted_type
    else:
        exceptions.error(f"type not recognized: {etype}")


def check_bam(bam):
    """check bam and bam index existence. Generate the latter if needed."""
    bam = check_file(bam, "bam")
    if not os.path.exists(f"{bam}.bai"):
        logger.info(f"Generating BAM index for {os.path.basename(bam)}")
        pysam.index(bam)
    return bam


def check_output(output, filename):
    """
    accepts an output file or directory
    creates the output directory
    return an output filepath
    """
    output = cleanpath(output)

    # obtain output dir and set output filename if unset
    if "." in os.path.splitext(output)[1]:
        output_dir = os.path.dirname(output)
    else:
        output_dir = output
        output = os.path.join(output_dir, filename)

    # create the output directory
    if not os.path.exists(output_dir):
        mkdir_p(output)

    return output


def check_arg(arg, name):
    if not arg:
        exceptions.error(f"argument missing: {name}")


def check_file(file, name):
    """return expanded path if it exists, else an error"""
    check_arg(file, name)
    file = cleanpath(file)
    if not os.path.exists(file):
        exceptions.error(f"file not found: {file}")
    return file


def cleanpath(path):
    """expand path"""
    return os.path.abspath(os.path.expanduser(path))
