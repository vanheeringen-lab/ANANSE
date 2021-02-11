import os
import warnings
import shutil
import tempfile

from pybedtools import BedTool
import pysam


def shhh_bedtool(func):
    """
    Decorator that silences pybedtools RuntimeWarnings such as
    `line buffering (buffering=1) isn't supported in binary mode`
    """
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            func(*args, **kwargs)
    return wrapper


@shhh_bedtool
def bed_sort(bed):
    """
    Sort a bed file.
    """
    tmpdir = tempfile.mkdtemp()
    tmpfile = os.path.join(tmpdir, "tmpfile")
    try:
        BedTool(bed).sort(output=tmpfile)
        shutil.copy2(tmpfile, bed)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@shhh_bedtool
def bed_merge(list_of_beds, merged_bed):
    """
    merge any number of bed files
    """
    bed = BedTool(list_of_beds[0])
    if list_of_beds[1:]:
        bed = bed.cat(*list_of_beds[1:])
    # bed.sort(output=merged_bed)
    bed.moveto(merged_bed)


@shhh_bedtool
def count_reads(bams, peakfile, bed_output):
    """
    Count bam reads in putative enhancer regions
    """
    bed = BedTool(peakfile)
    bam_list = bams if isinstance(bams, list) else [bams]
    bed.multi_bam_coverage(bams=bam_list, output=bed_output)


def bam_index(bam, force=True):
    if force or not os.path.exists(f"{bam}.bai"):
        pysam.index(bam)


def bam_sort(bam):
    tmpdir = tempfile.mkdtemp()
    try:
        sorted_bam = os.path.join(tmpdir, os.path.basename(bam))
        pysam.sort("-o", sorted_bam, bam)

        shutil.copy2(sorted_bam, bam)
        bam_index(bam)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def bam_merge(list_of_bams, merged_bam):
    """
    merge any number of bam files
    """
    if len(list_of_bams) > 1:
        pysam.merge(merged_bam, *list_of_bams)
        bam_index(merged_bam)
    else:
        # os.symlink() doesn't work with multi_bam_coverage()
        bam = list_of_bams[0]
        shutil.copy2(bam, merged_bam)
        shutil.copy2(f"{bam}.bai", f"{merged_bam}.bai")


def cleanpath(path):
    """Expand any path input to a literal path output"""
    return \
        os.path.abspath(  # expand relative paths (inc './' and '../')
            os.path.expanduser(  # expand '~'
                os.path.expandvars(  # expand '$VARIABLES'
                    path)))
