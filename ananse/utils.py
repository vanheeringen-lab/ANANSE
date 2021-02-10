import os
import warnings
import shutil
import tempfile

from pybedtools import BedTool


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
def sort_bed(bed):
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
def merge_bed_files(list_of_beds, merged_bed):
    """
    merge any number of bed files
    """
    bed = BedTool(list_of_beds[0])
    if list_of_beds[1:]:
        bed = bed.cat(*list_of_beds[1:])
    # bed.sort(output=merged_bed)
    bed.moveto(merged_bed)


@shhh_bedtool
def count_reads(bam_list, peakfile, bed_output):
    """
    Count bam reads in putative enhancer regions
    """
    bed = BedTool(peakfile)
    bed.multi_bam_coverage(bams=bam_list, output=bed_output)


def cleanpath(path):
    """Expand any path input to a literal path output"""
    return \
        os.path.abspath(  # expand relative paths (inc './' and '../')
            os.path.expanduser(  # expand '~'
                os.path.expandvars(  # expand '$VARIABLES'
                    path)))
