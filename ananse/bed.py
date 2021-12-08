import os
import tempfile
import shutil
import warnings

import genomepy
from loguru import logger
import pandas as pd
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
def bed_sort(bed):
    """
    Sort a bed file.
    """
    tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
    try:
        tmpfile = os.path.join(tmpdir, os.path.basename(bed))
        BedTool(bed).sort(output=tmpfile)
        shutil.copy2(tmpfile, bed)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@shhh_bedtool
def bed_merge(list_of_beds, merged_bed):
    """
    merge any number of bed files (merges overlapping regions)
    """
    bed = BedTool(list_of_beds[0])
    if list_of_beds[1:]:
        bed = bed.cat(*list_of_beds[1:])
    bed.saveas(merged_bed)


class CombineBedFiles:
    def __init__(self, genome, peakfiles, verbose=True):
        self.genome = genome
        self.list_of_peakfiles = (
            peakfiles if isinstance(peakfiles, list) else [peakfiles]
        )
        self.verbose = verbose

    @staticmethod
    def is_narrowpeak(bed, check_values=True):
        """
        Check BED type by column count.
        Check if peak values are not all zeroes unless check_values is False.

        Accepts a BED file (including narrowPeak, broadPeak, etc.)
        Returns bool
        """
        with open(bed) as b:
            for line in b:
                if line.startswith("#"):
                    continue
                line = line.split("\t")
                cols = len(line)
                break

        # narrowPeak has 10 columns
        # and the peak column is >= 0
        if cols != 10 or int(line[9]) < 0:
            return False

        if not check_values:
            return True

        # check if the peak values aren't all zeroes
        summit_values = 0
        sample_size = 20  # check an arbitrary number of lines
        with open(bed) as b:
            for n, line in enumerate(b):
                if line.startswith("#"):
                    continue
                line = line.split("\t")
                peak_val = int(line[9])

                # value must be >=0
                if peak_val < 0:
                    return False
                summit_values += peak_val
                if n >= sample_size:
                    break

        if summit_values > 0:
            return True
        return False

    @staticmethod
    def bed_resize(
        genome,
        bed_in,
        bed_out,
        width=200,
        narrowpeak=False,
        fix_outliers=False,
        output_bed3=True,
        verbose=True,
    ):
        """
        Set bed region width.

        If the input bed is a narrowPeak file (narrowpeak=True),
        center region on the summit (start+peak).
        Otherwise center on the middle of the region.

        If fix_outliers is set to True, shift regions to fit their chromosomes.
        Otherwise drop these regions.

        If output_bed3 is set to False, output the whole bed file.
        """
        half_seqlen = width // 2
        chrom_sizes = genomepy.Genome(genome).sizes
        missing_chrm = []

        if narrowpeak:

            def get_summit(_start, _, summit_offset):
                return _start + int(summit_offset)

            summit_col = 9
        else:

            def get_summit(_start, _end, _):
                return (_start + _end) // 2

            summit_col = 0  # unused

        with open(bed_in) as old, open(bed_out, "w") as new:
            for line in old:
                if line.startswith("#"):
                    continue
                line = line.split("\t")
                chrm = str(line[0])
                if chrm not in chrom_sizes.keys():
                    missing_chrm.append(chrm)
                    continue
                start = int(line[1])
                end = int(line[2])
                rest = line[3:] if not output_bed3 else []

                chrm_len = chrom_sizes[chrm]
                if width == end - start:
                    nstart = str(start)
                    nend = str(end)
                elif chrm_len <= width:
                    if not fix_outliers:
                        continue
                    nstart = str(0)
                    nend = str(chrm_len)
                else:
                    summit = get_summit(start, end, line[summit_col])
                    if not fix_outliers:
                        nstart = str(summit - half_seqlen)
                        nend = str(summit + half_seqlen)
                        if int(nstart) < 0 or int(nend) > chrm_len:
                            continue
                    else:
                        # adjust the summit for the chromosome boundaries
                        summit = max(summit, 0 + half_seqlen)
                        summit = min(summit, chrm_len - half_seqlen)

                        nstart = str(summit - half_seqlen)
                        nend = str(summit + half_seqlen)

                new.write("\t".join([chrm, nstart, nend] + rest) + "\n")

        if missing_chrm and verbose:
            logger.warning(
                "The following contigs were present in "
                + f"'{os.path.basename(bed_in)}', "
                + "but were missing in the genome file: "
                + f"{', '.join(list(set(missing_chrm)))}\n"
            )
        return bed_out

    def run(self, outfile, width=200, force=False):
        if force or not os.path.exists(outfile):
            if self.verbose:
                logger.info("Combining bed files")
            tmpdir = tempfile.mkdtemp(prefix="ANANSE_")
            try:
                list_of_beds = []
                for peakfile in self.list_of_peakfiles:
                    # use narrowPeak Peak location for region centering if possible
                    is_np = self.is_narrowpeak(peakfile)
                    resized_peakfile = os.path.join(tmpdir, os.path.basename(peakfile))

                    # resize each BED region to 200 BP
                    self.bed_resize(
                        genome=self.genome,
                        bed_in=peakfile,
                        bed_out=resized_peakfile,
                        width=width,
                        narrowpeak=is_np,
                        verbose=self.verbose,
                    )
                    bed_sort(resized_peakfile)
                    list_of_beds.append(resized_peakfile)

                # merge resized beds into one
                merged_bed = os.path.join(tmpdir, "merged")
                bed_merge(list_of_beds=list_of_beds, merged_bed=merged_bed)

                shutil.copy2(merged_bed, outfile)
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)


def df_to_bedtool(df_in):
    """sort, save and load a dataframe into a BedTool class"""
    tmpfile = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    df_in[1] = df_in[1].astype(int)
    df_in.sort_values([0, 1], inplace=True)
    df_in.to_csv(tmpfile, sep="\t", index=False, header=False)
    return BedTool(tmpfile)


def map_counts(regions: list, counts: pd.DataFrame) -> pd.DataFrame:
    """
    map regions in the countsfile to intersecting regions in the regions list.
    counts mapping to multiple regions are duplicated.
    missing regions are set to 0.
    """
    r = pd.Series(regions).str.split(r":|-", expand=True)
    r = df_to_bedtool(r)

    c = counts.index.str.split(r":|-", expand=True).to_frame()
    c = df_to_bedtool(c)

    # get overlapping regions, with the original location in the counts file
    intersect = r.intersect(c, sorted=True, wa=True, wb=True, loj=True).fn

    # map intersecting count locations to their corresponding region (can duplicate rows)
    df = pd.read_table(intersect, sep="\t", comment="#", header=None, dtype=str)
    df.columns = ["chrom", "start_r", "end_r", "_", "start_c", "end_c"]
    df["region"] = df["chrom"] + ":" + df["start_r"] + "-" + df["end_r"]
    df.index = df["chrom"] + ":" + df["start_c"] + "-" + df["end_c"]
    df = df[["region"]]

    # replace the locations with regions + add missing regions
    mapped_counts = df.join(counts).set_index("region", drop=True)
    # fix formatting
    mapped_counts.index.name = None
    mapped_counts.fillna(0.0, inplace=True)
    return mapped_counts
