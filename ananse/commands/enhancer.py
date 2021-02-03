import ananse.utils
import ananse.enhancer


def enhancer(**kwargs):
    """
    CLI parser for ananse.enhancer.Enhancer
    """
    # print(dict(**kwargs))
    # print(kwargs.get("genome"))

    genome = kwargs.get("genome")
    genome = ananse.utils.check_genome(genome)

    # etype = kwargs.get("type")
    # etype = ananse.utils.check_type(etype)

    bams = kwargs.get("bam")
    bams = ananse.utils.check_bam(bam)

    peaks = kwargs.get("peak")
    peaks = ananse.utils.check_file(peak, "peak")

    output = kwargs.get("output")
    output = ananse.utils.check_output(output, "enhancer.bed")

    width = kwargs.get("width")

    summit_column = kwargs.get("summit")

    b = ananse.enhancer.Enhancer(
        genome=genome,
        bams=bams,
        peaks=peaks,
        output=output,
        width=width,
        summit_column=summit_column
    )
    b.run_enhancer()
