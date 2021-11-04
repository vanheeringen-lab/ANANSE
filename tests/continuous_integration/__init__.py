# run tests locally with:
# pytest -vv --disable-pytest-warnings
# pytest -vv --disable-pytest-warnings tests/continuous_integration/test_01*
# pytest -vv --disable-pytest-warnings -k [substring]

# test_dir = os.path.dirname(os.path.dirname(__file__))


def write_file(filename, lines):
    with open(filename, "w") as f:
        for line in lines:
            if not line.endswith("\n"):
                line = line + "\n"
            f.write(line)


# def write_bam(filename, lines):
#     tmp_sam = os.path.join(outdir, "tmp.sam")
#     write_file(tmp_sam, lines)
#     pysam.view(tmp_sam, "-b", "-o", filename, catch_stdout=False)
#     genomepy.utils.rm_rf(tmp_sam)


def compare_contents(file1, file2, ftype="bed"):
    if ftype == "bed":
        with open(file1) as f:
            contents1 = f.readlines()
        with open(file2) as f:
            contents2 = f.readlines()
    else:
        # contents1 = pysam.view(file1)
        # contents2 = pysam.view(file2)
        raise NotImplementedError("only beds can be compared")
    return contents1 == contents2
