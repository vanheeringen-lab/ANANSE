import pytest

import ananse.peakpredictor


# _check_input_regions
# _check_input_files


def test_peakpredictor():
    with pytest.raises(ValueError):
        ananse.peakpredictor.PeakPredictor(reference="ananse/db/default_reference")

    p = ananse.peakpredictor.PeakPredictor(
        reference="ananse/db/default_reference",
        atac_bams=["tests/data/binding/bam1.bam"],
        regions=["chr1:1010-1020"],  # 1 (tiny) region
        genome="tests/data/binding/hg38_testgenome.fa",  # 1 chr
        pfmfile="tests/data/binding/test.pfm",  # 1 motif, 1 factor
        pfmscorefile="tests/data/binding/scan.tsv",  # precomputed
        ncore=1,
    )

    # boring stuff
    assert p.data_dir == "ananse/db/default_reference"
    assert p.genome == "tests/data/binding/hg38_testgenome.fa"
    assert p.pfmfile == "tests/data/binding/test.pfm"
    assert p.ncore == 1

    # set in test.motif2factors.txt
    assert p.factors() == ["SOX12"]
    assert p.f2m == {"SOX12": ["GM.5.0.Sox.0001"]}
    assert "GM.5.0.Sox.0001" in p.motifs

    # parsed from input
    assert p.regions == ["chr1:1010-1020"]
    assert p.region_type == "custom"
    assert p.species == "human"  # because hg38 is in the genome name

    assert p._histone_data is None
    assert "chr1:1010-1020" in p._atac_data.index
    assert p._X_columns == ["ATAC", "motif"]
    assert p._model_type == "ATAC_motif"

    assert len(p.motif_graph) == 0  # 1 TF -> no edges


# predict_peaks
