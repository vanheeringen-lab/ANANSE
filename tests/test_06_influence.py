import os
import pytest

import ananse.influence


@pytest.fixture
def influence_obj(outdir):
    i = ananse.influence.Influence(
        outfile=os.path.join(outdir, "influence.tsv"),
        degenes="",  # TODO
        GRN_target_file="",  # TODO
        edges=10,
    )
    return i


def test_read_top_interactions():
    pass  # TODO


def test_read_network():
    pass  # TODO


def test_difference():
    pass  # TODO


def test_targetScore():
    pass  # TODO


def test_filter_TF():
    pass  # TODO


def test_influence_init():
    pass  # TODO


def test_read_expression():
    # res = read_expression("tests/data/dge.tsv")
    # assert set(res.keys()) == {"ANPEP", "CD24", "COL6A3", "DAB2", "DMKN"}
    # assert res["ANPEP"].score - 7.44242618323665 < 0.001
    # assert res["ANPEP"].realfc - 7.44242618323665 < 0.001
    # assert res["ANPEP"].absfc - 7.44242618323665 < 0.001
    # assert res["COL6A3"].score == 0
    # assert res["COL6A3"].realfc - 11.0553152937569 < 0.001
    # assert res["COL6A3"].absfc - 11.0553152937569 < 0.001
    pass  # TODO


def test_save_reg_network():
    pass  # TODO


def test_run_target_score():
    pass  # TODO


def test_run_influence_score():
    pass  # TODO


def test_run_influence():
    pass  # TODO


def test_command_influence():
    pass  # TODO
