import os
import subprocess as sp

targets = ["setup.py", "ananse/", "tests/"]


def test_import_ananse():
    import ananse

    assert isinstance(ananse.__version__, str)
    assert isinstance(ananse.SEPARATOR, str)
    assert os.path.exists(ananse.PACKAGE_DIR)


def test_black_linting():
    sp.check_call(" ".join(["black --check"] + targets), shell=True)


def test_flake8_linting():
    ret = sp.check_call(" ".join(["flake8 --ignore=W503,B950"] + targets), shell=True)
    assert ret == 0
