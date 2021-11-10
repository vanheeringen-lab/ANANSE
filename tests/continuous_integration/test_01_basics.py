import subprocess as sp

targets = ["setup.py", "ananse/", "tests/"]


def test_import_ananse():
    import ananse

    assert str(ananse.__file__).endswith("ANANSE/ananse/__init__.py")


def test_black_linting():
    sp.check_call(" ".join(["black --check"] + targets), shell=True)


def test_flake8_linting():
    ret = sp.check_call(" ".join(["flake8"] + targets), shell=True)
    assert ret == 0
