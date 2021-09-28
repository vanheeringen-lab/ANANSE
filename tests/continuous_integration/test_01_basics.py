import subprocess as sp

# run tests locally with:
# pytest -vv --disable-pytest-warnings
# pytest -vv --disable-pytest-warnings tests/continuous_integration/test_01*
# pytest -vv --disable-pytest-warnings -k [substring]

targets = ["ananse/", "tests/"]


def test_import_ananse():
    import ananse

    assert str(ananse.__file__).endswith("ANANSE/ananse/__init__.py")


def test_black_linting():
    sp.check_call(" ".join(["black --check setup.py"] + targets), shell=True)


def test_flake8_linting():
    ret = sp.check_call(" ".join(["flake8 setup.py"] + targets), shell=True)
    assert ret == 0
