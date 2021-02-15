import subprocess as sp

# run tests locally with: pytest -vv --disable-pytest-warnings

# TODO: apply to all code --> targets = ["ananse/", "tests/"]
targets = [
    "ananse/enhancer_binding.py",
    "ananse/commands/enhancer_binding.py",
    "tests/",
]


def test_import_ananse():
    import ananse

    assert str(ananse.__file__).endswith("ANANSE/ananse/__init__.py")


def test_black_formatting():
    sp.check_call(" ".join(["black setup.py"] + targets), shell=True)


def test_flake8_formatting():
    ret = sp.check_call(" ".join(["flake8 setup.py"] + targets), shell=True)
    assert ret == 0


# class TestBasics:
#     def test_import_ananse(self):
#         import ananse
