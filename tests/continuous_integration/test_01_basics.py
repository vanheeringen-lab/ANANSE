# TODO: apply to all code --> targets = ["ananse/", "tests/"]
targets = [
    "ananse/enhancer_binding.py",
    "ananse/commands/enhancer_binding.py",
    "tests/",
]

class TestBasics:
    def test_import_ananse(self):
        import ananse
