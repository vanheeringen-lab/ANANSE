import atexit
from os import getpid
import shutil
from tempfile import mkdtemp


# TODO: switch to the version in ananse.utils
def mytmpdir():
    if not hasattr(mytmpdir, "dir") or not mytmpdir.dir:
        mytmpdir.dir = mkdtemp(
            prefix=f"ANANSE_{getpid()}."
        )  # can be cleaned by clean_tmp()
        atexit.register(shutil.rmtree, mytmpdir.dir, ignore_errors=True)
    return mytmpdir.dir


from ._version import get_versions  # noqa: actual versioneer method

__version__ = get_versions()["version"]
del get_versions
