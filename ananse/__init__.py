import os
import sys

from loguru import logger
from tqdm.auto import tqdm

from ._version import get_versions

# Remove default logger
logger.remove()
# Add logger
logger.add(sys.stderr, format="{time} | {level} | {message}", level="INFO")
# Combine tqdm output with loguru
logger.add(lambda msg: tqdm.write(msg, end=""))

# This is here to prevent very high memory usage on numpy import.
# On a machine with many cores, just importing numpy can result in up to
# 8GB of (virtual) memory. This wreaks havoc on management of the dask
# workers.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

__version__ = get_versions()["version"]
del get_versions

# we are using the em-dash (—) as a separator, since so far no one seems
# to be using that in their gene names (yet!). The em-dash is different
# from the en-dash (–) and hyphen (-).
SEPARATOR = "—"

PACKAGE_DIR = os.path.dirname(__file__)

from . import influence, network, peakpredictor, plot, utils, view  # noqa
