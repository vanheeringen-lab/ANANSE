from ._version import get_versions
import os

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
