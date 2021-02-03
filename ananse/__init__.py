import logging

from ananse import commands
# from ._version import get_versions

__version__ = "0.2.0"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.propagate = 0

# nice format
screen_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# Log to screen
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(screen_formatter)
logger.addHandler(sh)

# # versioneer
# __version__ = get_versions()["version"]
# del get_versions
