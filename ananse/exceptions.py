import logging
import sys

# use logger settings defined in __init__
logger = logging.getLogger(__name__)


def error(msg):
    logger.error(msg)
    sys.exit(1)
