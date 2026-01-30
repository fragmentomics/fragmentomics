"""
Logging configuration for fragmentomics.
"""

import logging
import sys


def setup_logging(
    level: int = logging.INFO,
    format_string: str | None = None,
) -> None:
    """
    Configure logging for fragmentomics.

    Parameters
    ----------
    level : int, default logging.INFO
        Logging level
    format_string : str, optional
        Custom format string
    """
    if format_string is None:
        format_string = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

    logging.basicConfig(
        level=level,
        format=format_string,
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger for a module.

    Parameters
    ----------
    name : str
        Logger name (usually __name__)

    Returns
    -------
    logging.Logger
    """
    return logging.getLogger(name)
