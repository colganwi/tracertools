from importlib.metadata import version

from . import cli, seq

__all__ = ["seq", "cli"]

__version__ = version("tracertools")
