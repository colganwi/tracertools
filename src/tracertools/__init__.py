from importlib.metadata import version

from . import cli, config, root, seq, solver, tree

__all__ = ["seq", "cli", "solver", "tree", "config", "root"]

__version__ = version("tracertools")
