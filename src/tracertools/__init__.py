from importlib.metadata import version

from . import cli, seq, solver, tree, config, root

__all__ = ["seq", "cli", "solver", "tree", "config", "root"]

__version__ = version("tracertools")
