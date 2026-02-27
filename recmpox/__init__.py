"""RecMpox: Classify consensus mpox genomes at diagnostic SNPs (recombinant calling)."""

from ._version import __version__

from .recmpox import main

__all__ = ["main"]
