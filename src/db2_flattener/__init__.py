"""db2-flattener: Flatten Lattice DB2 MatrixFileSet data to CSV."""

__version__ = "0.0.1"

from db2_flattener.flatten.flattener import DB2Flattener

__all__ = ["DB2Flattener", "__version__"]
