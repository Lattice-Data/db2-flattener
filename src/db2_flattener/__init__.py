"""db2-flattener: Flattener utilities for DB2."""

__version__ = "0.0.1"


def flatten(nested: dict, sep: str = ".") -> dict:
    """Flatten a nested dict. Placeholder implementation."""
    out: dict = {}

    def _walk(obj, prefix: str) -> None:
        if isinstance(obj, dict):
            for k, v in obj.items():
                key = f"{prefix}{sep}{k}" if prefix else str(k)
                _walk(v, key)
        else:
            out[prefix] = obj

    _walk(nested, "")
    return out
