"""Resolve BIOHUB samples/guide template version suffixes (v#)."""

from __future__ import annotations

import re
from pathlib import PurePosixPath

import s3fs

# samples_template_v1.csv, guide_template_set2_v3.csv, etc.
BIOHUB_TEMPLATE_NAME_RE = re.compile(
    r"^(?:samples|guide)_template(?:_set\d+)?_v(\d+)\.csv$",
    re.IGNORECASE,
)


def normalize_biohub_version(value: str) -> str:
    """Accept ``3`` or ``v3`` and return ``v3``."""
    text = str(value).strip()
    if not text:
        raise ValueError("biohub version must be a positive integer like 1 or v1")
    if text.lower().startswith("v"):
        text = text[1:].strip()
    if not text.isdigit() or int(text) < 1:
        raise ValueError(f"invalid biohub version {value!r}; expected e.g. 1 or v1")
    return f"v{int(text)}"


def versions_from_metadata_names(names) -> set[int]:
    """Collect v# integers from BIOHUB template basenames."""
    versions = set()
    for name in names:
        basename = PurePosixPath(str(name)).name
        match = BIOHUB_TEMPLATE_NAME_RE.fullmatch(basename)
        if match:
            versions.add(int(match.group(1)))
    return versions


def next_biohub_version_from_versions(versions: set[int]) -> str:
    """
    From versions found under metadata/, return the next template suffix.

    Empty → v1. Otherwise highest found v# + 1.
    """
    if not versions:
        return "v1"
    return f"v{max(versions) + 1}"


def resolve_biohub_version_from_s3(s3_order_dir: str, *, fs=None) -> str:
    """
    List ``{s3_order_dir}/metadata/`` and return the next BIOHUB v#.

    Missing/empty metadata or no matching template names → ``v1``.
    Order-dir access failures raise rather than defaulting.
    """
    order = str(s3_order_dir).rstrip("/")
    if not order.startswith("s3://"):
        raise ValueError(f"s3 order directory must be an s3:// URI, got {s3_order_dir!r}")

    filesystem = fs if fs is not None else s3fs.S3FileSystem()
    try:
        if not filesystem.exists(order):
            raise FileNotFoundError(f"S3 order directory not found: {order}")
    except FileNotFoundError:
        raise
    except Exception as exc:
        raise RuntimeError(f"Failed to access S3 order directory {order}: {exc}") from exc

    metadata = f"{order}/metadata"
    try:
        if not filesystem.exists(metadata):
            return "v1"
    except Exception as exc:
        raise RuntimeError(f"Failed to access {metadata}: {exc}") from exc

    try:
        entries = filesystem.ls(metadata, detail=False)
    except FileNotFoundError:
        return "v1"
    except Exception as exc:
        raise RuntimeError(f"Failed to list {metadata}: {exc}") from exc

    return next_biohub_version_from_versions(versions_from_metadata_names(entries))


def resolve_biohub_version(
    *, s3_order_dir: str | None = None, biohub_version: str | None = None, fs=None
) -> str:
    """Resolve from exactly one of ``biohub_version`` or ``s3_order_dir``."""
    if biohub_version is not None and s3_order_dir is not None:
        raise ValueError("provide only one of biohub_version or s3_order_dir")
    if biohub_version is not None:
        return normalize_biohub_version(biohub_version)
    if s3_order_dir is not None:
        return resolve_biohub_version_from_s3(s3_order_dir, fs=fs)
    raise ValueError("one of biohub_version or s3_order_dir is required")
