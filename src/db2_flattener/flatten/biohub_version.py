"""Resolve BIOHUB samples/guide template version suffixes (v#).

Versions are tracked independently per template type (``samples`` / ``guide``)
and per set number. Only basenames of the form
``{samples|guide}_template_set{N}_v{M}.csv`` count; anything else is ignored.
"""

from __future__ import annotations

import re
from pathlib import PurePosixPath

import s3fs

# samples_template_set1_v1.csv, guide_template_set2_v3.csv, etc.
BIOHUB_TEMPLATE_NAME_RE = re.compile(
    r"^(samples|guide)_template_set(\d+)_v(\d+)\.csv$",
    re.IGNORECASE,
)

BIOHUB_TEMPLATE_KINDS = ("samples", "guide")


def latest_versions_from_metadata_names(names) -> dict[tuple[str, int], int]:
    """
    Map ``(kind, set_number)`` to the highest v# found among template basenames.

    ``kind`` is lowercased (``samples`` or ``guide``). Names without a set
    number or that do not otherwise match are ignored.
    """
    latest: dict[tuple[str, int], int] = {}
    for name in names:
        basename = PurePosixPath(str(name)).name
        match = BIOHUB_TEMPLATE_NAME_RE.fullmatch(basename)
        if not match:
            continue
        key = (match.group(1).lower(), int(match.group(2)))
        version = int(match.group(3))
        if version > latest.get(key, 0):
            latest[key] = version
    return latest


def next_biohub_version(
    latest_versions: dict[tuple[str, int], int] | None, kind: str, set_number: int
) -> str:
    """
    Return the next ``v#`` for one template type and set.

    No prior version for ``(kind, set_number)`` → ``v1``; otherwise highest + 1.
    """
    if kind not in BIOHUB_TEMPLATE_KINDS:
        raise ValueError(f"unknown BIOHUB template kind {kind!r}")
    current = (latest_versions or {}).get((kind, int(set_number)), 0)
    return f"v{current + 1}"


def latest_biohub_versions_from_s3(s3_order_dir: str, *, fs=None) -> dict[tuple[str, int], int]:
    """
    List ``{s3_order_dir}/metadata/`` and return highest v# per ``(kind, set)``.

    Missing/empty metadata or no matching template names → ``{}`` (all v1).
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
            return {}
    except Exception as exc:
        raise RuntimeError(f"Failed to access {metadata}: {exc}") from exc

    try:
        entries = filesystem.ls(metadata, detail=False)
    except FileNotFoundError:
        return {}
    except Exception as exc:
        raise RuntimeError(f"Failed to list {metadata}: {exc}") from exc

    return latest_versions_from_metadata_names(entries)
