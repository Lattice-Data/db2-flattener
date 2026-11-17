"""Download a MatrixFileSet's RawMatrixFiles from S3, renamed to geo_file_name when set."""

import os
import posixpath

import s3fs

from db2_flattener.gather.lattice import get_report
from db2_flattener.schema.constants import BASE_URL_OVERHEAD, MAX_URL_LENGTH
from db2_flattener.utils import extract_uuid_from_id

RAW_MATRIX_FILE_FIELDS = ["@id", "uuid", "s3_uri", "geo_file_name"]
MAX_FILTER_LENGTH = MAX_URL_LENGTH - BASE_URL_OVERHEAD


def target_file_name(raw_matrix_file: dict) -> str:
    """Local file name: geo_file_name when present, else the S3 object's own name."""
    geo_file_name = raw_matrix_file.get("geo_file_name")
    if geo_file_name:
        return geo_file_name
    return posixpath.basename(raw_matrix_file["s3_uri"].rstrip("/"))


def existing_file_warning(raw_matrix_file: dict) -> str:
    message = f"====== WARNING ====== For {raw_matrix_file['s3_uri']}, there already exists a file"
    geo_file_name = raw_matrix_file.get("geo_file_name")
    if geo_file_name:
        message += f" {geo_file_name}"
    return message


def fetch_raw_matrix_files(matrix_file_set_uuid: str, connection) -> list[dict]:
    """Return the RawMatrixFiles (@id, uuid, s3_uri, geo_file_name) of one MatrixFileSet."""
    matrix_file_sets = get_report(
        "MatrixFileSet",
        f"&uuid={matrix_file_set_uuid}",
        ["@id", "uuid", "raw_matrix_files"],
        connection,
    )
    if not matrix_file_sets:
        raise ValueError(f"MatrixFileSet {matrix_file_set_uuid} not found")

    uuids = []
    for ref in matrix_file_sets[0].get("raw_matrix_files") or []:
        ref_id = ref.get("@id", "") if isinstance(ref, dict) else ref
        if ref_id:
            uuids.append(extract_uuid_from_id(ref_id))
    if not uuids:
        return []

    wanted = list(dict.fromkeys(uuids))
    raw_matrix_files = []
    for filter_url in uuid_filter_chunks(wanted):
        raw_matrix_files.extend(
            get_report("RawMatrixFile", filter_url, RAW_MATRIX_FILE_FIELDS, connection)
        )

    # Belt and braces: never download a file this MatrixFileSet doesn't list,
    # even if a report request somehow came back unfiltered.
    wanted_set = set(wanted)
    raw_matrix_files = [rmf for rmf in raw_matrix_files if rmf.get("uuid") in wanted_set]
    return sorted(raw_matrix_files, key=lambda rmf: rmf.get("s3_uri") or "")


def uuid_filter_chunks(uuids: list[str], max_length: int = MAX_FILTER_LENGTH) -> list[str]:
    """
    Split uuids into "&uuid=..." filters each at most max_length characters.

    get_report's own long-URL split only understands "&@id=" filters; an
    oversized "&uuid=" filter makes it send a request with no filter at all,
    so every filter must already fit before it gets there.
    """
    chunks, current = [], ""
    for uuid in uuids:
        term = f"&uuid={uuid}"
        if current and len(current) + len(term) > max_length:
            chunks.append(current)
            current = ""
        current += term
    if current:
        chunks.append(current)
    return chunks


def download_matrix_file_set(matrix_file_set_uuid: str, out_directory: str, connection, *, fs=None):
    """
    Download every RawMatrixFile of a MatrixFileSet into out_directory.

    A file whose target name already exists in out_directory is not downloaded;
    a warning is printed instead. Returns (downloaded, skipped, failed) lists of
    s3_uri strings.
    """
    raw_matrix_files = fetch_raw_matrix_files(matrix_file_set_uuid, connection)
    print(f"Found {len(raw_matrix_files)} raw matrix files")

    os.makedirs(out_directory, exist_ok=True)
    filesystem = fs if fs is not None else s3fs.S3FileSystem()

    downloaded, skipped, failed = [], [], []
    for rmf in raw_matrix_files:
        s3_uri = rmf.get("s3_uri")
        if not s3_uri:
            print(f"====== WARNING ====== {rmf.get('@id')} has no s3_uri; skipping")
            failed.append(rmf.get("@id"))
            continue

        name = target_file_name(rmf)
        if not name or name in (".", "..") or os.path.basename(name) != name:
            print(f"====== WARNING ====== For {s3_uri}, invalid file name {name!r}; skipping")
            failed.append(s3_uri)
            continue

        destination = os.path.join(out_directory, name)
        if os.path.lexists(destination):
            print(existing_file_warning(rmf))
            skipped.append(s3_uri)
            continue

        # Download to a temporary name first so an interrupted transfer never
        # leaves a partial file that a re-run would mistake for a finished one.
        partial = destination + ".part"
        print(f"Downloading {s3_uri} -> {destination}")
        try:
            filesystem.get(s3_uri, partial)
            os.replace(partial, destination)
        except Exception as exc:
            print(f"====== ERROR ====== Failed to download {s3_uri}: {exc}")
            if os.path.exists(partial):
                os.remove(partial)
            failed.append(s3_uri)
            continue
        downloaded.append(s3_uri)

    return downloaded, skipped, failed
