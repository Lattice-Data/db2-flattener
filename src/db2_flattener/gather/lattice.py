import io
import os
import sys
from collections import Counter
from urllib.parse import urljoin

import pandas as pd
import requests


class Connection:
    def __init__(self, mode: str):
        if not mode.upper().startswith("DB2_"):
            sys.exit(
                "ERROR: make sure your local variables start with DB2_ '(DB2_DEMO_SERVER, DB2_DEMO_KEY, etc...)'"
            )

        if not (
            os.environ.get(mode.upper() + "_KEY")
            and os.environ.get(mode.upper() + "_SECRET")
            and os.environ.get(mode.upper() + "_SERVER")
        ):
            sys.exit(
                "ERROR: "
                + mode.upper()
                + "_KEY "
                + mode.upper()
                + "_SECRET "
                + mode.upper()
                + "_SERVER "
                + "not all defined. "
                + "Try 'conda env config vars list' to list existing variables"
            )

        self.authid = os.environ.get(mode.upper() + "_KEY")
        self.authpw = os.environ.get(mode.upper() + "_SECRET")
        self.server = os.environ.get(mode.upper() + "_SERVER")
        if not self.server.endswith("/"):
            self.server += "/"
        self.headers = {
            "content-type": "application/json",
            "accept": "application/json",
        }
        self.auth = (self.authid, self.authpw)


MAX_REPORT_URL_LENGTH = 8000


def _report_urls(obj_type, filter_url, field_lst, connection, max_length=MAX_REPORT_URL_LENGTH):
    """
    Return report URLs that together cover filter_url, each at most max_length long.

    A URL that fits is returned as is. Otherwise one key's values are spread
    across as many URLs as needed, and every other filter term is repeated in
    each URL.

    Repeated values of one key are OR'd and different keys are AND'd, so
    splitting a single key while keeping the rest preserves the filter exactly:
    (a1 OR a2) AND b == (a1 AND b) OR (a2 AND b). The key split is the one with
    the most terms. A non-empty filter never produces an unfiltered URL; if the
    repeated terms plus any one split value can't fit, ValueError is raised
    instead of sending anything.
    """
    field_url = "".join(["&field=" + i for i in field_lst])

    def build(filter_part):
        return urljoin(
            connection.server,
            f"report/?type={obj_type}{filter_part}{field_url}&format=json&limit=all",
        )

    url = build(filter_url)
    if len(url) <= max_length:
        return [url]

    terms = [term for term in filter_url.split("&") if term]
    if not terms:
        raise ValueError(f"report URL for {obj_type} exceeds {max_length} characters")

    def key_of(term):
        # "status!=deleted" has key "status!", distinct from "status".
        return term.split("=", 1)[0]

    counts = Counter(key_of(term) for term in terms)
    split_key = max(counts, key=counts.get)  # ties: the key seen first
    fixed = "".join("&" + term for term in terms if key_of(term) != split_key)
    # A value repeated across chunks would return the same objects twice.
    split_pieces = list(dict.fromkeys("&" + term for term in terms if key_of(term) == split_key))

    budget = max_length - len(build(fixed))
    chunks, current = [], ""
    for piece in split_pieces:
        if len(piece) > budget:
            raise ValueError(
                f"report filter for {obj_type} can't fit in a {max_length}-character URL: "
                f"{piece[1:81]!r} plus the {len(fixed)} characters of other filter terms "
                "repeated in every request"
            )
        if current and len(current) + len(piece) > budget:
            chunks.append(current)
            current = ""
        current += piece
    chunks.append(current)
    return [build(fixed + chunk) for chunk in chunks]


def get_report(obj_type, filter_url, field_lst, connection):
    """
    Constructs a report url of fields in field_list for the objects determined by the filter and returns list of dictionaries from @graph
    URLs over MAX_REPORT_URL_LENGTH are split into as many filtered requests as needed (see _report_urls); results are returned as a single list of dictionaries from @graph
    """
    urls = _report_urls(obj_type, filter_url, field_lst, connection)
    graph = []
    for url in urls:
        obj = requests.get(url, auth=connection.auth)
        try:
            obj.raise_for_status()
        except requests.exceptions.HTTPError:
            if obj.status_code == 404 and "@graph" in obj.json():
                graph.extend(obj.json()["@graph"])
                continue
            raise
        graph.extend(obj.json().get("@graph"))
    if len(urls) > 1:
        graph = _dedupe_by_id(graph)
    return graph


def _dedupe_by_id(graph):
    """
    Drop repeat objects (same @id) from split requests, keeping the first.

    When the split key is a list field (e.g. samples), one object can match
    values in two chunks. Objects without an @id are kept as they are.
    """
    seen = set()
    unique = []
    for obj in graph:
        object_id = obj.get("@id") if isinstance(obj, dict) else None
        if object_id is not None:
            if object_id in seen:
                continue
            seen.add(object_id)
        unique.append(obj)
    return unique


def _separator_for_file(file_info):
    """Return the pandas read_csv separator, or None to sniff."""
    file_format = (file_info.get("file_format") or "").lower()
    if file_format == "tsv":
        return "\t"
    if file_format == "csv":
        return ","
    s3_uri = file_info.get("s3_uri") or ""
    if s3_uri.endswith(".tsv"):
        return "\t"
    if s3_uri.endswith(".csv"):
        return ","
    return None


def _read_csv_kwargs(file_info):
    sep = _separator_for_file(file_info)
    if sep is None:
        return {"sep": None, "engine": "python"}
    return {"sep": sep}


def download_file(object_id, connection):
    """
    Download a Lattice file object's bytes via @@download.

    object_id is a path like /tabular_files/<uuid>/.
    """
    url = urljoin(connection.server, object_id) + "@@download"
    response = requests.get(url, auth=connection.auth)
    response.raise_for_status()
    return response.content


def read_tabular_file(file_info, connection):
    """
    Read a TabularFile-like dict into a DataFrame without writing to disk.

    Tries fsspec/s3fs on s3_uri first, then Lattice @@download into memory.
    """
    read_kwargs = _read_csv_kwargs(file_info)
    s3_uri = file_info.get("s3_uri")
    object_id = file_info.get("@id")

    if s3_uri:
        try:
            return pd.read_csv(s3_uri, **read_kwargs)
        except Exception as exc:
            print(
                f"Warning: fsspec read of {s3_uri} failed ({exc}); "
                "falling back to Lattice @@download"
            )

    if not object_id:
        raise ValueError("guide RNA file has neither a readable s3_uri nor an @id")

    content = download_file(object_id, connection)
    return pd.read_csv(io.BytesIO(content), **read_kwargs)
