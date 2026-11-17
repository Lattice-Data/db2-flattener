import io
import os
import sys
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

    A URL that fits is returned as is. Otherwise the filter is split on its
    "&"-separated terms, whatever their key (uuid=, @id=, ...), and packed into
    as many URLs as needed. A non-empty filter never produces an unfiltered
    URL; a single term too long to fit raises ValueError instead of being sent.
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

    budget = max_length - len(build(""))
    chunks, current = [], ""
    for term in terms:
        piece = "&" + term
        if len(piece) > budget:
            raise ValueError(
                f"report filter term {term[:80]!r} is too long to fit in a "
                f"{max_length}-character URL for {obj_type}"
            )
        if current and len(current) + len(piece) > budget:
            chunks.append(current)
            current = ""
        current += piece
    chunks.append(current)
    return [build(chunk) for chunk in chunks]


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
    return graph


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
