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


def get_report(obj_type, filter_url, field_lst, connection):
    """
    Constructs a report url of fields in field_list for the objects determined by the filter and returns list of dictionaries from @graph
    Will split url into two requests if url is > 8000 characters, and still return a single list of dictionaries from @graph
    """
    field_url = "".join(["&field=" + i for i in field_lst])
    url1 = urljoin(
        connection.server, f"report/?type={obj_type}{filter_url}{field_url}&format=json&limit=all"
    )
    urls = []
    if len(url1) > 8000:
        filter_lst = filter_url.split("&@id=")
        filter_url1 = "".join(["&@id=" + i for i in filter_lst[: len(filter_lst) // 2]])
        filter_url2 = "".join(["&@id=" + i for i in filter_lst[len(filter_lst) // 2 :]])
        url1 = urljoin(
            connection.server,
            f"report/?type={obj_type}{filter_url1}{field_url}&format=json&limit=all",
        )
        url2 = urljoin(
            connection.server,
            f"report/?type={obj_type}{filter_url2}{field_url}&format=json&limit=all",
        )
        urls.append(url2)
    urls.append(url1)
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
    """Return the pandas read_csv separator for a TabularFile-like dict."""
    file_format = (file_info.get("file_format") or "").lower()
    if file_format == "tsv":
        return "\t"
    s3_uri = file_info.get("s3_uri") or ""
    if s3_uri.endswith(".tsv"):
        return "\t"
    return ","


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
    sep = _separator_for_file(file_info)
    s3_uri = file_info.get("s3_uri")
    object_id = file_info.get("@id")

    if s3_uri:
        try:
            return pd.read_csv(s3_uri, sep=sep)
        except Exception as exc:
            print(
                f"Warning: fsspec read of {s3_uri} failed ({exc}); "
                "falling back to Lattice @@download"
            )

    if not object_id:
        raise ValueError("guide RNA file has neither a readable s3_uri nor an @id")

    content = download_file(object_id, connection)
    return pd.read_csv(io.BytesIO(content), sep=sep)
