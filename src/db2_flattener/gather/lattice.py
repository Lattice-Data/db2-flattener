import os
import sys
from urllib.parse import urljoin

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
