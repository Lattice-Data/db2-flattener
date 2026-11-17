"""get_report's long-URL split: every request stays filtered and within the limit."""

import uuid as uuid_mod
from unittest import mock
from urllib.parse import parse_qsl, urlparse

import pytest

from db2_flattener.gather import lattice
from db2_flattener.gather.lattice import MAX_REPORT_URL_LENGTH, get_report


class Conn:
    server = "https://lattice.test/"
    auth = None


FIELDS = ["@id", "uuid", "status", "aliases"]


def _capture(monkeypatch):
    sent = []

    def fake_get(url, auth=None):
        sent.append(url)
        query = parse_qsl(urlparse(url).query, keep_blank_values=True)
        graph = [{key: value} for key, value in query if key in ("uuid", "@id")]
        response = mock.Mock()
        response.json.return_value = {"@graph": graph}
        return response

    monkeypatch.setattr(lattice.requests, "get", fake_get)
    return sent


def _filter_terms(url):
    query = parse_qsl(urlparse(url).query, keep_blank_values=True)
    return [(key, value) for key, value in query if key in ("uuid", "@id")]


@pytest.mark.parametrize(
    ("key", "values"),
    [
        ("uuid", [str(uuid_mod.UUID(int=i)) for i in range(400)]),
        ("@id", [f"/raw_matrix_files/{uuid_mod.UUID(int=i)}/" for i in range(400)]),
    ],
)
def test_long_filter_is_split_into_bounded_filtered_requests(monkeypatch, key, values):
    sent = _capture(monkeypatch)
    filter_url = "".join(f"&{key}={value}" for value in values)

    graph = get_report("RawMatrixFile", filter_url, FIELDS, Conn())

    assert len(sent) > 1
    requested = []
    for url in sent:
        assert len(url) <= MAX_REPORT_URL_LENGTH
        terms = _filter_terms(url)
        assert terms, f"unfiltered request sent: {url}"
        assert all(value for _, value in terms), f"empty filter term in {url}"
        requested.extend(terms)
    assert requested == [(key, value) for value in values]
    assert [next(iter(obj.values())) for obj in graph] == values


def test_short_filter_is_one_unchanged_request(monkeypatch):
    sent = _capture(monkeypatch)
    get_report("MatrixFileSet", "&uuid=abc", ["@id"], Conn())
    assert sent == [
        "https://lattice.test/report/?type=MatrixFileSet&uuid=abc&field=@id&format=json&limit=all"
    ]


def test_mixed_filter_keys_are_all_kept(monkeypatch):
    sent = _capture(monkeypatch)
    terms = [("uuid" if i % 2 else "@id", f"value-{i:04d}-" + "x" * 30) for i in range(400)]
    get_report("Sample", "".join(f"&{k}={v}" for k, v in terms), FIELDS, Conn())
    assert len(sent) > 1
    assert [t for url in sent for t in _filter_terms(url)] == terms


def test_single_oversized_term_raises_instead_of_sending(monkeypatch):
    sent = _capture(monkeypatch)
    with pytest.raises(ValueError, match="too long"):
        get_report("Sample", "&uuid=" + "x" * 9000, FIELDS, Conn())
    assert sent == []
