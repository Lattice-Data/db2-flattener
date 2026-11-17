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


def _matches(obj, query_terms):
    """snovault-style filter: same key OR'd, different keys AND'd, "key!" negated."""
    by_key = {}
    for key, value in query_terms:
        by_key.setdefault(key, set()).add(value)
    for key, values in by_key.items():
        if key.endswith("!"):
            if obj.get(key[:-1]) in values:
                return False
        elif obj.get(key) not in values:
            return False
    return True


FILTER_SKIP = {"type", "field", "format", "limit"}


def _query_filter(url):
    query = parse_qsl(urlparse(url).query, keep_blank_values=True)
    return [(k, v) for k, v in query if k not in FILTER_SKIP]


def _serve(monkeypatch, database):
    sent = []

    def fake_get(url, auth=None):
        sent.append(url)
        terms = _query_filter(url)
        response = mock.Mock()
        response.json.return_value = {"@graph": [o for o in database if _matches(o, terms)]}
        return response

    monkeypatch.setattr(lattice.requests, "get", fake_get)
    return sent


DATABASE = [
    {
        "@id": f"/samples/{uuid_mod.UUID(int=i)}/",
        "uuid": str(uuid_mod.UUID(int=i)),
        "status": ["released", "in progress", "deleted"][i % 3],
        "lab": ["/labs/a/", "/labs/b/"][i % 2],
    }
    for i in range(600)
]


def _expected(filter_url):
    terms = parse_qsl(filter_url.lstrip("&"), keep_blank_values=True)
    return sorted(o["uuid"] for o in DATABASE if _matches(o, terms))


@pytest.mark.parametrize(
    "filter_url",
    [
        pytest.param(
            "&status=released" + "".join(f"&@id={o['@id']}" for o in DATABASE[:400]),
            id="one-and-term-plus-many-ids",
        ),
        pytest.param(
            "".join(f"&@id={o['@id']}" for o in DATABASE[:400]) + "&status!=deleted",
            id="negated-term-after-ids",
        ),
        pytest.param(
            "&status=released&status=in progress&lab=/labs/a/"
            + "".join(f"&uuid={o['uuid']}" for o in DATABASE[:400]),
            id="two-multi-value-keys",
        ),
        pytest.param(
            "".join(f"&uuid={o['uuid']}" for o in DATABASE[:300])
            + "".join(f"&@id={o['@id']}" for o in DATABASE[150:170]),
            id="uuid-and-id-both-required",
        ),
    ],
)
def test_split_preserves_filter_meaning(monkeypatch, filter_url):
    """The union of the split requests returns exactly what the unsplit filter means."""
    sent = _serve(monkeypatch, DATABASE)

    graph = get_report("Sample", filter_url, FIELDS, Conn())

    assert len(sent) > 1
    assert all(len(url) <= MAX_REPORT_URL_LENGTH for url in sent)
    assert sorted(o["uuid"] for o in graph) == _expected(filter_url)


def test_split_repeats_every_other_term_in_each_request(monkeypatch):
    sent = _serve(monkeypatch, DATABASE)
    ids = [o["@id"] for o in DATABASE[:400]]
    filter_url = "&status=released&lab=/labs/a/&status!=deleted" + "".join(f"&@id={i}" for i in ids)

    get_report("Sample", filter_url, FIELDS, Conn())

    assert len(sent) > 1
    split_values = []
    for url in sent:
        terms = _query_filter(url)
        others = [t for t in terms if t[0] != "@id"]
        assert others == [("status", "released"), ("lab", "/labs/a/"), ("status!", "deleted")]
        split_values.extend(v for k, v in terms if k == "@id")
    assert split_values == ids


def test_single_oversized_term_raises_instead_of_sending(monkeypatch):
    sent = _capture(monkeypatch)
    with pytest.raises(ValueError, match="can't fit"):
        get_report("Sample", "&uuid=" + "x" * 9000, FIELDS, Conn())
    assert sent == []


def test_two_large_keys_raise_instead_of_changing_meaning(monkeypatch):
    """Only one key is split; when the repeated keys alone don't fit, nothing is sent."""
    sent = _capture(monkeypatch)
    filter_url = "".join(f"&uuid={o['uuid']}" for o in DATABASE[:300]) + "".join(
        f"&@id={o['@id']}" for o in DATABASE[:299]
    )
    with pytest.raises(ValueError, match="can't fit"):
        get_report("Sample", filter_url, FIELDS, Conn())
    assert sent == []
