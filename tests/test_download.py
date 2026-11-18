import uuid as uuid_mod
from unittest import mock
from urllib.parse import parse_qs, urlparse

import pytest

from db2_flattener.download import cli, downloader
from db2_flattener.download.downloader import (
    download_matrix_file_set,
    existing_file_warning,
    target_file_name,
    uuid_filter_chunks,
)
from db2_flattener.gather import lattice
from db2_flattener.schema.constants import MAX_URL_LENGTH

MFS_UUID = "mfs-uuid"


class FakeS3:
    def __init__(self, fail=()):
        self.calls = []
        self.fail = set(fail)

    def get(self, s3_uri, local_path):
        self.calls.append((s3_uri, local_path))
        if s3_uri in self.fail:
            with open(local_path, "w") as fh:
                fh.write("partial")
            raise OSError("boom")
        with open(local_path, "w") as fh:
            fh.write(s3_uri)


def _stub_report(monkeypatch, raw_matrix_files):
    def fake_get_report(obj_type, filter_url, field_lst, connection):
        if obj_type == "MatrixFileSet":
            assert filter_url == f"&uuid={MFS_UUID}"
            return [
                {"raw_matrix_files": [f"/raw_matrix_files/{r['uuid']}/" for r in raw_matrix_files]}
            ]
        assert obj_type == "RawMatrixFile"
        assert "geo_file_name" in field_lst and "s3_uri" in field_lst
        requested = {term.removeprefix("uuid=") for term in filter_url.split("&") if term}
        return [r for r in raw_matrix_files if r["uuid"] in requested]

    monkeypatch.setattr(downloader, "get_report", fake_get_report)


def _many_raw_matrix_files(n):
    return [
        {
            "uuid": str(uuid_mod.UUID(int=i)),
            "@id": f"/raw_matrix_files/{uuid_mod.UUID(int=i)}/",
            "s3_uri": f"s3://b/o/file{i:04d}.h5",
        }
        for i in range(n)
    ]


def test_many_raw_matrix_files_are_fetched_in_bounded_filtered_requests(monkeypatch):
    """Through the real get_report: every request is filtered and short enough."""
    files = _many_raw_matrix_files(200)
    sent = []

    def fake_requests_get(url, auth=None):
        sent.append(url)
        query = parse_qs(urlparse(url).query)
        if query["type"] == ["MatrixFileSet"]:
            graph = [{"raw_matrix_files": [f["@id"] for f in files]}]
        else:
            requested = set(query.get("uuid", []))
            # An unfiltered request would get every RawMatrixFile in the database.
            graph = (
                [f for f in files if f["uuid"] in requested]
                if requested
                else files
                + [
                    {
                        "uuid": "stranger",
                        "@id": "/raw_matrix_files/stranger/",
                        "s3_uri": "s3://x/y.h5",
                    }
                ]
            )
        response = mock.Mock()
        response.json.return_value = {"@graph": graph}
        return response

    class Conn:
        server = "https://lattice.test/"
        auth = None

    monkeypatch.setattr(lattice.requests, "get", fake_requests_get)
    result = downloader.fetch_raw_matrix_files(MFS_UUID, Conn())

    rmf_urls = [u for u in sent if "type=RawMatrixFile" in u]
    assert len(rmf_urls) > 1
    for url in rmf_urls:
        assert len(url) <= MAX_URL_LENGTH
        assert "uuid=" in url
        assert "@id=" not in url
    requested = [u for url in rmf_urls for u in parse_qs(urlparse(url).query)["uuid"]]
    assert sorted(requested) == sorted(f["uuid"] for f in files)
    assert sorted(r["uuid"] for r in result) == sorted(f["uuid"] for f in files)


def test_uuid_filter_chunks_respect_limit_and_keep_every_uuid():
    uuids = [str(uuid_mod.UUID(int=i)) for i in range(200)]
    chunks = uuid_filter_chunks(uuids, max_length=500)
    assert len(chunks) > 1
    assert all(0 < len(c) <= 500 for c in chunks)
    assert "".join(chunks) == "".join(f"&uuid={u}" for u in uuids)
    assert uuid_filter_chunks([]) == []


def test_unrequested_raw_matrix_files_are_never_returned(monkeypatch):
    listed = _many_raw_matrix_files(2)
    stranger = {"uuid": "stranger", "@id": "/raw_matrix_files/stranger/", "s3_uri": "s3://x/y.h5"}

    def fake_get_report(obj_type, filter_url, field_lst, connection):
        if obj_type == "MatrixFileSet":
            return [{"raw_matrix_files": [f["@id"] for f in listed]}]
        return [*listed, stranger]

    monkeypatch.setattr(downloader, "get_report", fake_get_report)
    result = downloader.fetch_raw_matrix_files(MFS_UUID, None)
    assert [r["uuid"] for r in result] == [f["uuid"] for f in listed]


def test_target_file_name_prefers_geo_file_name():
    assert target_file_name({"s3_uri": "s3://b/x/orig.h5", "geo_file_name": "GSM1.h5"}) == "GSM1.h5"
    assert target_file_name({"s3_uri": "s3://b/x/orig.h5"}) == "orig.h5"
    assert target_file_name({"s3_uri": "s3://b/x/orig.h5", "geo_file_name": ""}) == "orig.h5"


def test_existing_file_warning_only_names_geo_file_name_when_present():
    assert (
        existing_file_warning({"s3_uri": "s3://b/a.h5", "geo_file_name": "GSM1.h5"})
        == "====== WARNING ====== For s3://b/a.h5, there already exists a file GSM1.h5"
    )
    assert (
        existing_file_warning({"s3_uri": "s3://b/a.h5"})
        == "====== WARNING ====== For s3://b/a.h5, there already exists a file"
    )


def test_downloads_with_geo_name_or_s3_name(monkeypatch, tmp_path):
    _stub_report(
        monkeypatch,
        [
            {
                "uuid": "1",
                "@id": "/raw_matrix_files/1/",
                "s3_uri": "s3://b/o/one.h5",
                "geo_file_name": "GSM_one.h5",
            },
            {"uuid": "2", "@id": "/raw_matrix_files/2/", "s3_uri": "s3://b/o/two.h5"},
        ],
    )
    out = tmp_path / "new" / "dir"
    fs = FakeS3()
    downloaded, skipped, failed = download_matrix_file_set(MFS_UUID, str(out), None, fs=fs)

    assert downloaded == ["s3://b/o/one.h5", "s3://b/o/two.h5"]
    assert skipped == failed == []
    assert sorted(p.name for p in out.iterdir()) == ["GSM_one.h5", "two.h5"]
    assert (out / "GSM_one.h5").read_text() == "s3://b/o/one.h5"


def test_existing_file_is_not_downloaded_and_warns(monkeypatch, tmp_path, capsys):
    _stub_report(
        monkeypatch,
        [
            {
                "uuid": "1",
                "@id": "/raw_matrix_files/1/",
                "s3_uri": "s3://b/o/one.h5",
                "geo_file_name": "GSM_one.h5",
            },
            {"uuid": "2", "@id": "/raw_matrix_files/2/", "s3_uri": "s3://b/o/two.h5"},
        ],
    )
    (tmp_path / "GSM_one.h5").write_text("old")
    (tmp_path / "two.h5").write_text("old")
    fs = FakeS3()
    downloaded, skipped, failed = download_matrix_file_set(MFS_UUID, str(tmp_path), None, fs=fs)

    assert fs.calls == []
    assert downloaded == []
    assert skipped == ["s3://b/o/one.h5", "s3://b/o/two.h5"]
    assert (tmp_path / "GSM_one.h5").read_text() == "old"
    out = capsys.readouterr().out
    assert (
        "====== WARNING ====== For s3://b/o/one.h5, there already exists a file GSM_one.h5\n" in out
    )
    assert "====== WARNING ====== For s3://b/o/two.h5, there already exists a file\n" in out


def _rmf(uuid, s3_uri, geo_file_name=None):
    record = {"uuid": uuid, "@id": f"/raw_matrix_files/{uuid}/", "s3_uri": s3_uri}
    if geo_file_name is not None:
        record["geo_file_name"] = geo_file_name
    return record


@pytest.mark.parametrize(
    ("clashing", "saved_as"),
    [
        pytest.param(
            [
                _rmf("1", "s3://b/s1/outs/filtered_feature_bc_matrix.h5"),
                _rmf("2", "s3://b/s2/outs/filtered_feature_bc_matrix.h5"),
            ],
            "filtered_feature_bc_matrix.h5",
            id="same-s3-basename",
        ),
        pytest.param(
            [
                _rmf("1", "s3://b/s1/a.h5", "GSM1.h5"),
                _rmf("2", "s3://b/s2/b.h5", "GSM1.h5"),
            ],
            "GSM1.h5",
            id="same-geo-file-name",
        ),
        pytest.param(
            [
                _rmf("1", "s3://b/s1/a.h5", "GSM1.h5"),
                _rmf("2", "s3://b/s2/b.h5", "gsm1.h5"),
            ],
            "GSM1.h5 / gsm1.h5",
            id="differ-only-in-case",
        ),
        pytest.param(
            [
                _rmf("1", "s3://b/s1/GSM1.h5"),
                _rmf("2", "s3://b/s2/b.h5", "GSM1.h5"),
            ],
            "GSM1.h5",
            id="s3-basename-vs-geo-file-name",
        ),
    ],
)
def test_colliding_target_names_are_all_refused(monkeypatch, tmp_path, capsys, clashing, saved_as):
    ok = _rmf("9", "s3://b/s9/other.h5", "GSM9.h5")
    _stub_report(monkeypatch, [*clashing, ok])
    fs = FakeS3()

    downloaded, skipped, failed = download_matrix_file_set(MFS_UUID, str(tmp_path), None, fs=fs)

    assert downloaded == [ok["s3_uri"]]
    assert skipped == []
    assert sorted(failed) == sorted(r["s3_uri"] for r in clashing)
    assert [call[0] for call in fs.calls] == [ok["s3_uri"]]
    assert [p.name for p in tmp_path.iterdir()] == ["GSM9.h5"]
    out = capsys.readouterr().out
    assert f"2 raw matrix files would be saved as {saved_as}; none of them downloaded" in out
    for record in clashing:
        assert f"    {record['s3_uri']}\n" in out
    assert "already exists" not in out


def test_collision_is_refused_even_when_file_already_exists(monkeypatch, tmp_path):
    """A clash is a metadata problem; a file left from an earlier run doesn't settle it."""
    clashing = [
        _rmf("1", "s3://b/s1/a.h5", "GSM1.h5"),
        _rmf("2", "s3://b/s2/b.h5", "GSM1.h5"),
    ]
    _stub_report(monkeypatch, clashing)
    (tmp_path / "GSM1.h5").write_text("old")

    downloaded, skipped, failed = download_matrix_file_set(
        MFS_UUID, str(tmp_path), None, fs=FakeS3()
    )

    assert downloaded == skipped == []
    assert sorted(failed) == ["s3://b/s1/a.h5", "s3://b/s2/b.h5"]
    assert (tmp_path / "GSM1.h5").read_text() == "old"


def test_failed_download_leaves_no_partial_file(monkeypatch, tmp_path):
    _stub_report(
        monkeypatch,
        [{"uuid": "1", "@id": "/raw_matrix_files/1/", "s3_uri": "s3://b/o/one.h5"}],
    )
    downloaded, skipped, failed = download_matrix_file_set(
        MFS_UUID, str(tmp_path), None, fs=FakeS3(fail={"s3://b/o/one.h5"})
    )
    assert failed == ["s3://b/o/one.h5"]
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("bad_name", ["../escape.h5", "sub/dir.h5", ".."])
def test_geo_file_name_with_path_is_refused(monkeypatch, tmp_path, bad_name):
    _stub_report(
        monkeypatch,
        [
            {
                "uuid": "1",
                "@id": "/raw_matrix_files/1/",
                "s3_uri": "s3://b/o/one.h5",
                "geo_file_name": bad_name,
            }
        ],
    )
    fs = FakeS3()
    _, _, failed = download_matrix_file_set(MFS_UUID, str(tmp_path / "out"), None, fs=fs)
    assert failed == ["s3://b/o/one.h5"]
    assert fs.calls == []


def test_missing_s3_uri_counts_as_failed(monkeypatch, tmp_path):
    _stub_report(monkeypatch, [{"uuid": "1", "@id": "/raw_matrix_files/1/"}])
    _, _, failed = download_matrix_file_set(MFS_UUID, str(tmp_path), None, fs=FakeS3())
    assert failed == ["/raw_matrix_files/1/"]


def test_missing_matrix_file_set_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(downloader, "get_report", lambda *a, **k: [])
    with pytest.raises(ValueError, match="not found"):
        download_matrix_file_set(MFS_UUID, str(tmp_path), None, fs=FakeS3())


def _run_cli(monkeypatch, result):
    monkeypatch.setattr(cli, "Connection", lambda mode: object())
    seen = {}

    def fake_download(uuid, out_directory, connection):
        seen.update(uuid=uuid, out_directory=out_directory)
        return result

    monkeypatch.setattr(cli, "download_matrix_file_set", fake_download)
    monkeypatch.setattr(
        "sys.argv", ["db2-flattener-download", "-u", MFS_UUID, "--out-directory", "outdir"]
    )
    return seen


def test_cli_passes_args_and_succeeds(monkeypatch):
    seen = _run_cli(monkeypatch, (["a"], ["b"], []))
    cli.main()
    assert seen == {"uuid": MFS_UUID, "out_directory": "outdir"}


def test_cli_exits_nonzero_on_failures(monkeypatch):
    _run_cli(monkeypatch, ([], [], ["a"]))
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 1


def test_cli_requires_out_directory(monkeypatch):
    monkeypatch.setattr("sys.argv", ["db2-flattener-download", "-u", MFS_UUID])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2
