import pytest

from db2_flattener.download import cli, downloader
from db2_flattener.download.downloader import (
    download_matrix_file_set,
    existing_file_warning,
    target_file_name,
)

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
        return list(raw_matrix_files)

    monkeypatch.setattr(downloader, "get_report", fake_get_report)


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
