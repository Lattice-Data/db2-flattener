"""GEO/md5_checksums.tsv: geo_file_name and md5sum of each RawMatrixFile."""

from types import SimpleNamespace

import pytest
from test_audit_policy import _complete_data_for_write, _stub_downstream_writers, make_flattener

from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.schema.constants import (
    GEO_MD5_CHECKSUM_COL,
    GEO_MD5_FILE_NAME_COL,
    REQUIRED_FETCH_FIELDS,
    REQUIRED_FETCH_FIELDS_BY_TYPE,
    Configs,
)

WARNING = "====== WARNING ======"


def rmf(uuid, geo_file_name=None, md5sum=None, **extra):
    return {
        "@id": f"/raw_matrix_files/{uuid}/",
        "uuid": uuid,
        "geo_file_name": geo_file_name,
        "md5sum": md5sum,
        **extra,
    }


def data_with(*files):
    return {
        "resolved_objects": {"RawMatrixFile": {f["@id"]: f for f in files}},
        "libraries": {},
    }


def test_headers_follow_the_geo_template():
    assert GEO_MD5_FILE_NAME_COL == "file name"
    assert GEO_MD5_CHECKSUM_COL == "file checksum"


def test_rows_are_geo_file_name_and_md5sum_sorted_by_name(capsys):
    f = make_flattener()
    df = f.create_md5_checksums_dataframe(
        data_with(rmf("u2", "b.h5ad", "bbb"), rmf("u1", "a.h5ad", "aaa"))
    )

    assert list(df.columns) == ["file name", "file checksum"]
    assert df.to_dict("records") == [
        {"file name": "a.h5ad", "file checksum": "aaa"},
        {"file name": "b.h5ad", "file checksum": "bbb"},
    ]
    assert WARNING not in capsys.readouterr().out


def test_missing_geo_file_name_is_blank_with_a_warning_and_sorts_last(capsys):
    f = make_flattener()
    df = f.create_md5_checksums_dataframe(
        data_with(rmf("u1", None, "aaa"), rmf("u2", "b.h5ad", "bbb"))
    )

    assert df.to_dict("records") == [
        {"file name": "b.h5ad", "file checksum": "bbb"},
        {"file name": "", "file checksum": "aaa"},
    ]
    out = capsys.readouterr().out
    assert f"{WARNING} RawMatrixFile /raw_matrix_files/u1/ has no geo_file_name; left blank" in out
    assert "u2" not in out


def test_missing_md5sum_is_blank_with_a_warning(capsys):
    f = make_flattener()
    df = f.create_md5_checksums_dataframe(data_with(rmf("u1", "a.h5ad", None)))

    assert df.to_dict("records") == [{"file name": "a.h5ad", "file checksum": ""}]
    assert (
        f"{WARNING} RawMatrixFile /raw_matrix_files/u1/ has no md5sum; left blank"
        in capsys.readouterr().out
    )


def test_both_missing_warns_twice_and_keeps_the_row(capsys):
    f = make_flattener()
    df = f.create_md5_checksums_dataframe(data_with(rmf("u1", "", None)))

    assert len(df) == 1
    out = capsys.readouterr().out
    assert "has no geo_file_name" in out
    assert "has no md5sum" in out


@pytest.mark.parametrize(
    "complete_data",
    [{}, {"resolved_objects": {}, "libraries": {}}, data_with()],
)
def test_no_raw_matrix_files_gives_an_empty_frame(complete_data):
    df = make_flattener().create_md5_checksums_dataframe(complete_data)

    assert df.empty
    assert list(df.columns) == ["file name", "file checksum"]


def test_flatten_writes_md5_checksums_tsv_to_the_geo_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    complete_data = _complete_data_for_write(tissue_audit=None, sequence_audit=None)
    complete_data["resolved_objects"]["RawMatrixFile"] = {
        "/raw_matrix_files/r1/": rmf("r1", "GSM1.h5ad", "abc123"),
    }
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None
    _stub_downstream_writers(f)

    outdir = f.flatten_matrix_file_set("abc12345-uuid", audit_mode="all_ok")

    path = tmp_path / outdir / "GEO" / "md5_checksums.tsv"
    assert path.read_text().splitlines() == ["file name\tfile checksum", "GSM1.h5ad\tabc123"]
    assert not (tmp_path / outdir / "GEO" / "md5_checksums.csv").exists()


def test_flatten_tsv_collapses_tabs_and_quotes_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    complete_data = _complete_data_for_write(tissue_audit=None, sequence_audit=None)
    complete_data["resolved_objects"]["RawMatrixFile"] = {
        "/raw_matrix_files/r1/": rmf("r1", 'GSM1\t"odd",name.h5ad', "abc123"),
    }
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None
    _stub_downstream_writers(f)

    outdir = f.flatten_matrix_file_set("abc12345-uuid", audit_mode="all_ok")

    lines = (tmp_path / outdir / "GEO" / "md5_checksums.tsv").read_text().splitlines()
    assert lines[1] == 'GSM1 "odd",name.h5ad\tabc123'
    assert all(line.count("\t") == 1 for line in lines)


def test_flatten_omits_md5_checksums_tsv_without_raw_matrix_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    complete_data = _complete_data_for_write(tissue_audit=None, sequence_audit=None)
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None
    _stub_downstream_writers(f)

    outdir = f.flatten_matrix_file_set("abc12345-uuid", audit_mode="all_ok")

    assert (tmp_path / outdir / "GEO" / "geo_metadata.csv").is_file()
    assert not (tmp_path / outdir / "GEO" / "md5_checksums.tsv").exists()


def fetched_fields(monkeypatch, api_type, configured):
    seen = {}

    def fake_get_report(obj_type, filter_url, field_lst, connection):
        seen["fields"] = field_lst
        return []

    monkeypatch.setattr("db2_flattener.gather.gatherer.DB2lattice.get_report", fake_get_report)
    configs = Configs(
        FIELD_TYPES={},
        OBJECT_CONFIG={"x": {"api_type": api_type, "fields": configured, "references": {}}},
    )
    DB2Gatherer(None, configs).chunk_and_fetch(api_type, ["u1"])
    return seen["fields"]


@pytest.mark.parametrize("configured", [["@id"], ["@id", "md5sum"], ["geo_file_name", "uuid"]])
def test_raw_matrix_file_fetch_always_requests_geo_file_name_and_md5sum(monkeypatch, configured):
    fields = fetched_fields(monkeypatch, "RawMatrixFile", configured)

    for field in (*REQUIRED_FETCH_FIELDS, "geo_file_name", "md5sum"):
        assert fields.count(field) == 1
    assert fields[: len(configured)] == configured


def test_other_types_do_not_request_the_raw_matrix_file_fields(monkeypatch):
    fields = fetched_fields(monkeypatch, "SequenceFile", ["@id"])

    assert "geo_file_name" not in fields
    assert "md5sum" not in fields
    assert set(REQUIRED_FETCH_FIELDS_BY_TYPE) == {"RawMatrixFile"}
