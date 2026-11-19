"""Hard stop when two RawMatrixFiles share a geo_file_name (case-insensitively)."""

from types import SimpleNamespace

import pytest
from test_audit_policy import _complete_data_for_write, make_flattener


def rmf(uuid, geo_file_name, s3_uri="auto"):
    return {
        "@id": f"/raw_matrix_files/{uuid}/",
        "uuid": uuid,
        "geo_file_name": geo_file_name,
        "s3_uri": f"s3://bucket/{uuid}/matrix.h5" if s3_uri == "auto" else s3_uri,
    }


def data_with(*files):
    return {
        "resolved_objects": {"RawMatrixFile": {f["@id"]: f for f in files}},
        "libraries": {},
    }


def test_unique_names_pass(capsys):
    f = make_flattener()
    data = data_with(rmf("u1", "a.h5ad"), rmf("u2", "b.h5ad"))

    assert f._collect_duplicate_geo_file_names(data) == []
    f._enforce_unique_geo_file_names(data)
    assert "ERROR" not in capsys.readouterr().out


def test_identical_names_fail_and_list_both_files(capsys):
    f = make_flattener()
    data = data_with(rmf("u2", "GSM1.h5ad"), rmf("u1", "GSM1.h5ad"), rmf("u3", "other.h5ad"))

    with pytest.raises(ValueError, match="1 geo_file_name"):
        f._enforce_unique_geo_file_names(data)

    assert capsys.readouterr().out.splitlines() == [
        "====== ERROR ====== 2 raw matrix files share geo_file_name GSM1.h5ad:",
        "    /raw_matrix_files/u1/ | s3://bucket/u1/matrix.h5",
        "    /raw_matrix_files/u2/ | s3://bucket/u2/matrix.h5",
    ]


def test_names_differing_only_in_case_collide_and_show_each_spelling(capsys):
    f = make_flattener()
    data = data_with(rmf("u1", "GSM1.h5"), rmf("u2", "gsm1.h5"))

    with pytest.raises(ValueError):
        f._enforce_unique_geo_file_names(data)

    assert (
        "====== ERROR ====== 2 raw matrix files share geo_file_name GSM1.h5 / gsm1.h5:"
        in capsys.readouterr().out
    )


def test_three_files_sharing_a_name_are_one_entry(capsys):
    f = make_flattener()
    data = data_with(rmf("u1", "x.h5"), rmf("u2", "x.h5"), rmf("u3", "x.h5"))

    with pytest.raises(ValueError, match="1 geo_file_name"):
        f._enforce_unique_geo_file_names(data)

    out = capsys.readouterr().out
    assert out.count("====== ERROR ======") == 1
    assert "3 raw matrix files share geo_file_name x.h5:" in out


def test_two_duplicate_groups_are_both_listed(capsys):
    f = make_flattener()
    data = data_with(rmf("u1", "b.h5"), rmf("u2", "b.h5"), rmf("u3", "a.h5"), rmf("u4", "a.h5"))

    with pytest.raises(ValueError, match="2 geo_file_name"):
        f._enforce_unique_geo_file_names(data)

    out = capsys.readouterr().out
    assert out.index("geo_file_name a.h5") < out.index("geo_file_name b.h5")
    assert out.count("====== ERROR ======") == 2


@pytest.mark.parametrize("blank", [None, ""])
def test_blank_names_are_never_duplicates(blank):
    f = make_flattener()
    data = data_with(rmf("u1", blank), rmf("u2", blank), rmf("u3", "a.h5"))

    assert f._collect_duplicate_geo_file_names(data) == []
    f._enforce_unique_geo_file_names(data)


def test_missing_s3_uri_is_shown(capsys):
    f = make_flattener()
    data = data_with(rmf("u1", "x.h5", s3_uri=None), rmf("u2", "x.h5"))

    with pytest.raises(ValueError):
        f._enforce_unique_geo_file_names(data)

    assert "/raw_matrix_files/u1/ | (no s3_uri)" in capsys.readouterr().out


def test_no_raw_matrix_files_passes():
    f = make_flattener()

    for data in ({}, {"resolved_objects": {}, "libraries": {}}, data_with()):
        f._enforce_unique_geo_file_names(data)


@pytest.mark.parametrize("audit_mode", ["strict", "warnings_ok", "all_ok"])
def test_flatten_stops_before_outputs_regardless_of_audit_mode(
    tmp_path, monkeypatch, capsys, audit_mode
):
    monkeypatch.chdir(tmp_path)
    complete_data = _complete_data_for_write(tissue_audit=None, sequence_audit=None)
    complete_data["resolved_objects"]["RawMatrixFile"] = {
        "/raw_matrix_files/u1/": rmf("u1", "GSM1.h5"),
        "/raw_matrix_files/u2/": rmf("u2", "gsm1.h5"),
    }
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None

    with pytest.raises(ValueError, match="shared by more than one raw matrix file"):
        f.flatten_matrix_file_set("abc12345-uuid", audit_mode=audit_mode)

    assert "share geo_file_name GSM1.h5 / gsm1.h5" in capsys.readouterr().out
    assert not any(tmp_path.iterdir())
