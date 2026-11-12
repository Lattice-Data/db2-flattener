"""Independent hard stop when any gathered object has status 'deleted'."""

from types import SimpleNamespace

import pytest

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.schema.constants import REQUIRED_FETCH_FIELDS, Configs

TISSUE_UUID = "82d7ef44-5d12-4fd6-a67a-4006b04d69a2"
TISSUE_ID = f"/tissues/{TISSUE_UUID}/"


def make_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = Configs(
        FIELD_TYPES={},
        OBJECT_CONFIG={
            "tissues": {"api_type": "Tissue", "fields": ["@id", "aliases"], "references": {}},
            "sequence_files": {
                "api_type": "SequenceFile",
                "fields": ["@id"],
                "references": {},
            },
            "droplet_based_libraries": {
                "api_type": "DropletBasedLibrary",
                "fields": ["@id"],
                "references": {},
            },
        },
    )
    f.gatherer = None
    return f


def deleted_tissue(**extra):
    return {
        "@id": TISSUE_ID,
        "@type": ["Tissue"],
        "status": "deleted",
        "aliases": ["lab:t1", "lab:t1-alt"],
        **extra,
    }


def test_no_deleted_objects_passes():
    f = make_flattener()
    data = {
        "resolved_objects": {"Tissue": {TISSUE_ID: {"@id": TISSUE_ID, "status": "released"}}},
        "libraries": {},
    }
    assert f._collect_deleted_objects(data) == []
    f._enforce_no_deleted_objects(data)


def test_deleted_in_resolved_objects_with_all_raw_aliases():
    f = make_flattener()
    data = {"resolved_objects": {"Tissue": {TISSUE_ID: deleted_tissue()}}, "libraries": {}}
    assert f._collect_deleted_objects(data) == [("Tissue", TISSUE_UUID, ["lab:t1", "lab:t1-alt"])]


def test_deleted_nested_under_libraries():
    f = make_flattener()
    seq = {
        "@id": "/sequence_files/sf1/",
        "@type": ["SequenceFile"],
        "status": "deleted",
        "aliases": ["lab:sf1"],
    }
    data = {
        "resolved_objects": {},
        "libraries": {
            "uuid-a": {
                "library": {"@id": "/droplet_based_libraries/l1/", "status": "in progress"},
                "samples": [deleted_tissue()],
                "raw_matrix_files": [{"@id": "/raw_matrix_files/r1/", "sequence_files": [seq]}],
            }
        },
    }
    assert f._collect_deleted_objects(data) == [
        ("SequenceFile", "sf1", ["lab:sf1"]),
        ("Tissue", TISSUE_UUID, ["lab:t1", "lab:t1-alt"]),
    ]


def test_deleted_controlled_term_reported():
    f = make_flattener()
    term_id = "/controlled_terms/EFO:0001/"
    data = {
        "resolved_objects": {
            "ControlledTerm": {
                term_id: {"@id": term_id, "status": "deleted", "aliases": ["lab:term"]}
            }
        },
        "libraries": {},
    }
    assert f._collect_deleted_objects(data) == [("ControlledTerm", "EFO:0001", ["lab:term"])]


def test_deleted_in_matrix_file_set():
    f = make_flattener()
    data = {
        "matrix_file_set": {"@id": "/matrix_file_sets/m1/", "status": "deleted"},
        "resolved_objects": {},
        "libraries": {},
    }
    assert f._collect_deleted_objects(data) == [("Unknown", "m1", [])]


def test_deduped_by_id_and_stub_ignored():
    f = make_flattener()
    stub = {"@id": TISSUE_ID}
    data = {
        "resolved_objects": {
            "Tissue": {TISSUE_ID: deleted_tissue()},
            "RawMatrixFile": {
                "/raw_matrix_files/r1/": {"@id": "/raw_matrix_files/r1/", "samples": [stub]}
            },
        },
        "libraries": {"uuid-a": {"samples": [deleted_tissue()]}},
    }
    assert f._collect_deleted_objects(data) == [("Tissue", TISSUE_UUID, ["lab:t1", "lab:t1-alt"])]


def test_print_and_raise(capsys):
    f = make_flattener()
    data = {
        "resolved_objects": {
            "Tissue": {TISSUE_ID: deleted_tissue()},
            "SequenceFile": {
                "/sequence_files/sf1/": {"@id": "/sequence_files/sf1/", "status": "deleted"}
            },
        },
        "libraries": {},
    }
    with pytest.raises(ValueError, match="2 gathered object"):
        f._enforce_no_deleted_objects(data)
    assert capsys.readouterr().out.splitlines() == [
        "======",
        "DELETED objects:",
        "  SequenceFile | sf1 | (no aliases)",
        f"  Tissue | {TISSUE_UUID} | lab:t1, lab:t1-alt",
        "======",
    ]


@pytest.mark.parametrize("audit_mode", ["strict", "warnings_ok", "all_ok"])
def test_flatten_stops_before_outputs_regardless_of_audit_mode(
    tmp_path, monkeypatch, capsys, audit_mode
):
    monkeypatch.chdir(tmp_path)
    data = {"resolved_objects": {"Tissue": {TISSUE_ID: deleted_tissue()}}, "libraries": {}}
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: data)
    f._validate_genetic_modification_aliases = lambda d: None

    with pytest.raises(ValueError, match="status 'deleted'"):
        f.flatten_matrix_file_set("abc", audit_mode=audit_mode)

    out = capsys.readouterr().out
    assert f"Tissue | {TISSUE_UUID} | lab:t1, lab:t1-alt" in out
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("configured", [["@id"], ["@id", "status"], ["uuid", "aliases", "status"]])
def test_chunk_and_fetch_always_requests_status_and_aliases(monkeypatch, configured):
    seen = {}

    def fake_get_report(obj_type, filter_url, field_lst, connection):
        seen["fields"] = field_lst
        return []

    monkeypatch.setattr("db2_flattener.gather.gatherer.DB2lattice.get_report", fake_get_report)
    configs = Configs(
        FIELD_TYPES={},
        OBJECT_CONFIG={"tissues": {"api_type": "Tissue", "fields": configured, "references": {}}},
    )
    DB2Gatherer(None, configs).chunk_and_fetch("Tissue", ["u1"])

    for field in REQUIRED_FETCH_FIELDS:
        assert seen["fields"].count(field) == 1
    assert configured == seen["fields"][: len(configured)]
    assert configs.OBJECT_CONFIG["tissues"]["fields"] == configured  # config not mutated
