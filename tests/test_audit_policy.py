"""Audit category collection, printing, and CLI policy enforcement."""

from types import SimpleNamespace

import pytest

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import Configs


def make_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = Configs(
        FIELD_TYPES={},
        OBJECT_CONFIG={
            "droplet_based_libraries": {
                "api_type": "DropletBasedLibrary",
                "fields": ["@id", "audit"],
                "references": {},
            },
            "tissues": {
                "api_type": "Tissue",
                "fields": ["@id", "audit"],
                "references": {},
            },
            "sequence_files": {
                "api_type": "SequenceFile",
                "fields": ["@id", "audit"],
                "references": {},
            },
        },
    )
    f.gatherer = None
    return f


def test_collect_audit_category_counts_by_type_and_level():
    flattener = make_flattener()
    complete_data = {
        "resolved_objects": {
            "Tissue": {
                "/tissues/t1/": {
                    "@id": "/tissues/t1/",
                    "audit": {
                        "ERROR": [
                            {"category": "missing something"},
                            {"category": "missing something"},
                        ],
                        "WARNING": [{"category": "odd barcode"}],
                        "INTERNAL_ACTION": [{"category": "ignored"}],
                    },
                }
            },
            "DropletBasedLibrary": {
                "/droplet_based_libraries/l1/": {
                    "@id": "/droplet_based_libraries/l1/",
                    "audit": {"WARNING": [{"category": "odd barcode"}]},
                }
            },
        },
        "libraries": {},
    }

    counts = flattener._collect_audit_category_counts(complete_data)

    assert counts["ERROR"] == {("Tissue", "missing something"): 2}
    assert counts["WARNING"] == {
        ("DropletBasedLibrary", "odd barcode"): 1,
        ("Tissue", "odd barcode"): 1,
    }


def test_collect_includes_library_not_in_resolved_objects():
    flattener = make_flattener()
    complete_data = {
        "resolved_objects": {},
        "libraries": {
            "uuid-a": {
                "library": {
                    "@id": "/droplet_based_libraries/l1/",
                    "@type": ["DropletBasedLibrary"],
                    "audit": {"ERROR": [{"category": "lib-only"}]},
                }
            }
        },
    }

    counts = flattener._collect_audit_category_counts(complete_data)

    assert counts["ERROR"] == {("DropletBasedLibrary", "lib-only"): 1}


def test_collect_includes_nested_objects_under_libraries():
    """Objects only attached under library buckets (not resolved_objects) are assessed."""
    flattener = make_flattener()
    complete_data = {
        "resolved_objects": {},
        "libraries": {
            "uuid-a": {
                "library": {
                    "@id": "/droplet_based_libraries/l1/",
                    "@type": ["DropletBasedLibrary"],
                },
                "samples": [
                    {
                        "@id": "/tissues/t-nested/",
                        "@type": ["Tissue"],
                        "audit": {"WARNING": [{"category": "nested sample"}]},
                    }
                ],
                "raw_matrix_files": [
                    {
                        "@id": "/raw_matrix_files/r1/",
                        "sequence_files": [
                            {
                                "@id": "/sequence_files/s1/",
                                "@type": ["SequenceFile"],
                                "audit": {"ERROR": [{"category": "nested seq"}]},
                            }
                        ],
                    }
                ],
            }
        },
    }

    counts = flattener._collect_audit_category_counts(complete_data)

    assert counts["WARNING"] == {("Tissue", "nested sample"): 1}
    assert counts["ERROR"] == {("SequenceFile", "nested seq"): 1}


def test_collect_dedupes_object_in_resolved_and_libraries():
    flattener = make_flattener()
    sample = {
        "@id": "/tissues/t1/",
        "@type": ["Tissue"],
        "audit": {"ERROR": [{"category": "once"}]},
    }
    complete_data = {
        "resolved_objects": {"Tissue": {"/tissues/t1/": sample}},
        "libraries": {
            "uuid-a": {"library": {"@id": "/droplet_based_libraries/l1/"}, "samples": [sample]}
        },
    }

    counts = flattener._collect_audit_category_counts(complete_data)

    assert counts["ERROR"] == {("Tissue", "once"): 1}


def test_print_audit_category_counts_only_levels_with_rows(capsys):
    DB2Flattener._print_audit_category_counts(
        {
            "ERROR": {},
            "WARNING": {("Tissue", "odd barcode"): 2},
        }
    )
    out = capsys.readouterr().out
    assert "ERROR audits" not in out
    assert out.splitlines() == [
        "======",
        "WARNING audits:",
        "  Tissue | odd barcode: 2",
        "======",
    ]


@pytest.mark.parametrize(
    "audit_mode,counts,should_raise",
    [
        ("strict", {"ERROR": {("T", "a"): 1}, "WARNING": {}}, True),
        ("strict", {"ERROR": {}, "WARNING": {("T", "a"): 1}}, True),
        ("strict", {"ERROR": {}, "WARNING": {}}, False),
        ("warnings_ok", {"ERROR": {("T", "a"): 1}, "WARNING": {}}, True),
        ("warnings_ok", {"ERROR": {}, "WARNING": {("T", "a"): 1}}, False),
        ("all_ok", {"ERROR": {("T", "a"): 1}, "WARNING": {("T", "b"): 1}}, False),
    ],
)
def test_enforce_audit_policy(audit_mode, counts, should_raise):
    if should_raise:
        with pytest.raises(ValueError, match="Refusing to flatten"):
            DB2Flattener._enforce_audit_policy(counts, audit_mode)
    else:
        DB2Flattener._enforce_audit_policy(counts, audit_mode)


def test_flatten_stops_on_audits_before_outputs(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)

    complete_data = {
        "resolved_objects": {
            "GeneticModification": {},
            "Tissue": {
                "/tissues/t1/": {
                    "@id": "/tissues/t1/",
                    "audit": {"WARNING": [{"category": "odd barcode"}]},
                }
            },
        },
        "libraries": {},
    }
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None

    with pytest.raises(ValueError, match="audits present"):
        f.flatten_matrix_file_set("abc", biohub_version="v1", audit_mode="strict")

    out = capsys.readouterr().out
    assert "WARNING audits:" in out
    assert "Tissue | odd barcode: 1" in out
    assert not any(tmp_path.iterdir())


def test_sample_only_tissue_error_stops_strict_flatten(tmp_path, monkeypatch, capsys):
    """Regression: sample audits live under libraries.samples, not resolved_objects."""
    monkeypatch.chdir(tmp_path)

    complete_data = {
        "resolved_objects": {
            "GeneticModification": {},
            # Gatherer does not put samples here; only library references / files.
        },
        "libraries": {
            "uuid-a": {
                "library": {
                    "@id": "/droplet_based_libraries/l1/",
                    "@type": ["DropletBasedLibrary"],
                },
                "samples": [
                    {
                        "@id": "/tissues/82d7ef44-5d12-4fd6-a67a-4006b04d69a2/",
                        "@type": ["Tissue"],
                        "audit": {
                            "ERROR": [
                                {
                                    "path": "/tissues/82d7ef44-5d12-4fd6-a67a-4006b04d69a2/",
                                    "level_name": "ERROR",
                                    "level": 60,
                                    "name": "audit_biosample_object_dispatcher",
                                    "detail": "Tissue is missing developmental_stages.",
                                    "category": "missing developmental stage",
                                }
                            ]
                        },
                    }
                ],
            }
        },
    }
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None

    counts = f._collect_audit_category_counts(complete_data)
    assert counts["ERROR"] == {("Tissue", "missing developmental stage"): 1}

    with pytest.raises(ValueError, match="audits present"):
        f.flatten_matrix_file_set("abc", biohub_version="v1", audit_mode="strict")

    out = capsys.readouterr().out
    assert "ERROR audits:" in out
    assert "Tissue | missing developmental stage: 1" in out
    assert not any(tmp_path.iterdir())


def test_stub_before_full_sample_audit_still_stops_strict_flatten(tmp_path, monkeypatch, capsys):
    """Regression: link stubs without audit must not hide the full sample audit.

    RawMatrixFile / library.samples often embed ``{@id: ...}`` stubs that are
    walked before ``libraries[*].samples`` holds the fetched Tissue with audit.
    First-wins ``@id`` dedupe would skip the audited copy and let flatten proceed.
    """
    monkeypatch.chdir(tmp_path)

    tissue_id = "/tissues/82d7ef44-5d12-4fd6-a67a-4006b04d69a2/"
    full_tissue = {
        "@id": tissue_id,
        "@type": ["Tissue"],
        "audit": {
            "ERROR": [
                {
                    "path": tissue_id,
                    "level_name": "ERROR",
                    "level": 60,
                    "name": "audit_biosample_object_dispatcher",
                    "detail": "Tissue is missing developmental_stages.",
                    "category": "missing developmental stage",
                }
            ]
        },
    }
    stub = {"@id": tissue_id}

    complete_data = {
        "resolved_objects": {
            "GeneticModification": {},
            "RawMatrixFile": {
                "/raw_matrix_files/r1/": {
                    "@id": "/raw_matrix_files/r1/",
                    "samples": [stub],
                }
            },
        },
        "libraries": {
            "uuid-a": {
                "library": {
                    "@id": "/droplet_based_libraries/l1/",
                    "@type": ["DropletBasedLibrary"],
                    "samples": [stub],
                },
                "samples": [full_tissue],
            }
        },
    }
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None

    counts = f._collect_audit_category_counts(complete_data)
    assert counts["ERROR"] == {("Tissue", "missing developmental stage"): 1}

    with pytest.raises(ValueError, match="audits present"):
        f.flatten_matrix_file_set("abc", biohub_version="v1", audit_mode="strict")

    out = capsys.readouterr().out
    assert "ERROR audits:" in out
    assert "Tissue | missing developmental stage: 1" in out
    assert not any(tmp_path.iterdir())
