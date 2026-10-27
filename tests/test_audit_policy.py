"""Audit category collection, printing, and CLI policy enforcement."""

from types import SimpleNamespace

import pandas as pd
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
                "fields": ["@id", "aliases", "feature_types", "audit"],
                "references": {},
            },
            "tissues": {
                "api_type": "Tissue",
                "fields": ["@id", "aliases", "audit"],
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


def _stub_downstream_writers(flattener):
    """Skip heavy sheet builders so flatten can write MAIN from real create_dataframe."""
    flattener.biohub_source_dataframe = lambda main, samples: main
    flattener.create_biohub_dataframe = lambda main: pd.DataFrame({"sample_name": ["s1"]})
    flattener.create_geo_dataframe = lambda main: pd.DataFrame({"library_name": ["lib1"]})
    flattener.create_sra_biosample_dataframe = lambda main: pd.DataFrame({"sample_name": ["s1"]})
    flattener.create_sequence_file_set_dataframe = lambda complete_data: pd.DataFrame()
    flattener.create_sra_files_dataframe = lambda main, sequence_file_set_df=None: pd.DataFrame(
        {"sample_name": ["s1"], "library_ID": ["lib1"], "file_size": [0]}
    )
    flattener._biohub_guide_sets = lambda complete_data: []
    flattener.create_guide_metadata_dataframe = lambda file_info: None
    flattener._attach_plate_library_id = lambda main, sample: sample
    flattener.create_samples_dataframe = lambda sample: sample


def _complete_data_for_write(*, tissue_audit, sequence_audit):
    tissue_id = "/tissues/t1/"
    return {
        "resolved_objects": {"GeneticModification": {}},
        "libraries": {
            "uuid-a": {
                "library": {
                    "@id": "/droplet_based_libraries/l1/",
                    "@type": ["DropletBasedLibrary"],
                    "aliases": ["lab:gex"],
                    "feature_types": ["Gene Expression"],
                },
                "samples": [
                    {
                        "@id": tissue_id,
                        "@type": ["Tissue"],
                        "aliases": ["lab:sample1"],
                        "audit": tissue_audit,
                    }
                ],
                "raw_matrix_files": [
                    {
                        "@id": "/raw_matrix_files/r1/",
                        "aliases": ["lab:rmf1"],
                        "samples": [tissue_id],
                        "sequence_files": [
                            {
                                "@id": "/sequence_files/sf1/",
                                "audit": sequence_audit,
                            }
                        ],
                    }
                ],
            }
        },
    }


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


def test_flatten_all_ok_writes_sample_and_sequence_audits_to_main(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    tissue_audit = {"ERROR": [{"category": "missing developmental stage"}]}
    sequence_audit = {"WARNING": [{"category": "odd barcode"}]}
    complete_data = _complete_data_for_write(
        tissue_audit=tissue_audit,
        sequence_audit=sequence_audit,
    )

    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None
    _stub_downstream_writers(f)

    outdir = f.flatten_matrix_file_set("abc12345-uuid", biohub_version="v1", audit_mode="all_ok")

    assert outdir is not None
    main_path = tmp_path / outdir / "MAIN.csv"
    assert main_path.is_file()
    main_df = pd.read_csv(main_path)
    assert "tissues_audit" in main_df.columns
    assert "sequence_files_audit" in main_df.columns
    assert "missing developmental stage" in str(main_df["tissues_audit"].iloc[0])
    assert "odd barcode" in str(main_df["sequence_files_audit"].iloc[0])


@pytest.mark.parametrize("empty_audit", [None, {}, ""])
def test_flatten_all_ok_skips_empty_sample_and_sequence_audits(tmp_path, monkeypatch, empty_audit):
    monkeypatch.chdir(tmp_path)

    complete_data = _complete_data_for_write(
        tissue_audit=empty_audit,
        sequence_audit=empty_audit,
    )

    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None
    _stub_downstream_writers(f)

    outdir = f.flatten_matrix_file_set("abc12345-uuid", biohub_version="v1", audit_mode="all_ok")

    assert outdir is not None
    main_df = pd.read_csv(tmp_path / outdir / "MAIN.csv")
    assert "tissues_audit" not in main_df.columns
    assert "sequence_files_audit" not in main_df.columns


def test_flatten_warnings_ok_writes_warning_audits_to_main(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    tissue_audit = {"WARNING": [{"category": "odd barcode"}]}
    sequence_audit = {"WARNING": [{"category": "low read count"}]}
    complete_data = _complete_data_for_write(
        tissue_audit=tissue_audit,
        sequence_audit=sequence_audit,
    )

    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f._validate_genetic_modification_aliases = lambda data: None
    _stub_downstream_writers(f)

    outdir = f.flatten_matrix_file_set(
        "abc12345-uuid", biohub_version="v1", audit_mode="warnings_ok"
    )

    assert outdir is not None
    main_df = pd.read_csv(tmp_path / outdir / "MAIN.csv")
    assert "tissues_audit" in main_df.columns
    assert "sequence_files_audit" in main_df.columns
    assert "odd barcode" in str(main_df["tissues_audit"].iloc[0])
