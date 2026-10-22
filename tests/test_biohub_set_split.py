"""BIOHUB samples/guide templates split by GeneticModification set."""

from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import Configs

UUID = "9b1c06e5-aaaa-bbbb-cccc-dddddddddddd"
TIMESTAMP = "20260814_095000"
GM_AAA = "/genetic_modifications/aaa/"
GM_BBB = "/genetic_modifications/bbb/"
LIB_A = "/droplet_based_libraries/lib-a/"
LIB_B = "/droplet_based_libraries/lib-b/"
GUIDE_A = {"@id": "/tabular_files/guide-a/", "s3_uri": "s3://bucket/a.tsv", "file_format": "tsv"}
GUIDE_B = {"@id": "/tabular_files/guide-b/", "s3_uri": "s3://bucket/b.tsv", "file_format": "tsv"}


def _freeze_time(monkeypatch):
    class FrozenDateTime:
        @staticmethod
        def now():
            return datetime(2026, 8, 14, 9, 50, 0)

    monkeypatch.setattr("db2_flattener.flatten.flattener.datetime", FrozenDateTime)


def _base_flattener(complete_data, main_df, sample_df=None):
    if sample_df is None:
        sample_df = pd.DataFrame()
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f.create_dataframe = lambda data: (main_df.copy(), sample_df.copy())
    f.create_geo_dataframe = lambda main: pd.DataFrame({"library_name": ["lib1"]})
    f.create_sra_biosample_dataframe = lambda main: pd.DataFrame({"sample_name": ["s1"]})
    f.create_sequence_file_set_dataframe = lambda data: pd.DataFrame()
    f.create_sra_files_dataframe = lambda main, sequence_file_set_df=None: pd.DataFrame(
        {"library_ID": ["lib1"]}
    )
    f._attach_plate_library_id = lambda main, sample: sample
    f.create_samples_dataframe = lambda sample: sample
    return f


def test_write_splits_samples_and_guides_by_gm_set(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    main_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1", "rmf2"],
            "raw_file_samples": ["s1", "s2"],
            "droplet_based_libraries_@id": [LIB_A, LIB_B],
            "sequence_file_sets_CRO_order": ["ABC-123", "ABC-123"],
        }
    )
    complete_data = {
        "libraries": {
            "uuid-a": {
                "library": {"@id": LIB_A},
                "genetic_modifications": [{"@id": GM_AAA, "guide_rna_files": GUIDE_A}],
            },
            "uuid-b": {
                "library": {"@id": LIB_B},
                "genetic_modifications": [{"@id": GM_BBB, "guide_rna_files": GUIDE_B}],
            },
        },
        "resolved_objects": {
            "GeneticModification": {
                GM_AAA: {"@id": GM_AAA, "guide_rna_files": GUIDE_A},
                GM_BBB: {"@id": GM_BBB, "guide_rna_files": GUIDE_B},
            },
            "TabularFile": {GUIDE_A["@id"]: GUIDE_A, GUIDE_B["@id"]: GUIDE_B},
        },
    }

    flattener = _base_flattener(complete_data, main_df)
    flattener.create_biohub_dataframe = lambda source: pd.DataFrame(
        {
            "sample_name": source["raw_file_samples"].tolist()
            if "raw_file_samples" in source.columns and not source.empty
            else []
        }
    )
    flattener.create_guide_metadata_dataframe = lambda file_info: pd.DataFrame(
        {"id": [file_info["@id"].rsplit("/", 2)[-2]]}
    )

    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")
    outdir = Path(result)
    biohub = outdir / "BIOHUB"

    assert (biohub / "samples_template_set1_v1.csv").is_file()
    assert (biohub / "samples_template_set2_v1.csv").is_file()
    assert (biohub / "guide_template_set1_v1.csv").is_file()
    assert (biohub / "guide_template_set2_v1.csv").is_file()
    assert not (biohub / "samples_template_v1.csv").exists()

    set1 = pd.read_csv(biohub / "samples_template_set1_v1.csv")
    set2 = pd.read_csv(biohub / "samples_template_set2_v1.csv")
    assert list(set1["sample_name"]) == ["s1"]
    assert list(set2["sample_name"]) == ["s2"]

    guide1 = pd.read_csv(biohub / "guide_template_set1_v1.csv")
    guide2 = pd.read_csv(biohub / "guide_template_set2_v1.csv")
    assert list(guide1["id"]) == ["guide-a"]
    assert list(guide2["id"]) == ["guide-b"]


def test_no_gms_writes_samples_template_v1(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    main_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1"],
            "raw_file_samples": ["s1"],
            "sequence_file_sets_CRO_order": ["ABC-123"],
        }
    )
    complete_data = {"libraries": {}, "resolved_objects": {}}
    flattener = _base_flattener(complete_data, main_df)
    flattener.create_biohub_dataframe = lambda source: pd.DataFrame({"sample_name": ["s1"]})

    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")
    biohub = Path(result) / "BIOHUB"
    assert (biohub / "samples_template_v1.csv").is_file()
    assert not list(biohub.glob("samples_template_set*_v1.csv"))
    assert not list(biohub.glob("guide_template_set*_v1.csv"))


def test_unassigned_rows_warn(tmp_path, monkeypatch, capsys):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    main_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1", "rmf2"],
            "raw_file_samples": ["s1", "s2"],
            "droplet_based_libraries_@id": [LIB_A, "/droplet_based_libraries/orphan/"],
            "sequence_file_sets_CRO_order": ["ABC-123", "ABC-123"],
        }
    )
    complete_data = {
        "libraries": {
            "uuid-a": {
                "library": {"@id": LIB_A},
                "genetic_modifications": [{"@id": GM_AAA}],
            },
        },
        "resolved_objects": {
            "GeneticModification": {
                GM_AAA: {"@id": GM_AAA, "guide_rna_files": GUIDE_A},
            },
            "TabularFile": {GUIDE_A["@id"]: GUIDE_A},
        },
    }
    flattener = _base_flattener(complete_data, main_df)
    flattener.create_biohub_dataframe = lambda source: pd.DataFrame(
        {
            "sample_name": source["raw_file_samples"].tolist()
            if "raw_file_samples" in source.columns and not source.empty
            else []
        }
    )
    flattener.create_guide_metadata_dataframe = lambda file_info: pd.DataFrame({"id": ["g1"]})

    flattener.flatten_matrix_file_set(UUID, biohub_version="v1")
    assert "1 BIOHUB source row(s) have no matching GeneticModification set" in capsys.readouterr().out
