from datetime import datetime
from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import Configs

UUID = "9b1c06e5-aaaa-bbbb-cccc-dddddddddddd"
TIMESTAMP = "20260814_095000"


def make_flattener(sample_df, *, cro_orders=None):
    n = len(cro_orders) if cro_orders is not None else 1
    main_data = {
        "raw_matrix_file_alias": [f"rmf{i}" for i in range(n)],
        "sample_alias": [f"s{i}" for i in range(n)],
        "droplet_based_libraries_aliases": [["alex-marson:lib1"]] * n,
    }
    if cro_orders is not None:
        main_data["sequence_file_sets_CRO_order"] = cro_orders
    main_df = pd.DataFrame(main_data)

    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: {"libraries": {}})
    f.create_dataframe = lambda complete_data: (main_df.copy(), sample_df.copy())
    f.biohub_source_dataframe = lambda main, samples: main
    f.create_biohub_dataframe = lambda main: pd.DataFrame({"sample_name": ["s1"]})
    f.create_geo_dataframe = lambda main: pd.DataFrame({"library_name": ["lib1"]})
    f.create_sra_biosample_dataframe = lambda main: pd.DataFrame({"sample_name": ["s1"]})
    f.create_sequence_file_set_dataframe = lambda complete_data: pd.DataFrame()
    f.create_sra_files_dataframe = lambda main, sequence_file_set_df=None: pd.DataFrame(
        {"sample_name": ["s1"], "library_ID": ["lib1"], "file_size": [0]}
    )
    # No GeneticModifications → samples_template_v1.csv fallback
    f._biohub_guide_sets = lambda complete_data: []
    f.create_guide_metadata_dataframe = lambda file_info: None
    f._attach_plate_library_id = lambda main, sample: sample
    f.create_samples_dataframe = lambda sample: sample
    return f


def _freeze_time(monkeypatch):
    class FrozenDateTime:
        @staticmethod
        def now():
            return datetime(2026, 8, 14, 9, 50, 0)

    monkeypatch.setattr("db2_flattener.flatten.flattener.datetime", FrozenDateTime)


def test_single_cro_order_names_output_dir(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    flattener = make_flattener(
        pd.DataFrame({"sample_alias": ["s1"]}),
        cro_orders=["ABC-123"],
    )
    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")

    outdir_name = f"ABC-123_{TIMESTAMP}"
    outdir = tmp_path / outdir_name
    assert result == outdir_name
    assert (outdir / "MAIN.csv").is_file()
    assert (outdir / "GEO" / "geo_metadata.csv").is_file()
    assert (outdir / "GEO" / "demultiplexing_samples.csv").is_file()
    assert (outdir / "SRA" / "sra_biosample_set1.csv").is_file()
    assert (outdir / "SRA" / "sra_file_set1.csv").is_file()
    assert (outdir / "BIOHUB" / "samples_template_v1.csv").is_file()


def test_multiple_cro_orders_concatenated_sorted(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    # One cell with joined orders plus a second row with a duplicate; uniques sort to ABC then DEF.
    flattener = make_flattener(
        pd.DataFrame({"sample_alias": ["s1", "s2"]}),
        cro_orders=["DEF-456; ABC-123", "ABC-123"],
    )
    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")

    outdir_name = f"ABC-123_DEF-456_{TIMESTAMP}"
    assert result == outdir_name
    assert (tmp_path / outdir_name).is_dir()


def test_missing_cro_order_falls_back_to_uuid(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    flattener = make_flattener(pd.DataFrame({"sample_alias": ["s1"]}))
    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")

    outdir_name = f"MatrixFileSet_{UUID[:8]}_{TIMESTAMP}"
    outdir = tmp_path / outdir_name
    assert result == outdir_name
    assert (outdir / "MAIN.csv").is_file()
    assert (outdir / "GEO" / "geo_metadata.csv").is_file()
    assert (outdir / "SRA" / "sra_biosample_set1.csv").is_file()
    assert (outdir / "SRA" / "sra_file_set1.csv").is_file()
    assert (outdir / "BIOHUB" / "samples_template_v1.csv").is_file()


def test_empty_sample_df_skips_demultiplexing_samples(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    flattener = make_flattener(pd.DataFrame(), cro_orders=["ABC-123"])
    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")

    outdir_name = f"ABC-123_{TIMESTAMP}"
    outdir = tmp_path / outdir_name
    assert result == outdir_name
    assert (outdir / "MAIN.csv").is_file()
    assert (outdir / "GEO" / "geo_metadata.csv").is_file()
    assert not (outdir / "GEO" / "demultiplexing_samples.csv").exists()
    assert not (outdir / "BIOHUB" / "guide_template_set1_v1.csv").exists()


def test_guide_template_set_written_when_present(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    flattener = make_flattener(
        pd.DataFrame({"sample_alias": ["s1"]}),
        cro_orders=["ABC-123"],
    )
    flattener._biohub_guide_sets = lambda complete_data: [
        {
            "set_number": 1,
            "gm_id": "/genetic_modifications/aaa/",
            "file_info": {"@id": "/tabular_files/x/"},
            "library_ids": set(),
        }
    ]
    flattener.create_guide_metadata_dataframe = lambda file_info: pd.DataFrame({"id": ["g1"]})

    result = flattener.flatten_matrix_file_set(UUID, biohub_version="v1")

    outdir_name = f"ABC-123_{TIMESTAMP}"
    outdir = tmp_path / outdir_name
    assert result == outdir_name
    assert (outdir / "BIOHUB" / "samples_template_set1_v1.csv").is_file()
    assert (outdir / "BIOHUB" / "guide_template_set1_v1.csv").is_file()
    assert not (outdir / "BIOHUB" / "samples_template_v1.csv").exists()


def test_empty_guide_df_skips_guide_template(tmp_path, monkeypatch):
    _freeze_time(monkeypatch)
    monkeypatch.chdir(tmp_path)

    flattener = make_flattener(
        pd.DataFrame({"sample_alias": ["s1"]}),
        cro_orders=["ABC-123"],
    )
    flattener._biohub_guide_sets = lambda complete_data: [
        {
            "set_number": 1,
            "gm_id": "/genetic_modifications/aaa/",
            "file_info": {"@id": "/tabular_files/x/"},
            "library_ids": set(),
        }
    ]
    flattener.create_guide_metadata_dataframe = lambda file_info: pd.DataFrame()

    flattener.flatten_matrix_file_set(UUID, biohub_version="v1")

    outdir = tmp_path / f"ABC-123_{TIMESTAMP}"
    assert (outdir / "BIOHUB" / "samples_template_set1_v1.csv").is_file()
    assert not (outdir / "BIOHUB" / "guide_template_set1_v1.csv").exists()
