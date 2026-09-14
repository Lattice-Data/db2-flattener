from datetime import datetime
from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import Configs

UUID = "9b1c06e5-aaaa-bbbb-cccc-dddddddddddd"


def make_flattener(sample_df):
    main_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1"],
            "sample_alias": ["s1"],
            "droplet_based_libraries_aliases": [["alex-marson:lib1"]],
        }
    )
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: {"libraries": {}})
    f.create_dataframe = lambda complete_data: (main_df.copy(), sample_df.copy())
    f.create_biohub_dataframe = lambda main: pd.DataFrame({"sample_name": ["s1"]})
    f.create_geo_dataframe = lambda main: pd.DataFrame({"library_name": ["lib1"]})
    return f


def test_output_prefix_names_all_csvs(tmp_path):
    flattener = make_flattener(pd.DataFrame({"sample_alias": ["s1"]}))
    prefix = str(tmp_path / "myrun")

    result = flattener.flatten_matrix_file_set(UUID, output_prefix=prefix)

    assert result == f"{prefix}_MAIN.csv"
    for suffix in ("MAIN", "BIOHUB", "GEO", "SRA_BIOSAMPLE", "SAMPLES"):
        assert (tmp_path / f"myrun_{suffix}.csv").is_file()


def test_default_prefix_uses_uuid_and_timestamp(tmp_path, monkeypatch):
    class FrozenDateTime:
        @staticmethod
        def now():
            return datetime(2026, 8, 14, 9, 50, 0)

    monkeypatch.setattr("db2_flattener.flatten.flattener.datetime", FrozenDateTime)
    monkeypatch.chdir(tmp_path)

    flattener = make_flattener(pd.DataFrame({"sample_alias": ["s1"]}))
    result = flattener.flatten_matrix_file_set(UUID)

    prefix = f"MatrixFileSet_{UUID[:8]}_20260814_095000"
    assert result == f"{prefix}_MAIN.csv"
    for suffix in ("MAIN", "BIOHUB", "GEO", "SRA_BIOSAMPLE", "SAMPLES"):
        assert (tmp_path / f"{prefix}_{suffix}.csv").is_file()


def test_empty_sample_df_skips_samples_csv(tmp_path):
    flattener = make_flattener(pd.DataFrame())
    prefix = str(tmp_path / "myrun")

    result = flattener.flatten_matrix_file_set(UUID, output_prefix=prefix)

    assert result == f"{prefix}_MAIN.csv"
    assert (tmp_path / "myrun_MAIN.csv").is_file()
    assert (tmp_path / "myrun_BIOHUB.csv").is_file()
    assert (tmp_path / "myrun_GEO.csv").is_file()
    assert (tmp_path / "myrun_SRA_BIOSAMPLE.csv").is_file()
    assert not (tmp_path / "myrun_SAMPLES.csv").exists()
    assert not (tmp_path / "myrun_GUIDE_METADATA.csv").exists()


def test_guide_metadata_csv_written_when_present(tmp_path):
    flattener = make_flattener(pd.DataFrame({"sample_alias": ["s1"]}))
    flattener.create_guide_metadata_dataframe = lambda file_info: pd.DataFrame({"guide_id": ["g1"]})
    prefix = str(tmp_path / "myrun")

    result = flattener.flatten_matrix_file_set(UUID, output_prefix=prefix)

    assert result == f"{prefix}_MAIN.csv"
    assert (tmp_path / "myrun_GUIDE_METADATA.csv").is_file()


def test_empty_guide_df_skips_guide_metadata_csv(tmp_path):
    flattener = make_flattener(pd.DataFrame({"sample_alias": ["s1"]}))
    flattener.create_guide_metadata_dataframe = lambda file_info: pd.DataFrame()
    prefix = str(tmp_path / "myrun")

    flattener.flatten_matrix_file_set(UUID, output_prefix=prefix)

    assert not (tmp_path / "myrun_GUIDE_METADATA.csv").exists()
