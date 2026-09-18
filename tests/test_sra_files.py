"""
SRA file sheet: one row per CRO group × library, keyed on the library's
CRO_group_identifier under the name 'sample_name' and the cleaned library
alias under 'library_ID'. library_strategy is derived from feature_types.
Every library is kept, including non-GEX.
"""

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import PROP_MAP_SRA_FILE, Configs

LIBRARY_ID_COLUMN = "library_ID"
LIBRARY_STRATEGY_COLUMN = "library_strategy"


def make_flattener():
    """Empty configs on purpose: this sheet is built from main_df alone."""
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = None
    flattener.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    return flattener


def test_paired_libraries_each_get_a_row():
    """A group's GEX and CRISPR libraries share one sample_name, two library_IDs."""
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_A"],
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_A_CRI"],
            ],
            "droplet_based_libraries_@id": [
                "/droplet_based_libraries/gex/",
                "/droplet_based_libraries/cri/",
            ],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
            ],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df.columns) == ["sample_name", LIBRARY_ID_COLUMN, LIBRARY_STRATEGY_COLUMN]
    assert list(sra_df["sample_name"]) == ["LIB_A", "LIB_A"]
    assert list(sra_df[LIBRARY_ID_COLUMN]) == ["LIB_A_GEX", "LIB_A_CRI"]
    assert list(sra_df[LIBRARY_STRATEGY_COLUMN]) == ["RNA-Seq", "OTHER"]


def test_crispr_only_group_is_kept():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_B"],
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_B_CRI"],
            ],
            "droplet_based_libraries_@id": [
                "/droplet_based_libraries/gex/",
                "/droplet_based_libraries/cri/",
            ],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
            ],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["LIB_A", "LIB_B"]
    assert list(sra_df[LIBRARY_ID_COLUMN]) == ["LIB_A_GEX", "LIB_B_CRI"]
    assert list(sra_df[LIBRARY_STRATEGY_COLUMN]) == ["RNA-Seq", "OTHER"]


def test_plate_library_columns_are_used_when_droplet_is_absent():
    main_df = pd.DataFrame(
        {
            "plate_based_libraries_CRO_group_identifier": ["PLATE_1"],
            "plate_based_libraries_aliases": [["alex-marson:PLATE_1_GEX"]],
            "plate_based_libraries_@id": ["/plate_based_libraries/p1/"],
            "plate_based_libraries_feature_types": [["Gene Expression"]],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["PLATE_1"]
    assert list(sra_df[LIBRARY_ID_COLUMN]) == ["PLATE_1_GEX"]
    assert list(sra_df[LIBRARY_STRATEGY_COLUMN]) == ["RNA-Seq"]


def test_duplicate_rmf_rows_for_the_same_library_collapse_to_one():
    main_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1", "rmf2"],
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_A"],
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_A_GEX"],
            ],
            "droplet_based_libraries_@id": [
                "/droplet_based_libraries/gex/",
                "/droplet_based_libraries/gex/",
            ],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert len(sra_df) == 1
    assert sra_df.loc[0, "sample_name"] == "LIB_A"
    assert sra_df.loc[0, LIBRARY_ID_COLUMN] == "LIB_A_GEX"
    assert pd.isna(sra_df.loc[0, LIBRARY_STRATEGY_COLUMN])


def test_missing_library_group_column_returns_empty_frame(capsys):
    main_df = pd.DataFrame({"droplet_based_libraries_aliases": [["alex-marson:LIB_A"]]})

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert sra_df.empty
    assert "no library CRO group identifier column" in capsys.readouterr().out


def test_rows_with_no_library_group_are_dropped(capsys):
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", None],
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_B_GEX"],
            ],
            "droplet_based_libraries_@id": [
                "/droplet_based_libraries/a/",
                "/droplet_based_libraries/b/",
            ],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["LIB_A"]
    assert list(sra_df[LIBRARY_ID_COLUMN]) == ["LIB_A_GEX"]
    assert "dropping 1 of 2 MAIN row(s)" in capsys.readouterr().out


def test_library_id_strips_the_lab_prefix():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["some-other-lab:LIB_A_GEX"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df[LIBRARY_ID_COLUMN]) == ["LIB_A_GEX"]


def test_missing_aliases_column_omits_library_id(capsys):
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["LIB_A"]
    assert LIBRARY_ID_COLUMN not in sra_df.columns
    assert LIBRARY_STRATEGY_COLUMN in sra_df.columns
    assert "no library aliases column" in capsys.readouterr().out


def test_library_strategy_maps_feature_types():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": [
                "LIB_GEX",
                "LIB_CRI",
                "LIB_MUX",
                "LIB_ATAC",
            ],
            "droplet_based_libraries_aliases": [
                ["lab:LIB_GEX"],
                ["lab:LIB_CRI"],
                ["lab:LIB_MUX"],
                ["lab:LIB_ATAC"],
            ],
            "droplet_based_libraries_@id": [
                "/droplet_based_libraries/gex/",
                "/droplet_based_libraries/cri/",
                "/droplet_based_libraries/mux/",
                "/droplet_based_libraries/atac/",
            ],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
                ["Multiplexing Capture"],
                ["ATAC"],
            ],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df[LIBRARY_STRATEGY_COLUMN]) == ["RNA-Seq", "OTHER", "OTHER", "ATAC-seq"]
    assert "droplet_based_libraries_feature_types" not in sra_df.columns


def test_library_strategy_maps_string_feature_types():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["lab:LIB_A_GEX"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
            "droplet_based_libraries_feature_types": ["Gene Expression"],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df[LIBRARY_STRATEGY_COLUMN]) == ["RNA-Seq"]


def test_library_strategy_is_blank_when_feature_types_are_missing_or_unmapped():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_B"],
            "droplet_based_libraries_aliases": [["lab:LIB_A"], ["lab:LIB_B"]],
            "droplet_based_libraries_@id": [
                "/droplet_based_libraries/a/",
                "/droplet_based_libraries/b/",
            ],
            "droplet_based_libraries_feature_types": [None, ["Antibody Capture"]],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert list(sra_df.columns) == ["sample_name", LIBRARY_ID_COLUMN, LIBRARY_STRATEGY_COLUMN]
    assert sra_df[LIBRARY_STRATEGY_COLUMN].isna().all()


def test_library_strategy_is_blank_when_feature_types_column_is_absent():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["lab:LIB_A_GEX"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
        }
    )

    sra_df = make_flattener().create_sra_files_dataframe(main_df)

    assert LIBRARY_STRATEGY_COLUMN in sra_df.columns
    assert pd.isna(sra_df.loc[0, LIBRARY_STRATEGY_COLUMN])


def test_prop_map_sends_both_library_types_to_sample_name():
    droplet_key = "droplet_based_libraries_CRO_group_identifier"
    assert PROP_MAP_SRA_FILE[droplet_key] == "sample_name"
    assert PROP_MAP_SRA_FILE["plate_based_libraries_CRO_group_identifier"] == "sample_name"
    assert "sample_alias" not in PROP_MAP_SRA_FILE
    assert "droplet_based_libraries_aliases" not in PROP_MAP_SRA_FILE
    assert "droplet_based_libraries_feature_types" not in PROP_MAP_SRA_FILE
    assert not PROP_MAP_SRA_FILE[droplet_key].startswith("*")
    assert LIBRARY_ID_COLUMN not in PROP_MAP_SRA_FILE.values()
    assert LIBRARY_STRATEGY_COLUMN not in PROP_MAP_SRA_FILE.values()


def test_empty_main_df_returns_empty_frame():
    main_df = pd.DataFrame(columns=["droplet_based_libraries_CRO_group_identifier"])

    assert make_flattener().create_sra_files_dataframe(main_df).empty
