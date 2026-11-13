"""Refuse to flatten multiome and ATAC libraries."""

from types import SimpleNamespace

import pandas as pd
import pytest

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import Configs

TERM_ID_COL = "droplet_based_libraries_library_construction_technology_term_id"


def make_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    f.gatherer = None
    return f


def library_row(group, feature_types="Gene Expression", term_id="EFO:0009922", **overrides):
    row = {
        "droplet_based_libraries_@id": f"/droplet_based_libraries/{group}-{feature_types}/",
        "droplet_based_libraries_CRO_group_identifier": group,
        "droplet_based_libraries_feature_types": feature_types,
        TERM_ID_COL: term_id,
    }
    row.update(overrides)
    return row


@pytest.mark.parametrize("term_id", ["EFO:0030059", "EFO:0920134", "EFO:0920135"])
def test_refuses_multiome_library(term_id, capsys):
    main_df = pd.DataFrame([library_row("libA", term_id=term_id), library_row("libB")])

    with pytest.raises(ValueError, match="Refusing to flatten: multiome libraries"):
        DB2Flattener._refuse_multiome_libraries(main_df)

    out = capsys.readouterr().out
    assert "====WARNING==== 1 library group(s) use a multiome assay:" in out
    assert f"  libA ({term_id})" in out
    assert "libB" not in out


def test_refuses_multiome_term_id_inside_a_list():
    main_df = pd.DataFrame([library_row("libA", term_id=["EFO:0009922", "EFO:0030059"])])

    with pytest.raises(ValueError, match="multiome"):
        DB2Flattener._refuse_multiome_libraries(main_df)


def test_refuses_plate_multiome_library():
    main_df = pd.DataFrame(
        [
            {
                "plate_based_libraries_CRO_group_identifier": "plateA",
                "plate_based_libraries_library_construction_technology_term_id": "EFO:0920135",
            }
        ]
    )

    with pytest.raises(ValueError, match="multiome"):
        DB2Flattener._refuse_multiome_libraries(main_df)


def test_allows_non_multiome_and_missing_term_ids():
    main_df = pd.DataFrame([library_row("libA"), library_row("libB", term_id=None)])
    DB2Flattener._refuse_multiome_libraries(main_df)
    DB2Flattener._refuse_multiome_libraries(pd.DataFrame({"raw_matrix_file_alias": ["a.h5"]}))


def test_refuses_atac_library_even_alongside_gex(capsys):
    main_df = pd.DataFrame(
        [
            library_row("libA"),
            library_row("libA", feature_types="ATAC"),
            library_row("libB"),
        ]
    )

    with pytest.raises(ValueError, match="Refusing to flatten: ATAC libraries"):
        DB2Flattener._refuse_atac_libraries(main_df)

    out = capsys.readouterr().out
    assert "====WARNING==== 1 library group(s) have an ATAC library:" in out
    assert "  libA" in out
    assert "libB" not in out


def test_refuses_atac_in_a_feature_type_list():
    main_df = pd.DataFrame([library_row("libA", feature_types=["Gene Expression", "ATAC"])])

    with pytest.raises(ValueError, match="ATAC"):
        DB2Flattener._refuse_atac_libraries(main_df)


def test_refuses_plate_atac_library():
    main_df = pd.DataFrame(
        [
            {
                "plate_based_libraries_CRO_group_identifier": "plateA",
                "plate_based_libraries_feature_types": "ATAC",
            }
        ]
    )

    with pytest.raises(ValueError, match="ATAC"):
        DB2Flattener._refuse_atac_libraries(main_df)


def test_refuses_atac_library_without_group_id(capsys):
    main_df = pd.DataFrame([library_row(None, feature_types="ATAC")])

    with pytest.raises(ValueError, match="ATAC"):
        DB2Flattener._refuse_atac_libraries(main_df)

    assert "  (no GroupID)" in capsys.readouterr().out


def test_allows_libraries_without_atac():
    main_df = pd.DataFrame(
        [
            library_row("libA"),
            library_row("libA", feature_types="CRISPR Guide Capture"),
            library_row("libB", feature_types=None),
        ]
    )
    DB2Flattener._refuse_atac_libraries(main_df)
    DB2Flattener._refuse_atac_libraries(pd.DataFrame({"raw_matrix_file_alias": ["a.h5"]}))


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([library_row("libA", term_id="EFO:0030059")], "multiome libraries"),
        ([library_row("libA", feature_types="ATAC")], "ATAC libraries"),
    ],
)
def test_flatten_refuses_before_writing_outputs(rows, message, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    complete_data = {"resolved_objects": {"GeneticModification": {}}, "libraries": {}}
    f = make_flattener()
    f.gatherer = SimpleNamespace(gather_complete_library_data=lambda uuid: complete_data)
    f.create_dataframe = lambda data: (pd.DataFrame(rows), None)

    with pytest.raises(ValueError, match=f"Refusing to flatten: .*{message}"):
        f.flatten_matrix_file_set("abc")

    assert "====WARNING====" in capsys.readouterr().out
    assert not any(tmp_path.iterdir())
