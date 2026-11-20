"""A donor with no ethnicities pooled with a known one must not read as one ethnicity.

HumanDonor.ethnicities can be empty. Such a donor adds nothing to the merged
terms, so create_dataframe records it in human_donors_without_ethnicities and
GEO, SRA biosample and SAMPLES add 'not provided' for it. BIOHUB is untouched:
its samples are pre-pooled.
"""

import pandas as pd
import pytest
import test_controlled_terms
from test_controlled_terms import EUROPEAN, _main_df, _tissue
from test_controlled_terms import make_flattener as make_terms_flattener
from test_geo_dataframe import gex_row
from test_geo_dataframe import make_flattener as make_geo_flattener
from test_samples_dataframe import make_flattener as make_samples_flattener
from test_sra_biosample import ETHNICITY_COLUMN, donor_main_df
from test_sra_biosample import make_flattener as make_sra_flattener

from db2_flattener.schema.constants import Configs

D1 = "/human_donors/d1/"
D2 = "/human_donors/d2/"
BOTH_DONORS = f"{D1}; {D2}"
UNKNOWN_COL = "human_donors_without_ethnicities"

KNOWN = {"@id": D1, "ethnicities": [EUROPEAN["@id"]]}
ALSO_KNOWN = {"@id": D2, "ethnicities": [EUROPEAN["@id"]]}
UNKNOWN = {"@id": D2, "ethnicities": []}


# --- create_dataframe records the unknown donors ---


def main_row(donors):
    tissue = _tissue(donors=[d["@id"] for d in donors])
    return _main_df([tissue], extra_lib_data={"human_donors": donors}).iloc[0]


def library_cell(row):
    """GEO/SRA's per-library cell for a one-row MAIN."""
    main_df = pd.DataFrame([row])
    f = make_terms_flattener()
    key = pd.Series(["LIB"], index=main_df.index)
    return f._ethnicity_by_library(main_df, key, f._tissue_row_mask(main_df))["LIB"]


def test_main_records_the_donor_with_no_ethnicities():
    row = main_row([KNOWN, UNKNOWN])

    assert row[UNKNOWN_COL] == D2
    assert row["human_donors_ethnicities_term_name"] == "European"


def test_main_leaves_the_column_empty_when_every_donor_is_known():
    row = main_row([KNOWN, ALSO_KNOWN])

    assert pd.isna(row[UNKNOWN_COL])


def test_known_plus_unknown_donor_in_one_sample_is_pooled_with_the_default():
    assert library_cell(main_row([KNOWN, UNKNOWN])) == "pooled: European, not provided"


def test_two_known_donors_with_the_same_ethnicity_stay_bare():
    assert library_cell(main_row([KNOWN, ALSO_KNOWN])) == "European"


def test_nothing_is_recorded_when_ethnicities_is_not_a_fetched_field(monkeypatch):
    """Without the field every donor would look unknown, so nothing is recorded."""
    object_config = dict(test_controlled_terms.CONFIGS.OBJECT_CONFIG)
    object_config["human_donors"] = {
        "api_type": "HumanDonor",
        "fields": ["@id"],
        "references": {},
    }
    monkeypatch.setattr(
        test_controlled_terms,
        "CONFIGS",
        Configs(FIELD_TYPES=test_controlled_terms.CONFIGS.FIELD_TYPES, OBJECT_CONFIG=object_config),
    )

    row = main_row([KNOWN, UNKNOWN])

    assert UNKNOWN_COL not in row.index


# --- GEO ---


def geo_cell(*rows):
    """rows: (terms, donor ids, donors without ethnicities)"""
    main_df = pd.DataFrame(
        [
            gex_row(
                **{
                    "tissues_@id": f"/tissues/t{i}/",
                    "human_donors_ethnicities_term_name": terms,
                    "human_donors_@id": ids,
                    UNKNOWN_COL: unknown,
                }
            )
            for i, (terms, ids, unknown) in enumerate(rows, start=1)
        ]
    ).dropna(axis=1, how="all")
    return list(make_geo_flattener().create_geo_dataframe(main_df)["donor_ethnicity"])


def test_geo_known_plus_unknown_donor_adds_the_default():
    assert geo_cell((["European"], BOTH_DONORS, D2)) == ["pooled: European, not provided"]


def test_geo_without_an_unknown_donor_nothing_is_added():
    assert geo_cell((["European"], BOTH_DONORS, None)) == ["European"]


def test_geo_the_default_is_added_once_per_library():
    cells = geo_cell(
        (["European"], D1, None),
        (None, D2, D2),
        (["European"], BOTH_DONORS, D2),
    )

    assert cells == ["pooled: European, not provided"]


def test_geo_several_terms_and_an_unknown_donor():
    cells = geo_cell((["Asian", "European"], BOTH_DONORS, D2))

    assert cells == ["pooled: Asian, European, not provided"]


# --- SRA biosample ---


def sra_cell(**columns):
    main_df = donor_main_df(**columns)
    return make_sra_flattener().create_sra_biosample_dataframe(main_df).loc[0, ETHNICITY_COLUMN]


def test_sra_known_plus_unknown_donor_adds_the_default():
    cell = sra_cell(
        **{
            "human_donors_ethnicities_term_name": [["European"]] * 4,
            "human_donors_@id": [BOTH_DONORS] * 4,
            UNKNOWN_COL: [D2] * 4,
        }
    )

    assert cell == "pooled: European, not provided"


def test_sra_without_an_unknown_donor_nothing_is_added():
    cell = sra_cell(
        **{
            "human_donors_ethnicities_term_name": [["European"]] * 4,
            "human_donors_@id": [BOTH_DONORS] * 4,
        }
    )

    assert "not provided" not in cell


def test_sra_the_default_is_added_once_when_a_row_also_lacks_terms():
    cell = sra_cell(
        **{
            "human_donors_ethnicities_term_name": [["European"], None, ["European"], None],
            "human_donors_@id": [D1, D2, D1, D2],
            UNKNOWN_COL: [None, D2, None, D2],
        }
    )

    assert cell == "pooled: European, not provided"


# --- SAMPLES ---


def samples_result(rows):
    """rows: (terms, donor ids, donors without ethnicities)"""
    sample_df = pd.DataFrame(
        {
            "sample_alias": [f"s{i}" for i in range(len(rows))],
            "tissues_@id": [f"/tissues/t{i}/" for i in range(len(rows))],
            "human_donors_ethnicities_term_name": [terms for terms, _, _ in rows],
            "human_donors_@id": [ids for _, ids, _ in rows],
            UNKNOWN_COL: [unknown for _, _, unknown in rows],
        }
    )
    return make_samples_flattener().create_samples_dataframe(sample_df)


def test_samples_known_plus_unknown_donor_adds_the_default():
    result = samples_result([(["European"], BOTH_DONORS, D2)])

    assert list(result["donor_ethnicity"]) == ["pooled: European, not provided"]
    assert UNKNOWN_COL not in result.columns
    assert "human_donors_@id" not in result.columns


def test_samples_per_row():
    result = samples_result(
        [
            (["European"], BOTH_DONORS, D2),
            (["European"], D1, None),
            (None, D2, D2),
        ]
    )

    assert list(result["donor_ethnicity"]) == [
        "pooled: European, not provided",
        "European",
        "not provided",
    ]


@pytest.mark.parametrize("unknown", [None, ""])
def test_samples_no_unknown_donor_changes_nothing(unknown):
    result = samples_result([(["Asian", "European"], BOTH_DONORS, unknown)])

    assert list(result["donor_ethnicity"]) == ["pooled: Asian, European"]
