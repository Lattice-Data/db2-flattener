"""
SRA/BioSample sheet: one row per library, keyed on the library alias under the
name 'sample_name' because the submitted 'sample' is the sequencing library.
"""

import pandas as pd
import pytest

from db2_flattener.flatten.flattener import (
    SAMPLE_URL_PREFIXES,
    DB2Flattener,
    age_with_units,
)
from db2_flattener.schema.constants import (
    PROP_MAP_SRA_BIOSAMPLE,
    TISSUE_TYPE_MAP,
    Configs,
)
from db2_flattener.utils import age_from_developmental_stage

# Derived columns, named in the flattener rather than in PROP_MAP_SRA_BIOSAMPLE.
BARCODE_COLUMN = "sample_name: sample_probe_barcode"
ISOLATE_COLUMN = "*isolate"
AGE_COLUMN = "*age"
SEX_COLUMN = "*sex"
PROVIDER_COLUMN = "*biomaterial_provider"
DATE_COLUMN = "*collection_date"
GEO_COLUMN = "*geo_loc_name"
TISSUE_COLUMN = "*tissue"
# Optional columns: present only when the data is
AGE_LOWER_COLUMN = "age_lower_bound"
AGE_UPPER_COLUMN = "age_upper_bound"
ETHNICITY_COLUMN = "ethnicity"
SUSPENSION_COLUMN = "suspension_type"
PERTURBATION_COLUMN = "experimental_perturbation"
FACTORS_COLUMN = "experimental_perturbation_factors"
# Read off the map so a rename there, e.g. dropping the SRA '*' required marker,
# does not have to be chased through every assertion below.
ORGANISM_COLUMN = PROP_MAP_SRA_BIOSAMPLE["human_donors_taxa"]

LAB = {"@id": "/labs/alex-marson/", "title": "Alex Marson, UCSF"}
OTHER_LAB = {"@id": "/labs/other/", "title": "Other Lab"}
SOURCE = {"@id": "/sources/abcam/", "title": "Abcam"}


def make_flattener():
    """Empty configs on purpose: this sheet is built from main_df alone."""
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = None
    flattener.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    return flattener


def sample(alias, barcodes, lab="alex-marson"):
    return {
        "@id": f"/tissues/{alias}/",
        "aliases": [f"{lab}:{alias}"],
        "multiplexing_barcodes": barcodes,
    }


# Deliberately not in alphabetical order: the builder is expected to sort.
FOUR_SAMPLES = [
    sample("TregD5_Rest", ["BC001+CR001", "BC002+CR002"]),
    sample("TregD6_Rest", ["BC009+CR009", "BC010+CR010"]),
    sample("TregD5_Stim8hr", ["BC005+CR005", "BC006+CR006"]),
    sample("TregD6_Stim8hr", ["BC013+CR013", "BC014+CR014"]),
]

FOUR_SAMPLES_MAP = (
    "TregD5_Rest : BC001+CR001|BC002+CR002, "
    "TregD5_Stim8hr : BC005+CR005|BC006+CR006, "
    "TregD6_Rest : BC009+CR009|BC010+CR010, "
    "TregD6_Stim8hr : BC013+CR013|BC014+CR014"
)


# _sample_probe_barcode_map


def test_barcode_map_sorts_by_alias_and_strips_lab_prefix():
    assert make_flattener()._sample_probe_barcode_map(FOUR_SAMPLES) == FOUR_SAMPLES_MAP


def test_barcode_map_is_independent_of_source_order():
    flattener = make_flattener()
    forward = flattener._sample_probe_barcode_map(FOUR_SAMPLES)
    reversed_ = flattener._sample_probe_barcode_map(list(reversed(FOUR_SAMPLES)))
    assert forward == reversed_


def test_barcode_map_single_sample():
    samples = [sample("s1", ["BC001+CR001", "BC002+CR002"])]
    assert make_flattener()._sample_probe_barcode_map(samples) == "s1 : BC001+CR001|BC002+CR002"


def test_barcode_map_keeps_sample_with_no_barcodes_when_another_has_some():
    samples = [sample("s2", []), sample("s1", ["BC001+CR001"])]
    assert make_flattener()._sample_probe_barcode_map(samples) == "s1 : BC001+CR001, s2 : "


@pytest.mark.parametrize(
    "library_samples",
    [
        pytest.param([], id="empty-list"),
        pytest.param(None, id="none"),
        pytest.param(float("nan"), id="nan"),
        pytest.param([sample("s1", []), sample("s2", None)], id="no-sample-has-barcodes"),
        pytest.param(["/tissues/s1/", "/tissues/s2/"], id="bare-id-strings"),
        pytest.param([{"multiplexing_barcodes": ["BC001+CR001"]}], id="no-alias"),
    ],
)
def test_barcode_map_returns_none(library_samples):
    assert make_flattener()._sample_probe_barcode_map(library_samples) is None


def test_barcode_map_accepts_a_bare_barcode_string():
    samples = [sample("s1", "BC001+CR001")]
    assert make_flattener()._sample_probe_barcode_map(samples) == "s1 : BC001+CR001"


# create_sra_biosample_dataframe


def droplet_main_df():
    """Four samples pooled into one droplet library, one MAIN row per sample."""
    aliases = [s["aliases"][0].split(":", 1)[1] for s in FOUR_SAMPLES]
    return pd.DataFrame(
        {
            "sample_alias": aliases,
            "droplet_based_libraries_aliases": [["alex-marson:TregR3_L13_L05_GEX"]] * 4,
            "droplet_based_libraries_samples": [FOUR_SAMPLES] * 4,
            "human_donors_taxa": ["Homo sapiens"] * 4,
        }
    )


def test_droplet_run_is_one_row_per_library(capsys):
    sra_df = make_flattener().create_sra_biosample_dataframe(droplet_main_df())

    # No donor id or lab column here, so only the always-on columns are derived
    assert list(sra_df.columns) == [
        "sample_name",
        ORGANISM_COLUMN,
        BARCODE_COLUMN,
        DATE_COLUMN,
        GEO_COLUMN,
    ]
    assert "no donor id column" in capsys.readouterr().out
    assert len(sra_df) == 1
    assert sra_df.loc[0, "sample_name"] == "TregR3_L13_L05_GEX"
    assert sra_df.loc[0, ORGANISM_COLUMN] == "Homo sapiens"
    assert sra_df.loc[0, DATE_COLUMN] == "not provided"
    assert sra_df.loc[0, BARCODE_COLUMN] == FOUR_SAMPLES_MAP


def test_plate_library_columns_are_used_when_droplet_is_absent():
    main_df = pd.DataFrame(
        {
            "plate_based_libraries_aliases": [["alex-marson:PLATE_1"]],
            "plate_based_libraries_samples": [[sample("s1", ["BC001+CR001"])]],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["PLATE_1"]
    assert list(sra_df[BARCODE_COLUMN]) == ["s1 : BC001+CR001"]


def test_paired_libraries_get_their_own_rows_sharing_the_same_map():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_A_CRI"]] * 2,
            "droplet_based_libraries_samples": [FOUR_SAMPLES] * 4,
            "human_donors_taxa": ["Homo sapiens"] * 4,
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sorted(sra_df.index) == ["LIB_A_CRI", "LIB_A_GEX"]
    assert sra_df.loc["LIB_A_GEX", BARCODE_COLUMN] == FOUR_SAMPLES_MAP
    assert sra_df.loc["LIB_A_CRI", BARCODE_COLUMN] == FOUR_SAMPLES_MAP


def test_library_rows_disagreeing_on_the_map_collapse_to_a_list():
    other_samples = [sample("s1", ["BC099+CR099"])]
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2,
            "droplet_based_libraries_samples": [FOUR_SAMPLES, other_samples],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert len(sra_df) == 1
    assert sra_df.loc[0, BARCODE_COLUMN] == [FOUR_SAMPLES_MAP, "s1 : BC099+CR099"]


def test_rows_with_no_library_alias_are_dropped(capsys):
    main_df = droplet_main_df()
    main_df.loc[1, "droplet_based_libraries_aliases"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["TregR3_L13_L05_GEX"]
    assert "dropping 1 of 4 MAIN row(s) with no library alias" in capsys.readouterr().out


def test_sample_alias_is_not_the_grouping_key():
    """The library alias keys the sheet, so nulling sample_alias changes nothing."""
    expected = make_flattener().create_sra_biosample_dataframe(droplet_main_df())

    main_df = droplet_main_df()
    main_df["sample_alias"] = None

    actual = make_flattener().create_sra_biosample_dataframe(main_df)

    pd.testing.assert_frame_equal(actual, expected)


def test_missing_library_aliases_column_returns_empty_frame(capsys):
    main_df = pd.DataFrame({"sample_alias": ["s1"]})

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.empty
    assert "no library aliases column" in capsys.readouterr().out


def test_alias_only_main_df_collapses_without_aggregating():
    """groupby().agg({}) raises, so a frame of nothing but the key dedupes instead."""
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_B_GEX"],
            ]
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df.columns) == ["sample_name", DATE_COLUMN, GEO_COLUMN]
    assert list(sra_df["sample_name"]) == ["LIB_A_GEX", "LIB_B_GEX"]
    assert list(sra_df[DATE_COLUMN]) == ["not provided", "not provided"]
    assert list(sra_df[GEO_COLUMN]) == ["not provided", "not provided"]


def test_library_without_barcodes_leaves_the_column_empty():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]],
            "droplet_based_libraries_samples": [[sample("s1", [])]],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert pd.isna(sra_df.loc[0, BARCODE_COLUMN])


def test_empty_main_df_returns_empty_frame():
    main_df = pd.DataFrame({"droplet_based_libraries_aliases": pd.Series(dtype=object)})

    assert make_flattener().create_sra_biosample_dataframe(main_df).empty


# _clean_alias_cell


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(["alex-marson:LIB_A_GEX"], "LIB_A_GEX", id="list"),
        pytest.param("alex-marson:LIB_A_GEX", "LIB_A_GEX", id="bare-string"),
        pytest.param("LIB_A_GEX", "LIB_A_GEX", id="string-without-prefix"),
        pytest.param([], None, id="empty-list"),
        pytest.param("", "", id="empty-string"),
        pytest.param(None, None, id="none"),
    ],
)
def test_clean_alias_cell(value, expected):
    assert make_flattener()._clean_alias_cell(value) == expected


# isolate and age


def test_age_from_developmental_stage_parses_years():
    assert age_from_developmental_stage("29-year-old stage") == "29 years"


@pytest.mark.parametrize(
    ("term_name", "expected"),
    [
        pytest.param("29-year-old stage", "29 years", id="years"),
        pytest.param("1-year-old stage", "1 year", id="singular-year"),
        pytest.param("6-month-old stage", "6 months", id="months"),
        pytest.param("1-week-old stage", "1 week", id="singular-week"),
        pytest.param(
            ["adult stage", "42-year-old stage"], "42 years", id="numeric-wins-over-qualitative"
        ),
        pytest.param("adult stage", "adult", id="qualitative"),
        pytest.param("newborn stage", "newborn", id="newborn"),
        pytest.param("adult", "adult", id="no-stage-suffix"),
        pytest.param(
            "10th week post-fertilization human stage",
            "10th week post-fertilization",
            id="qualitative-containing-a-number",
        ),
        pytest.param("adult human stage", "adult", id="human-stage-suffix"),
        pytest.param("mouse adult stage", "mouse adult", id="only-human-is-trimmed"),
        pytest.param(["adult stage", "newborn stage"], "adult", id="first-qualitative-in-list"),
        pytest.param(None, None, id="none"),
        pytest.param(float("nan"), None, id="nan"),
        pytest.param([], None, id="empty-list"),
        pytest.param("   ", None, id="blank"),
    ],
)
def test_age_from_developmental_stage(term_name, expected):
    assert age_from_developmental_stage(term_name) == expected


def donor_main_df(**columns):
    """
    One library, four samples, two donors. Mirrors the real Treg run: 889023040 is
    the 32-year-old female (D1, the lower id), 889081306 the 29-year-old male (D2).
    """
    base = {
        "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 4,
        "droplet_based_libraries_samples": [FOUR_SAMPLES] * 4,
        "human_donors_cxg_donor_id": [889081306, 889023040, 889081306, 889023040],
        "human_donors_sex": ["male", "female", "male", "female"],
        "tissues_developmental_stages_term_name": [
            "29-year-old stage",
            "32-year-old stage",
            "29-year-old stage",
            "32-year-old stage",
        ],
    }
    base.update(columns)
    return pd.DataFrame(base)


def test_two_donors_are_pooled_and_numbered_by_donor_id():
    sra_df = make_flattener().create_sra_biosample_dataframe(donor_main_df())

    # D1 is the lower id, 889023040, whose samples are the 32-year-old ones.
    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"
    assert sra_df.loc[0, AGE_COLUMN] == "pooled: D1 - 32 years, D2 - 29 years"


def test_single_donor_has_no_pooled_prefix():
    main_df = donor_main_df()
    main_df["human_donors_cxg_donor_id"] = 889023040
    main_df["tissues_developmental_stages_term_name"] = "32-year-old stage"

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "889023040"
    assert sra_df.loc[0, AGE_COLUMN] == "32 years"


def test_numeric_donor_ids_sort_numerically_not_lexically():
    main_df = donor_main_df()
    main_df["human_donors_cxg_donor_id"] = [9, 10, 9, 10]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # Lexical order would put 10 before 9.
    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 9, D2 - 10"


def test_non_numeric_donor_ids_sort_lexically():
    main_df = donor_main_df()
    main_df["human_donors_cxg_donor_id"] = ["CE0010866", "CE0008162", "CE0010866", "CE0008162"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - CE0008162, D2 - CE0010866"


def test_donor_id_is_not_rendered_as_a_float():
    """A null in the column makes pandas store the ids as float64."""
    main_df = donor_main_df()
    main_df.loc[3, "human_donors_cxg_donor_id"] = None
    assert main_df["human_donors_cxg_donor_id"].dtype == "float64"

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"


def test_qualitative_stage_lands_in_age_alongside_a_numeric_one():
    main_df = donor_main_df()
    # The 889081306 rows, D2, have no numeric stage.
    main_df.loc[[0, 2], "tissues_developmental_stages_term_name"] = "adult stage"

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_COLUMN] == "pooled: D1 - 32 years, D2 - adult"


def test_donor_without_an_age_keeps_its_label_in_the_other_column():
    main_df = donor_main_df()
    # Blank the 889081306 rows, which are D2.
    main_df.loc[[0, 2], "tissues_developmental_stages_term_name"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"
    # D2 drops out of age, but D1 keeps its number so the labels still agree.
    assert sra_df.loc[0, AGE_COLUMN] == "pooled: D1 - 32 years"


def test_no_stage_at_all_leaves_the_age_column_empty():
    main_df = donor_main_df()
    main_df["tissues_developmental_stages_term_name"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert pd.isna(sra_df.loc[0, AGE_COLUMN])
    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"


def test_missing_stage_column_still_gives_isolate():
    main_df = donor_main_df().drop(columns=["tissues_developmental_stages_term_name"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"
    assert pd.isna(sra_df.loc[0, AGE_COLUMN])


def test_no_donor_column_omits_all_donor_columns(capsys):
    main_df = donor_main_df().drop(columns=["human_donors_cxg_donor_id"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert ISOLATE_COLUMN not in sra_df.columns
    assert AGE_COLUMN not in sra_df.columns
    assert SEX_COLUMN not in sra_df.columns
    assert "no donor id column" in capsys.readouterr().out


# *sex


@pytest.mark.parametrize(
    ("sexes", "expected"),
    [
        pytest.param(["male"], "male", id="male-only"),
        pytest.param(["female"], "female", id="female-only"),
        pytest.param(["female", "male"], "pooled male and female", id="male-listed-first"),
        pytest.param(["male", "female"], "pooled male and female", id="order-independent"),
        pytest.param(["male", "male", "female"], "pooled male and female", id="deduplicated"),
        pytest.param(["unknown", "male"], "pooled male and unknown", id="male-then-unknown"),
        pytest.param(["unknown", "female"], "pooled female and unknown", id="female-then-unknown"),
        pytest.param(
            ["unknown", "female", "male"],
            "pooled male, female and unknown",
            id="three-values",
        ),
        pytest.param(["unknown"], "unknown", id="unknown-alone-passes-through"),
        pytest.param(["  male  "], "male", id="whitespace-stripped"),
        pytest.param([None, "male"], "male", id="nulls-ignored"),
        pytest.param([], None, id="empty"),
        pytest.param([None, float("nan"), ""], None, id="all-empty"),
    ],
)
def test_format_pooled_sex(sexes, expected):
    assert make_flattener()._format_pooled_sex(sexes) == expected


def test_mixed_sex_pool_reads_male_first():
    sra_df = make_flattener().create_sra_biosample_dataframe(donor_main_df())

    assert sra_df.loc[0, SEX_COLUMN] == "pooled male and female"


def test_single_sex_library_is_the_bare_value():
    main_df = donor_main_df()
    main_df["human_donors_sex"] = "female"

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "female"


def test_sex_is_summarised_over_donors_not_rows():
    main_df = donor_main_df()
    main_df["human_donors_cxg_donor_id"] = 889023040
    main_df["human_donors_sex"] = ["female", "female", "female", "female"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "female"
    assert sra_df.loc[0, ISOLATE_COLUMN] == "889023040"


def test_missing_sex_column_leaves_the_sex_column_empty():
    main_df = donor_main_df().drop(columns=["human_donors_sex"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert pd.isna(sra_df.loc[0, SEX_COLUMN])
    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"


def test_non_human_sex_column_is_used_when_human_is_absent():
    main_df = donor_main_df().rename(columns={"human_donors_sex": "non_human_donors_sex"})

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "pooled male and female"


def test_sex_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "human_donors_cxg_donor_id": [11, 22, 33],
            "human_donors_sex": ["male", "female", "male"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", SEX_COLUMN] == "pooled male and female"
    assert sra_df.loc["LIB_B_GEX", SEX_COLUMN] == "male"


def test_non_human_donor_column_is_used_when_human_is_absent():
    main_df = donor_main_df().rename(
        columns={"human_donors_cxg_donor_id": "non_human_donors_cxg_donor_id"}
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"


def test_donor_cells_are_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "human_donors_cxg_donor_id": [11, 22, 33],
            "tissues_developmental_stages_term_name": [
                "40-year-old stage",
                "50-year-old stage",
                "60-year-old stage",
            ],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", ISOLATE_COLUMN] == "pooled: D1 - 11, D2 - 22"
    assert sra_df.loc["LIB_A_GEX", AGE_COLUMN] == "pooled: D1 - 40 years, D2 - 50 years"
    # A single-donor library gets the bare form even when a sibling is pooled.
    assert sra_df.loc["LIB_B_GEX", ISOLATE_COLUMN] == "33"
    assert sra_df.loc["LIB_B_GEX", AGE_COLUMN] == "60 years"


# *tissue


def tissue_main_df(prefix, term="blood", **columns):
    """One library, two samples of the given sample type."""
    base = {
        "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2,
        f"{prefix}_@id": [f"/{prefix}/s1/", f"/{prefix}/s2/"],
    }
    if term is not None:
        base[f"{prefix}_sample_terms_term_name"] = [term, term]
    base.update(columns)
    return pd.DataFrame(base)


@pytest.mark.parametrize("prefix", ["cell_lines", "primary_cell_cultures"])
def test_tissueless_sample_types_report_not_available(prefix):
    sra_df = make_flattener().create_sra_biosample_dataframe(tissue_main_df(prefix, "HeLa"))

    assert sra_df.loc[0, TISSUE_COLUMN] == "not available"


@pytest.mark.parametrize("prefix", ["tissues", "organoids"])
def test_tissue_and_organoid_report_their_sample_term(prefix):
    sra_df = make_flattener().create_sra_biosample_dataframe(tissue_main_df(prefix))

    assert sra_df.loc[0, TISSUE_COLUMN] == "blood"


def test_tissue_flattens_a_multi_term_sample():
    """sample_terms is an array, so two terms give two entries not a stringified list."""
    main_df = tissue_main_df("tissues", term=None)
    main_df["tissues_sample_terms_term_name"] = [["blood", "lung"], ["blood", "lung"]]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, TISSUE_COLUMN] == "blood; lung"


def test_tissue_mixes_a_real_term_with_not_available():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2,
            "tissues_@id": ["/tissues/s1/", None],
            "tissues_sample_terms_term_name": ["blood", None],
            "cell_lines_@id": [None, "/cell_lines/s2/"],
            "cell_lines_sample_terms_term_name": [None, "HeLa"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, TISSUE_COLUMN] == "blood; not available"


def test_tissue_with_no_sample_term_stays_empty():
    """Should not happen - the schema requires it - but must not invent a value."""
    main_df = tissue_main_df("tissues", term=None)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert TISSUE_COLUMN not in sra_df.columns


def test_tissue_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]]
            + [["alex-marson:LIB_B_GEX"]],
            "tissues_@id": ["/tissues/s1/", None],
            "tissues_sample_terms_term_name": ["blood", None],
            "cell_lines_@id": [None, "/cell_lines/s2/"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", TISSUE_COLUMN] == "blood"
    assert sra_df.loc["LIB_B_GEX", TISSUE_COLUMN] == "not available"


def test_tissue_skips_rows_with_no_library_alias():
    main_df = tissue_main_df("tissues")
    main_df["tissues_sample_terms_term_name"] = ["blood", "lung"]
    main_df.loc[1, "droplet_based_libraries_aliases"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # The unkeyed row's term must not leak into the surviving library.
    assert list(sra_df[TISSUE_COLUMN]) == ["blood"]


def test_tissue_is_no_longer_a_prop_map_rename():
    assert "*tissue" not in PROP_MAP_SRA_BIOSAMPLE.values()
    assert "tissues_sample_terms_term_name" not in PROP_MAP_SRA_BIOSAMPLE


# *biomaterial_provider


def test_sample_url_prefixes_are_derived_from_the_tissue_type_map():
    assert sorted(SAMPLE_URL_PREFIXES) == [
        "cell_lines",
        "organoids",
        "primary_cell_cultures",
        "tissues",
    ]
    assert len(SAMPLE_URL_PREFIXES) == len(TISSUE_TYPE_MAP)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param(LAB, ["Alex Marson, UCSF"], id="embedded-dict"),
        pytest.param([SOURCE], ["Abcam"], id="list-of-one-dict"),
        pytest.param([SOURCE, LAB], ["Abcam", "Alex Marson, UCSF"], id="list-of-dicts"),
        pytest.param({"name": "abcam"}, ["abcam"], id="name-when-no-title"),
        pytest.param({"title": "  Abcam  "}, ["Abcam"], id="whitespace-stripped"),
        pytest.param({"title": ""}, [], id="blank-title"),
        pytest.param("/labs/alex-marson/", [], id="bare-id-string-has-no-title"),
        pytest.param(["/sources/abcam/"], [], id="list-of-bare-id-strings"),
        pytest.param({"@id": "/labs/x/"}, [], id="dict-without-title-or-name"),
        pytest.param(None, [], id="none"),
        pytest.param(float("nan"), [], id="nan"),
        pytest.param([], [], id="empty-list"),
    ],
)
def test_provider_titles(value, expected):
    assert make_flattener()._provider_titles(value) == expected


def provider_main_df(**columns):
    base = {
        "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2,
        "tissues_lab": [LAB, LAB],
    }
    base.update(columns)
    return pd.DataFrame(base)


def test_provider_falls_back_to_sample_lab_when_sources_is_absent():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Alex Marson, UCSF"


def test_provider_prefers_sources_over_lab():
    main_df = provider_main_df(tissues_sources=[[SOURCE], [SOURCE]])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Abcam"


def test_provider_falls_back_per_row_not_per_column():
    main_df = provider_main_df(tissues_sources=[[SOURCE], None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Abcam; Alex Marson, UCSF"


def test_provider_falls_back_when_sources_is_an_unresolved_id_path():
    main_df = provider_main_df(tissues_sources=[["/sources/abcam/"], ["/sources/abcam/"]])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Alex Marson, UCSF"


@pytest.mark.parametrize(
    "sample_type", ["tissues", "cell_lines", "organoids", "primary_cell_cultures"]
)
def test_provider_reads_any_sample_type(sample_type):
    main_df = provider_main_df().drop(columns=["tissues_lab"])
    main_df[f"{sample_type}_lab"] = [LAB, LAB]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Alex Marson, UCSF"


def test_provider_ignores_lab_columns_on_non_sample_objects():
    """'lab' is on 19 object types - libraries, files, donors - only samples count."""
    main_df = provider_main_df().drop(columns=["tissues_lab"])
    main_df["droplet_based_libraries_lab"] = [LAB, LAB]
    main_df["human_donors_lab"] = [LAB, LAB]
    main_df["sequence_files_lab"] = [LAB, LAB]
    main_df["treatments_lab"] = [LAB, LAB]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert PROVIDER_COLUMN not in sra_df.columns


def test_provider_skips_rows_with_no_library_alias():
    main_df = provider_main_df(tissues_lab=[LAB, OTHER_LAB])
    main_df.loc[1, "droplet_based_libraries_aliases"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # The unkeyed row's lab must not leak into the surviving library.
    assert list(sra_df[PROVIDER_COLUMN]) == ["Alex Marson, UCSF"]


def test_no_sources_or_lab_column_omits_the_provider(capsys):
    main_df = provider_main_df().drop(columns=["tissues_lab"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert PROVIDER_COLUMN not in sra_df.columns
    assert "no sample sources or lab column" in capsys.readouterr().out


def test_provider_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "tissues_lab": [LAB, LAB, OTHER_LAB],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", PROVIDER_COLUMN] == "Alex Marson, UCSF"
    assert sra_df.loc["LIB_B_GEX", PROVIDER_COLUMN] == "Other Lab"


# *collection_date


def test_collection_date_uses_date_obtained():
    main_df = provider_main_df(tissues_date_obtained=["2023-01-05", "2023-01-05"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05"


def test_collection_date_defaults_when_the_column_is_absent():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert sra_df.loc[0, DATE_COLUMN] == "not provided"


def test_collection_date_defaults_when_the_column_is_all_null():
    main_df = provider_main_df(tissues_date_obtained=[None, None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "not provided"


def test_collection_date_mixes_a_date_with_the_default():
    main_df = provider_main_df(tissues_date_obtained=["2023-01-05", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05; not provided"


def test_collection_date_joins_distinct_dates_in_order():
    main_df = provider_main_df(tissues_date_obtained=["2023-06-01", "2023-01-05"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05; 2023-06-01"


def test_collection_date_default_sorts_last_even_before_a_letter_date():
    """The literal is appended, not sorted in, so a non-ISO date cannot displace it."""
    main_df = provider_main_df(tissues_date_obtained=["osmotic-era", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "osmotic-era; not provided"


def test_collection_date_strips_whitespace_and_dedupes():
    main_df = provider_main_df(tissues_date_obtained=["  2023-01-05  ", "2023-01-05"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05"


@pytest.mark.parametrize(
    "sample_type", ["tissues", "cell_lines", "organoids", "primary_cell_cultures"]
)
def test_collection_date_reads_any_sample_type(sample_type):
    main_df = provider_main_df()
    main_df[f"{sample_type}_date_obtained"] = ["2023-01-05", "2023-01-05"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05"


def test_collection_date_ignores_date_obtained_on_non_sample_objects():
    main_df = provider_main_df()
    main_df["droplet_based_libraries_date_obtained"] = ["2023-01-05", "2023-01-05"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "not provided"


def test_collection_date_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "tissues_lab": [LAB, LAB, LAB],
            "tissues_date_obtained": ["2023-01-05", None, "2024-02-02"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", DATE_COLUMN] == "2023-01-05; not provided"
    assert sra_df.loc["LIB_B_GEX", DATE_COLUMN] == "2024-02-02"


def test_collection_date_skips_rows_with_no_library_alias():
    main_df = provider_main_df(tissues_date_obtained=["2023-01-05", "2024-02-02"])
    main_df.loc[1, "droplet_based_libraries_aliases"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # The unkeyed row's date must not leak into the surviving library.
    assert list(sra_df[DATE_COLUMN]) == ["2023-01-05"]


# *geo_loc_name


def test_geo_loc_name_uses_collection_geographical_location():
    main_df = provider_main_df(tissues_collection_geographical_location=["USA", "USA"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "USA"


def test_geo_loc_name_defaults_when_the_column_is_absent():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert sra_df.loc[0, GEO_COLUMN] == "not provided"


def test_geo_loc_name_mixes_a_location_with_the_default():
    main_df = provider_main_df(tissues_collection_geographical_location=["USA", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "USA; not provided"


def test_geo_loc_name_joins_distinct_locations():
    main_df = provider_main_df(tissues_collection_geographical_location=["USA", "Canada"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "Canada; USA"


@pytest.mark.parametrize(
    "sample_type", ["tissues", "cell_lines", "organoids", "primary_cell_cultures"]
)
def test_geo_loc_name_reads_any_sample_type(sample_type):
    main_df = provider_main_df()
    main_df[f"{sample_type}_collection_geographical_location"] = ["USA", "USA"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "USA"


def test_geo_loc_name_ignores_non_sample_objects():
    main_df = provider_main_df()
    main_df["droplet_based_libraries_collection_geographical_location"] = ["USA", "USA"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "not provided"


def test_geo_loc_name_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "tissues_collection_geographical_location": ["USA", None, "Canada"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", GEO_COLUMN] == "USA; not provided"
    assert sra_df.loc["LIB_B_GEX", GEO_COLUMN] == "Canada"


# ethnicity


def test_ethnicity_is_pooled_and_labelled_like_isolate():
    main_df = donor_main_df(
        human_donors_ethnicity_term_name=[
            "African American",
            "European American",
            "African American",
            "European American",
        ]
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # D1 is 889023040, the European American donor - same labels as *isolate.
    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 889023040, D2 - 889081306"
    assert sra_df.loc[0, ETHNICITY_COLUMN] == (
        "pooled: D1 - European American, D2 - African American"
    )


def test_single_donor_ethnicity_has_no_pooled_prefix():
    main_df = donor_main_df(
        human_donors_cxg_donor_id=[889023040] * 4,
        human_donors_ethnicity_term_name=["European American"] * 4,
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ETHNICITY_COLUMN] == "European American"


def test_donor_without_an_ethnicity_keeps_its_label():
    main_df = donor_main_df(
        human_donors_ethnicity_term_name=[None, "European American", None, "European American"]
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # D2 drops out, D1 keeps its number so it still lines up with *isolate.
    assert sra_df.loc[0, ETHNICITY_COLUMN] == "pooled: D1 - European American"


def test_no_ethnicity_column_omits_the_ethnicity_column():
    sra_df = make_flattener().create_sra_biosample_dataframe(donor_main_df())

    assert ETHNICITY_COLUMN not in sra_df.columns


def test_all_null_ethnicity_omits_the_column():
    main_df = donor_main_df(human_donors_ethnicity_term_name=[None] * 4)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert ETHNICITY_COLUMN not in sra_df.columns


def test_ethnicity_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "human_donors_cxg_donor_id": [11, 22, 33],
            "human_donors_ethnicity_term_name": ["Asian", "European American", "Asian"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", ETHNICITY_COLUMN] == (
        "pooled: D1 - Asian, D2 - European American"
    )
    assert sra_df.loc["LIB_B_GEX", ETHNICITY_COLUMN] == "Asian"


def test_ethnicity_has_no_required_marker():
    """Optional columns carry no '*'; the required donor columns do."""
    main_df = donor_main_df(human_donors_ethnicity_term_name=["Asian"] * 4)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert ETHNICITY_COLUMN in sra_df.columns
    assert f"*{ETHNICITY_COLUMN}" not in sra_df.columns
    assert ISOLATE_COLUMN.startswith("*")


# experimental_perturbation


@pytest.mark.parametrize(
    ("lower", "upper", "units", "expected"),
    [
        pytest.param(8, 8, "hour", "8 hour", id="equal-bounds"),
        pytest.param(8, 24, "hour", "8-24 hour", id="unequal-bounds"),
        pytest.param(8.0, 8.0, "hour", "8 hour", id="float-loses-the-dot-zero"),
        pytest.param(8, None, "hour", "8 hour", id="lower-only"),
        pytest.param(None, 24, "hour", "24 hour", id="upper-only"),
        pytest.param(8, 8, None, "8", id="no-units"),
        # verbatim, not pluralised - '8 hour stimulation' reads adjectivally
        pytest.param(2, 2, "hour", "2 hour", id="units-not-pluralised"),
        pytest.param(None, None, "hour", "", id="no-bounds"),
        pytest.param(None, None, None, "", id="nothing"),
    ],
)
def test_duration_text(lower, upper, units, expected):
    assert make_flattener()._duration_text(lower, upper, units) == expected


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        pytest.param(["8 hour stimulation"], "8 hour stimulation", id="one-is-bare"),
        pytest.param(["a", "b"], "pooled: a, b", id="several-are-prefixed"),
        pytest.param(["a", "", None], "a", id="empties-dropped"),
        pytest.param([], None, id="empty"),
        pytest.param([None, ""], None, id="all-empty"),
    ],
)
def test_format_pooled_values(values, expected):
    assert make_flattener()._format_pooled_values(values) == expected


def perturbation_main_df(rows=2, **columns):
    base = {
        "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * rows,
        "treatments_lower_bound_duration": [8] * rows,
        "treatments_upper_bound_duration": [8] * rows,
        "treatments_duration_units": ["hour"] * rows,
        "treatments_description": ["stimulation"] * rows,
    }
    base.update(columns)
    return pd.DataFrame(base)


def test_perturbation_single_value_is_bare():
    sra_df = make_flattener().create_sra_biosample_dataframe(perturbation_main_df())

    assert sra_df.loc[0, PERTURBATION_COLUMN] == "8 hour stimulation"


def test_perturbation_partly_treated_library_marks_the_gap():
    """The Treg shape: half the samples treated, half not."""
    main_df = perturbation_main_df(
        treatments_lower_bound_duration=[8, None],
        treatments_upper_bound_duration=[8, None],
        treatments_duration_units=["hour", None],
        treatments_description=["stimulation", None],
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PERTURBATION_COLUMN] == "pooled: 8 hour stimulation, not provided"


def test_perturbation_pools_two_treatments():
    main_df = perturbation_main_df(
        treatments_lower_bound_duration=[8, 24],
        treatments_upper_bound_duration=[8, 24],
        treatments_description=["stimulation", "fasting"],
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PERTURBATION_COLUMN] == "pooled: 24 hour fasting, 8 hour stimulation"


def test_perturbation_unequal_bounds():
    main_df = perturbation_main_df(treatments_upper_bound_duration=[24, 24])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PERTURBATION_COLUMN] == "8-24 hour stimulation"


def test_perturbation_description_alone():
    main_df = perturbation_main_df(
        treatments_lower_bound_duration=[None, None],
        treatments_upper_bound_duration=[None, None],
        treatments_duration_units=[None, None],
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PERTURBATION_COLUMN] == "stimulation"


def test_no_treatments_columns_omits_the_perturbation():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert PERTURBATION_COLUMN not in sra_df.columns


def test_all_null_treatments_omits_the_perturbation():
    main_df = perturbation_main_df(
        treatments_lower_bound_duration=[None, None],
        treatments_upper_bound_duration=[None, None],
        treatments_duration_units=[None, None],
        treatments_description=[None, None],
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert PERTURBATION_COLUMN not in sra_df.columns


def test_perturbation_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "treatments_lower_bound_duration": [8, None, 24],
            "treatments_upper_bound_duration": [8, None, 24],
            "treatments_duration_units": ["hour", None, "hour"],
            "treatments_description": ["stimulation", None, "fasting"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", PERTURBATION_COLUMN] == (
        "pooled: 8 hour stimulation, not provided"
    )
    assert sra_df.loc["LIB_B_GEX", PERTURBATION_COLUMN] == "24 hour fasting"


def test_perturbation_skips_rows_with_no_library_alias():
    main_df = perturbation_main_df(treatments_description=["stimulation", "fasting"])
    main_df.loc[1, "droplet_based_libraries_aliases"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # The unkeyed row's treatment must not leak into the surviving library.
    assert list(sra_df[PERTURBATION_COLUMN]) == ["8 hour stimulation"]


# experimental_perturbation_factors


def factors_main_df(rows=2, **columns):
    base = {
        "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * rows,
        "treatments_ontological_term_term_name": [["anti-CD2_HUMAN", "IL2_HUMAN"]] * rows,
    }
    base.update(columns)
    return pd.DataFrame(base)


def test_factors_bracket_a_multi_factor_sample():
    sra_df = make_flattener().create_sra_biosample_dataframe(factors_main_df())

    assert sra_df.loc[0, FACTORS_COLUMN] == "[IL2_HUMAN, anti-CD2_HUMAN]"


def test_factors_single_factor_is_not_bracketed():
    main_df = factors_main_df(treatments_ontological_term_term_name=["IL2_HUMAN"] * 2)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, FACTORS_COLUMN] == "IL2_HUMAN"


def test_factors_untreated_sample_contributes_na():
    """The Treg shape: half the samples treated, half not."""
    main_df = factors_main_df(
        treatments_ontological_term_term_name=[["anti-CD2_HUMAN", "IL2_HUMAN"], None]
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, FACTORS_COLUMN] == "pooled: [IL2_HUMAN, anti-CD2_HUMAN], na"


def test_factors_pool_two_distinct_sets():
    main_df = factors_main_df(
        treatments_ontological_term_term_name=[["anti-CD2_HUMAN", "IL2_HUMAN"], ["IL6_HUMAN"]]
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, FACTORS_COLUMN] == "pooled: IL6_HUMAN, [IL2_HUMAN, anti-CD2_HUMAN]"


def test_factors_deduplicate_matching_sets():
    """Two samples with the same factors give one entry, not a pooled pair."""
    sra_df = make_flattener().create_sra_biosample_dataframe(factors_main_df(rows=4))

    assert sra_df.loc[0, FACTORS_COLUMN] == "[IL2_HUMAN, anti-CD2_HUMAN]"


def test_factors_sort_is_case_sensitive():
    """Deliberate: uppercase sorts ahead of lowercase, matching plain sorted()."""
    main_df = factors_main_df(
        treatments_ontological_term_term_name=[["anti-CD2_HUMAN", "IL2_HUMAN"]] * 2
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, FACTORS_COLUMN].startswith("[IL2_HUMAN,")


def test_no_ontological_term_column_omits_the_factors():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert FACTORS_COLUMN not in sra_df.columns


def test_all_null_ontological_terms_omits_the_factors():
    main_df = factors_main_df(treatments_ontological_term_term_name=[None, None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert FACTORS_COLUMN not in sra_df.columns


def test_factors_are_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "treatments_ontological_term_term_name": [
                ["anti-CD2_HUMAN", "IL2_HUMAN"],
                None,
                ["IL6_HUMAN"],
            ],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", FACTORS_COLUMN] == "pooled: [IL2_HUMAN, anti-CD2_HUMAN], na"
    assert sra_df.loc["LIB_B_GEX", FACTORS_COLUMN] == "IL6_HUMAN"


def test_factors_skip_rows_with_no_library_alias():
    main_df = factors_main_df(treatments_ontological_term_term_name=[["IL2_HUMAN"], ["IL6_HUMAN"]])
    main_df.loc[1, "droplet_based_libraries_aliases"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # The unkeyed row's factor must not leak into the surviving library.
    assert list(sra_df[FACTORS_COLUMN]) == ["IL2_HUMAN"]


# suspension_type


def test_suspension_type_uses_the_sample_field():
    main_df = provider_main_df(tissues_suspension_type=["cell", "cell"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SUSPENSION_COLUMN] == "cell"


def test_suspension_type_joins_distinct_values():
    main_df = provider_main_df(tissues_suspension_type=["nucleus", "cell"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SUSPENSION_COLUMN] == "cell; nucleus"


def test_suspension_type_fills_a_gap_like_a_required_column():
    """Being optional decides only whether the column exists, not how gaps fill."""
    main_df = provider_main_df(tissues_suspension_type=["cell", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SUSPENSION_COLUMN] == "cell; not provided"


def test_a_library_with_no_suspension_type_still_gets_a_cell():
    """One library has values so the column exists; the other is filled, not blank."""
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]]
            + [["alex-marson:LIB_B_GEX"]],
            "tissues_suspension_type": ["cell", None],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", SUSPENSION_COLUMN] == "cell"
    assert sra_df.loc["LIB_B_GEX", SUSPENSION_COLUMN] == "not provided"


def test_no_suspension_type_column_omits_it():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert SUSPENSION_COLUMN not in sra_df.columns


def test_all_null_suspension_type_omits_it():
    main_df = provider_main_df(tissues_suspension_type=[None, None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert SUSPENSION_COLUMN not in sra_df.columns


@pytest.mark.parametrize(
    "sample_type", ["tissues", "cell_lines", "organoids", "primary_cell_cultures"]
)
def test_suspension_type_reads_any_sample_type(sample_type):
    main_df = provider_main_df()
    main_df[f"{sample_type}_suspension_type"] = ["cell", "cell"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SUSPENSION_COLUMN] == "cell"


def test_suspension_type_is_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "tissues_suspension_type": ["cell", "nucleus", "cell"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", SUSPENSION_COLUMN] == "cell; nucleus"
    assert sra_df.loc["LIB_B_GEX", SUSPENSION_COLUMN] == "cell"


# age_lower_bound / age_upper_bound


@pytest.mark.parametrize(
    ("value", "units", "expected"),
    [
        pytest.param(29, "year", "29 years", id="int"),
        pytest.param(29.0, "year", "29 years", id="float-loses-the-dot-zero"),
        pytest.param(1, "year", "1 year", id="singular"),
        pytest.param(29, "years", "29 years", id="already-plural-units"),
        pytest.param(6, "month", "6 months", id="months"),
        pytest.param(1, "day", "1 day", id="singular-day"),
        pytest.param(29.5, "year", "29.5 years", id="non-integer-kept"),
        pytest.param(29, None, "29", id="no-units"),
        pytest.param(29, "", "29", id="blank-units"),
        pytest.param(None, "year", None, id="no-value"),
        pytest.param(float("nan"), "year", None, id="nan-value"),
    ],
)
def test_age_with_units(value, units, expected):
    assert age_with_units(value, units) == expected


def age_main_df(**columns):
    """
    donor_main_df plus age bounds, keyed by donor so the ordering is checkable:
    889023040 is D1 (the lower id) at 32-35, 889081306 is D2 at 29-30.
    """
    frame = donor_main_df()
    frame["tissues_lower_bound_age"] = frame["human_donors_cxg_donor_id"].map(
        {889023040: 32, 889081306: 29}
    )
    frame["tissues_upper_bound_age"] = frame["human_donors_cxg_donor_id"].map(
        {889023040: 35, 889081306: 30}
    )
    frame["tissues_age_units"] = "year"
    for column, value in columns.items():
        frame[column] = value
    return frame


def test_age_bounds_follow_donor_order_so_entry_n_matches_dn():
    sra_df = make_flattener().create_sra_biosample_dataframe(age_main_df())

    # D1 is 889023040, the 32-year-old, so it leads all three columns.
    assert sra_df.loc[0, AGE_COLUMN] == "pooled: D1 - 32 years, D2 - 29 years"
    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "32 years; 29 years"
    assert sra_df.loc[0, AGE_UPPER_COLUMN] == "35 years; 30 years"


def test_age_bounds_are_not_sorted():
    """A younger second donor must not be reordered to the front."""
    main_df = age_main_df()
    main_df["tissues_lower_bound_age"] = main_df["human_donors_cxg_donor_id"].map(
        {889023040: 29, 889081306: 3}
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "29 years; 3 years"


def test_age_bounds_do_not_dedupe_matching_donors():
    """Two donors at the same age must still give two entries."""
    main_df = age_main_df(tissues_lower_bound_age=29)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "29 years; 29 years"


def test_age_bounds_fill_a_donor_with_no_value():
    """Three donors, only the third with a bound - the slots hold."""
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 3,
            "human_donors_cxg_donor_id": [1, 2, 3],
            "tissues_lower_bound_age": [None, None, 40],
            "tissues_age_units": ["year"] * 3,
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: D1 - 1, D2 - 2, D3 - 3"
    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "not provided; not provided; 40 years"
    # No donor has an upper bound, so that column is dropped.
    assert AGE_UPPER_COLUMN not in sra_df.columns


def test_no_age_bounds_omits_both_columns():
    sra_df = make_flattener().create_sra_biosample_dataframe(donor_main_df())

    assert AGE_LOWER_COLUMN not in sra_df.columns
    assert AGE_UPPER_COLUMN not in sra_df.columns


def test_lower_bound_alone_omits_only_the_upper_column():
    main_df = age_main_df().drop(columns=["tissues_upper_bound_age"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "32 years; 29 years"
    assert AGE_UPPER_COLUMN not in sra_df.columns


def test_all_null_bound_omits_that_column():
    main_df = age_main_df(tissues_lower_bound_age=None)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert AGE_LOWER_COLUMN not in sra_df.columns
    assert sra_df.loc[0, AGE_UPPER_COLUMN] == "35 years; 30 years"


def test_age_bounds_read_primary_cell_cultures_too():
    main_df = age_main_df().rename(
        columns={
            "tissues_lower_bound_age": "primary_cell_cultures_lower_bound_age",
            "tissues_upper_bound_age": "primary_cell_cultures_upper_bound_age",
            "tissues_age_units": "primary_cell_cultures_age_units",
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "32 years; 29 years"


def test_age_bounds_without_units():
    main_df = age_main_df(tissues_age_units=None)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "32; 29"


def test_age_bounds_are_scoped_to_each_library():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_aliases": [["alex-marson:LIB_A_GEX"]] * 2
            + [["alex-marson:LIB_B_GEX"]],
            "human_donors_cxg_donor_id": [11, 22, 33],
            "tissues_lower_bound_age": [29, 22, 40],
            "tissues_age_units": ["year"] * 3,
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A_GEX", AGE_LOWER_COLUMN] == "29 years; 22 years"
    assert sra_df.loc["LIB_B_GEX", AGE_LOWER_COLUMN] == "40 years"


def test_age_bounds_need_a_donor_id_column():
    """The bounds are keyed per donor, so with no donor id there is nothing to key."""
    main_df = age_main_df().drop(columns=["human_donors_cxg_donor_id"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert AGE_LOWER_COLUMN not in sra_df.columns
    assert AGE_UPPER_COLUMN not in sra_df.columns


def test_prop_map_sends_both_library_types_to_sample_name():
    """The rename relies on both library keys sharing one output name."""
    assert PROP_MAP_SRA_BIOSAMPLE["droplet_based_libraries_aliases"] == "sample_name"
    assert PROP_MAP_SRA_BIOSAMPLE["plate_based_libraries_aliases"] == "sample_name"
    # The per-sample alias is deliberately absent: the key is the library.
    assert "sample_alias" not in PROP_MAP_SRA_BIOSAMPLE
    # sample_name is the one output column without the SRA required marker.
    assert not PROP_MAP_SRA_BIOSAMPLE["droplet_based_libraries_aliases"].startswith("*")
    assert ORGANISM_COLUMN.startswith("*")


def test_prop_map_holds_only_renames():
    """The derived columns are named in the flattener, not mapped from MAIN."""
    assert "isolate" not in PROP_MAP_SRA_BIOSAMPLE
    assert "age" not in PROP_MAP_SRA_BIOSAMPLE
    assert ISOLATE_COLUMN not in PROP_MAP_SRA_BIOSAMPLE.values()
    assert AGE_COLUMN not in PROP_MAP_SRA_BIOSAMPLE.values()
    assert BARCODE_COLUMN not in PROP_MAP_SRA_BIOSAMPLE.values()
