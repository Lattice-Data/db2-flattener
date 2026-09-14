"""
SRA/BioSample sheet: one row per CRO group, keyed on the library's
CRO_group_identifier under the name 'sample_name'. That identifier spans a
group's libraries, so only the GEX ones reach the sheet.
"""

import pandas as pd
import pytest

from db2_flattener.flatten.flattener import (
    SAMPLE_URL_PREFIXES,
    DB2Flattener,
    age_with_units,
)
from db2_flattener.schema.constants import (
    GENETIC_PERTURBATION_MAP,
    PROP_MAP_SRA_BIOSAMPLE,
    TISSUE_TYPE_MAP,
    Configs,
)
from db2_flattener.utils import ages_from_developmental_stages

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
INTENDED_CELL_TYPE_COLUMN = "intended_cell_type"
ENRICHED_COLUMN = "suspension_enriched_cell_types"
DEPLETED_COLUMN = "suspension_depleted_cell_types"
KITS_COLUMN = "suspension_selection_kits"
ENRICHMENT_FACTORS_COLUMN = "suspension_enrichment_factors"
STRATEGY_COLUMN = "genetic_perturbation_strategy"
PRESERVATION_COLUMN = "preservation_method"
# Read off the map so a rename there, e.g. dropping the SRA '*' required marker,
# does not have to be chased through every assertion below.
ORGANISM_COLUMN = PROP_MAP_SRA_BIOSAMPLE["human_donors_taxa"]

LAB = {"@id": "/labs/alex-marson/", "title": "Alex Marson, UCSF"}
SOURCE = {"@id": "/sources/abcam/", "title": "Abcam"}


def make_flattener():
    """Empty configs on purpose: this sheet is built from main_df alone."""
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = None
    flattener.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    return flattener


def main_frame(columns):
    """
    A MAIN frame whose droplet libraries declare themselves GEX.

    The sheet takes GEX libraries only, and a droplet library that leaves
    feature_types unset is read as non-GEX, so a fixture that did not say so
    would be filtered out before any cell was built. A plate library needs no
    such marker: unset feature_types is read as GEX there.
    """
    frame = pd.DataFrame(columns)
    if "droplet_based_libraries_CRO_group_identifier" in frame.columns:
        frame["droplet_based_libraries_feature_types"] = [["Gene Expression"]] * len(frame)
    return frame


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


# create_sra_biosample_dataframe


def droplet_main_df():
    """Four samples pooled into one droplet library, one MAIN row per sample."""
    aliases = [s["aliases"][0].split(":", 1)[1] for s in FOUR_SAMPLES]
    return main_frame(
        {
            "sample_alias": aliases,
            "droplet_based_libraries_CRO_group_identifier": ["TregR3_L13_L05"] * 4,
            "droplet_based_libraries_samples": [FOUR_SAMPLES] * 4,
            "human_donors_taxa": ["Homo sapiens"] * 4,
        }
    )


def test_droplet_run_is_one_row_per_library(capsys):
    sra_df = make_flattener().create_sra_biosample_dataframe(droplet_main_df())

    # No optional column has data here, so only the required ones are derived,
    # and every one of those is present even where MAIN has nothing to fill it
    assert list(sra_df.columns) == [
        "sample_name",
        ORGANISM_COLUMN,
        BARCODE_COLUMN,
        ISOLATE_COLUMN,
        AGE_COLUMN,
        SEX_COLUMN,
        TISSUE_COLUMN,
        PROVIDER_COLUMN,
        DATE_COLUMN,
        GEO_COLUMN,
    ]
    assert "no sample sources or lab column" in capsys.readouterr().out
    assert len(sra_df) == 1
    assert sra_df.loc[0, "sample_name"] == "TregR3_L13_L05"
    assert sra_df.loc[0, ORGANISM_COLUMN] == "Homo sapiens"
    assert sra_df.loc[0, DATE_COLUMN] == "not provided"
    assert sra_df.loc[0, BARCODE_COLUMN] == FOUR_SAMPLES_MAP


def test_plate_library_columns_are_used_when_droplet_is_absent():
    main_df = main_frame(
        {
            "plate_based_libraries_CRO_group_identifier": ["PLATE_1"],
            "plate_based_libraries_samples": [[sample("s1", ["BC001+CR001"])]],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["PLATE_1"]
    assert list(sra_df[BARCODE_COLUMN]) == ["s1 : BC001+CR001"]


def test_paired_libraries_share_one_row_because_the_crispr_half_is_dropped(capsys):
    """
    A group's GEX and CRISPR libraries carry one CRO id, so they are one row.

    This is what taking GEX libraries only buys: without the filter the pair
    would collapse together anyway, but every cell would be built from twice
    the rows.
    """
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 4,
            "droplet_based_libraries_feature_types": [["Gene Expression"]] * 2
            + [["CRISPR Guide Capture"]] * 2,
            "droplet_based_libraries_samples": [FOUR_SAMPLES] * 4,
            "human_donors_taxa": ["Homo sapiens"] * 4,
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["LIB_A"]
    assert sra_df.loc[0, BARCODE_COLUMN] == FOUR_SAMPLES_MAP
    assert "filtered to 2 GEX rows out of 4 MAIN rows" in capsys.readouterr().out


def test_group_with_no_gex_library_is_dropped_entirely(capsys):
    """A CRISPR-only group has no GEX library, so it leaves the sheet."""
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_B"],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
            ],
            "human_donors_taxa": ["Homo sapiens"] * 2,
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["LIB_A"]
    assert "filtered to 1 GEX rows out of 2 MAIN rows" in capsys.readouterr().out


def test_all_libraries_non_gex_returns_empty_frame():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_feature_types": [["CRISPR Guide Capture"]],
        }
    )

    assert make_flattener().create_sra_biosample_dataframe(main_df).empty


# library_id


def library_id_main_df():
    """Two groups, each a GEX library paired with a CRISPR one."""
    return pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_A", "LIB_B", "LIB_B"],
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_A_CRI"],
                ["alex-marson:LIB_B_GEX"],
                ["alex-marson:LIB_B_CRI"],
            ],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
            ]
            * 2,
        }
    )


def test_library_id_lists_the_gex_library_before_the_crispr_one():
    """
    Sorting alone would lead with the CRISPR library, so GEX is put first.

    The whole column comes from the unfiltered frame: LIB_A_CRI is on a row the
    GEX filter drops, so a builder reading the filtered frame could not see it.
    """
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
            "droplet_based_libraries_aliases": [
                ["alex-marson:LIB_A_GEX"],
                ["alex-marson:LIB_A_CRI"],
            ],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
            ],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["LIB_A"]
    assert list(sra_df["library_id"]) == ["LIB_A_GEX, LIB_A_CRI"]


def test_library_id_is_scoped_to_each_group():
    sra_df = make_flattener().create_sra_biosample_dataframe(library_id_main_df())

    assert list(sra_df["library_id"]) == ["LIB_A_GEX, LIB_A_CRI", "LIB_B_GEX, LIB_B_CRI"]


def test_library_id_strips_the_lab_prefix():
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["some-other-lab:LIB_A_GEX"]],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["library_id"]) == ["LIB_A_GEX"]


def test_missing_aliases_column_omits_library_id(capsys):
    sra_df = make_flattener().create_sra_biosample_dataframe(droplet_main_df())

    assert "library_id" not in sra_df.columns
    assert "no library aliases column" in capsys.readouterr().out


def test_library_rows_disagreeing_on_the_map_collapse_to_a_list():
    other_samples = [sample("s1", ["BC099+CR099"])]
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
            "droplet_based_libraries_samples": [FOUR_SAMPLES, other_samples],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert len(sra_df) == 1
    assert sra_df.loc[0, BARCODE_COLUMN] == [FOUR_SAMPLES_MAP, "s1 : BC099+CR099"]


def test_rows_with_no_library_group_are_dropped(capsys):
    main_df = droplet_main_df()
    main_df.loc[1, "droplet_based_libraries_CRO_group_identifier"] = None

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df["sample_name"]) == ["TregR3_L13_L05"]
    assert (
        "dropping 1 of 4 MAIN row(s) with no library CRO group identifier"
        in capsys.readouterr().out
    )


def test_missing_library_group_column_returns_empty_frame(capsys):
    main_df = main_frame({"sample_alias": ["s1"]})

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.empty
    assert "no library CRO group identifier column" in capsys.readouterr().out


def test_key_only_main_df_collapses_without_aggregating():
    """groupby().agg({}) raises, so a frame of nothing but the key dedupes instead."""
    main_df = main_frame(
        {"droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_A", "LIB_B"]}
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert list(sra_df.columns) == [
        "sample_name",
        ISOLATE_COLUMN,
        AGE_COLUMN,
        SEX_COLUMN,
        TISSUE_COLUMN,
        PROVIDER_COLUMN,
        DATE_COLUMN,
        GEO_COLUMN,
    ]
    assert list(sra_df["sample_name"]) == ["LIB_A", "LIB_B"]
    # A frame of nothing but the key still fills every required column
    for column in (ISOLATE_COLUMN, AGE_COLUMN, SEX_COLUMN, TISSUE_COLUMN, PROVIDER_COLUMN):
        assert list(sra_df[column]) == ["not provided", "not provided"], column
    assert list(sra_df[DATE_COLUMN]) == ["not provided", "not provided"]
    assert list(sra_df[GEO_COLUMN]) == ["not provided", "not provided"]


def test_empty_main_df_returns_empty_frame():
    main_df = main_frame({"droplet_based_libraries_CRO_group_identifier": pd.Series(dtype=object)})

    assert make_flattener().create_sra_biosample_dataframe(main_df).empty


# isolate and age


@pytest.mark.parametrize(
    ("term_name", "expected"),
    [
        pytest.param("29-year-old stage", ["29 years"], id="years"),
        pytest.param("1-year-old stage", ["1 year"], id="singular-year"),
        pytest.param("6-month-old stage", ["6 months"], id="months"),
        pytest.param("1-week-old stage", ["1 week"], id="singular-week"),
        # Every stage contributes: nothing is dropped for being coarser
        pytest.param(
            ["adult stage", "42-year-old stage"],
            ["adult", "42 years"],
            id="qualitative-and-numeric-both-kept",
        ),
        pytest.param(
            ["29-year-old stage", "32-year-old stage"],
            ["29 years", "32 years"],
            id="two-numeric-stages",
        ),
        pytest.param(
            "29-year-old stage; 32-year-old stage",
            ["29 years", "32 years"],
            id="two-numeric-stages-joined-in-one-cell",
        ),
        pytest.param(
            ["29-year-old stage", "29-year-old stage"],
            ["29 years"],
            id="matching-stages-dedupe",
        ),
        pytest.param("adult stage", ["adult"], id="qualitative"),
        pytest.param("newborn stage", ["newborn"], id="newborn"),
        pytest.param("adult", ["adult"], id="no-stage-suffix"),
        pytest.param(
            "10th week post-fertilization human stage",
            ["10th week post-fertilization"],
            id="qualitative-containing-a-number",
        ),
        pytest.param("adult human stage", ["adult"], id="human-stage-suffix"),
        pytest.param("mouse adult stage", ["mouse adult"], id="only-human-is-trimmed"),
        pytest.param(
            ["adult stage", "newborn stage"],
            ["adult", "newborn"],
            id="two-qualitative-stages",
        ),
        pytest.param(None, [], id="none"),
        pytest.param(float("nan"), [], id="nan"),
        pytest.param([], [], id="empty-list"),
        pytest.param("   ", [], id="blank"),
    ],
)
def test_ages_from_developmental_stages(term_name, expected):
    assert ages_from_developmental_stages(term_name) == expected


def donor_main_df(**columns):
    """
    One library, four samples, two donors - one donor per sample row.
    889023040 is the 32-year-old female, 889081306 the 29-year-old male.
    """
    base = {
        "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 4,
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
    return main_frame(base)


def test_donor_columns_are_unordered_pooled_sets():
    """'pooled:' marks the cell as an unordered set, nothing more."""
    sra_df = make_flattener().create_sra_biosample_dataframe(donor_main_df())

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: 889023040, 889081306"
    assert sra_df.loc[0, AGE_COLUMN] == "pooled: 29 years, 32 years"


def test_a_sample_pooling_several_donors_is_split_apart():
    """
    create_dataframe collapses a multi-donor sample into one '; '-joined cell.
    Because the cell is an unordered set there is no pairing to recover, so
    splitting it is enough.
    """
    main_df = donor_main_df(human_donors_cxg_donor_id=["889023040; 889081306"] * 4)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: 889023040, 889081306"


def test_age_pools_alongside_isolate_and_sex_on_a_multi_donor_row():
    """
    A row claiming two donors must not claim one age.

    *age goes through a transform, which bypasses split_joined, so this is the
    case where it used to keep the first stage and disagree with its own row.
    """
    main_df = donor_main_df(
        human_donors_cxg_donor_id=["889023040; 889081306"] * 4,
        human_donors_sex=["female; male"] * 4,
        tissues_developmental_stages_term_name=["29-year-old stage; 32-year-old stage"] * 4,
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: 889023040, 889081306"
    assert sra_df.loc[0, SEX_COLUMN] == "pooled male and female"
    assert sra_df.loc[0, AGE_COLUMN] == "pooled: 29 years, 32 years"


def test_age_keeps_a_qualitative_stage_alongside_a_numeric_one():
    """Nothing is dropped for being coarser: both stages reach the cell."""
    main_df = donor_main_df(
        human_donors_cxg_donor_id=["D1"] * 4,
        tissues_developmental_stages_term_name=[["adult stage", "42-year-old stage"]] * 4,
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_COLUMN] == "pooled: 42 years, adult"


def test_donor_id_is_not_rendered_as_a_float():
    """A null in the column makes pandas store the ids as float64."""
    main_df = donor_main_df()
    main_df.loc[3, "human_donors_cxg_donor_id"] = None
    assert main_df["human_donors_cxg_donor_id"].dtype == "float64"

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert "889023040.0" not in sra_df.loc[0, ISOLATE_COLUMN]


def test_no_donor_id_column_still_fills_isolate():
    main_df = donor_main_df().drop(columns=["human_donors_cxg_donor_id"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    # isolate is required, so it stays as a filled column rather than vanishing
    assert sra_df.loc[0, ISOLATE_COLUMN] == "not provided"


# *sex


@pytest.mark.parametrize(
    ("sexes", "expected"),
    [
        pytest.param(["male"], "male", id="male-only"),
        pytest.param(["female", "male"], "pooled male and female", id="male-listed-first"),
        pytest.param(["male", "male", "female"], "pooled male and female", id="deduplicated"),
        pytest.param(["unknown", "male"], "pooled male and unknown", id="male-then-unknown"),
        pytest.param(
            ["unknown", "female", "male"],
            "pooled male, female and unknown",
            id="three-values",
        ),
        pytest.param(["unknown"], "unknown", id="unknown-alone-passes-through"),
        pytest.param([], None, id="empty"),
    ],
)
def test_format_pooled_sex(sexes, expected):
    assert make_flattener()._format_pooled_sex(sexes) == expected


def test_a_sample_pooling_several_donors_splits_the_sexes():
    """The reviewer's case: 'female; male' in one cell must not read as one sex."""
    main_df = donor_main_df(human_donors_sex=["female; male"] * 4)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "pooled male and female"


def test_missing_sex_column_still_fills_the_required_cell():
    """*sex is required, so it is never blank, even with no sex data at all."""
    main_df = donor_main_df().drop(columns=["human_donors_sex"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "not provided"


def test_partly_known_sex_pool_names_the_gap():
    """One donor sexed and one not is a pool of two states, not one."""
    main_df = donor_main_df(human_donors_sex=["male", None, None, None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "pooled male and not provided"


def test_non_human_donor_columns_are_used_when_human_is_absent():
    main_df = donor_main_df().rename(
        columns={
            "human_donors_sex": "non_human_donors_sex",
            "human_donors_cxg_donor_id": "non_human_donors_cxg_donor_id",
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SEX_COLUMN] == "pooled male and female"
    assert sra_df.loc[0, ISOLATE_COLUMN] == "pooled: 889023040, 889081306"


def test_donor_columns_are_scoped_to_each_library():
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2 + ["LIB_B"],
            "human_donors_cxg_donor_id": [11, 22, 33],
            "human_donors_sex": ["male", "female", "male"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df).set_index("sample_name")

    assert sra_df.loc["LIB_A", ISOLATE_COLUMN] == "pooled: 11, 22"
    assert sra_df.loc["LIB_A", SEX_COLUMN] == "pooled male and female"
    assert sra_df.loc["LIB_B", ISOLATE_COLUMN] == "33"
    assert sra_df.loc["LIB_B", SEX_COLUMN] == "male"


# *tissue


def tissue_main_df(prefix, term="blood", **columns):
    """One library, two samples of the given sample type."""
    base = {
        "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
        f"{prefix}_@id": [f"/{prefix}/s1/", f"/{prefix}/s2/"],
    }
    if term is not None:
        base[f"{prefix}_sample_terms_term_name"] = [term, term]
    base.update(columns)
    return main_frame(base)


@pytest.mark.parametrize("prefix", ["cell_lines", "primary_cell_cultures"])
def test_tissueless_sample_types_report_not_applicable(prefix):
    sra_df = make_flattener().create_sra_biosample_dataframe(tissue_main_df(prefix, "HeLa"))

    assert sra_df.loc[0, TISSUE_COLUMN] == "not applicable"


@pytest.mark.parametrize("prefix", ["tissues", "organoids"])
def test_tissue_and_organoid_report_their_sample_term(prefix):
    sra_df = make_flattener().create_sra_biosample_dataframe(tissue_main_df(prefix))

    assert sra_df.loc[0, TISSUE_COLUMN] == "blood"


def test_tissue_mixes_a_real_term_with_not_applicable():
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
            "tissues_@id": ["/tissues/s1/", None],
            "tissues_sample_terms_term_name": ["blood", None],
            "cell_lines_@id": [None, "/cell_lines/s2/"],
            "cell_lines_sample_terms_term_name": [None, "HeLa"],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, TISSUE_COLUMN] == "blood; not applicable"


def test_tissue_with_no_sample_term_column_still_fills_the_required_cell():
    """Should not happen - the schema requires a term - but *tissue is required."""
    main_df = tissue_main_df("tissues", term=None)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, TISSUE_COLUMN] == "not provided"


# SAMPLE_URL_PREFIXES


def test_sample_url_prefixes_are_derived_from_the_tissue_type_map():
    assert sorted(SAMPLE_URL_PREFIXES) == [
        "cell_lines",
        "organoids",
        "primary_cell_cultures",
        "tissues",
    ]
    assert len(SAMPLE_URL_PREFIXES) == len(TISSUE_TYPE_MAP)


# *biomaterial_provider


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
        "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
        "tissues_lab": [LAB, LAB],
    }
    base.update(columns)
    return main_frame(base)


def test_provider_prefers_sources_over_lab():
    main_df = provider_main_df(tissues_sources=[[SOURCE], [SOURCE]])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Abcam"


def test_provider_falls_back_per_row_not_per_column():
    main_df = provider_main_df(tissues_sources=[[SOURCE], None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "Abcam; Alex Marson, UCSF"


def test_provider_ignores_lab_columns_on_non_sample_objects():
    """'lab' is on 19 object types - libraries, files, donors - only samples count."""
    main_df = provider_main_df().drop(columns=["tissues_lab"])
    main_df["droplet_based_libraries_lab"] = [LAB, LAB]
    main_df["human_donors_lab"] = [LAB, LAB]
    main_df["sequence_files_lab"] = [LAB, LAB]
    main_df["treatments_lab"] = [LAB, LAB]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "not provided"


def test_no_sources_or_lab_column_still_fills_the_required_cell(capsys):
    """Warned about, but a required column is filled rather than dropped."""
    main_df = provider_main_df().drop(columns=["tissues_lab"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PROVIDER_COLUMN] == "not provided"
    assert "no sample sources or lab column" in capsys.readouterr().out


# *collection_date


def test_collection_date_uses_date_obtained():
    main_df = provider_main_df(tissues_date_obtained=["2023-01-05", "2023-01-05"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05"


def test_collection_date_mixes_a_date_with_the_default():
    main_df = provider_main_df(tissues_date_obtained=["2023-01-05", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "2023-01-05; not provided"


def test_collection_date_ignores_date_obtained_on_non_sample_objects():
    main_df = provider_main_df()
    main_df["droplet_based_libraries_date_obtained"] = ["2023-01-05", "2023-01-05"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, DATE_COLUMN] == "not provided"


# *geo_loc_name


def test_geo_loc_name_uses_collection_geographical_location():
    main_df = provider_main_df(tissues_collection_geographical_location=["USA", "USA"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "USA"


def test_geo_loc_name_ignores_non_sample_objects():
    main_df = provider_main_df()
    main_df["droplet_based_libraries_collection_geographical_location"] = ["USA", "USA"]

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, GEO_COLUMN] == "not provided"


# ethnicity


def test_a_donor_without_an_ethnicity_adds_the_gap():
    main_df = donor_main_df(
        human_donors_ethnicity_term_name=[None, "European American", None, "European American"]
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ETHNICITY_COLUMN] == "pooled: European American, not provided"


def test_no_ethnicity_column_omits_the_ethnicity_column():
    sra_df = make_flattener().create_sra_biosample_dataframe(donor_main_df())

    assert ETHNICITY_COLUMN not in sra_df.columns


def test_ethnicity_is_human_only():
    """The field is on HumanDonor, so a non-human run has no such column."""
    main_df = donor_main_df(non_human_donors_ethnicity_term_name=["Asian"] * 4)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert ETHNICITY_COLUMN not in sra_df.columns


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
        "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * rows,
        "treatments_lower_bound_duration": [8] * rows,
        "treatments_upper_bound_duration": [8] * rows,
        "treatments_duration_units": ["hour"] * rows,
        "treatments_description": ["stimulation"] * rows,
    }
    base.update(columns)
    return main_frame(base)


def test_perturbation_partly_treated_library_marks_the_gap():
    """The Treg shape: half the samples treated, half not."""
    main_df = perturbation_main_df(
        treatments_lower_bound_duration=[8, None],
        treatments_upper_bound_duration=[8, None],
        treatments_duration_units=["hour", None],
        treatments_description=["stimulation", None],
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PERTURBATION_COLUMN] == "pooled: 8 hour stimulation, no treatment"


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


def test_perturbation_warns_when_a_sample_has_several_treatments(capsys):
    """
    create_dataframe collapses them into one cell with the duration and the
    description sorted independently, so the pairing is gone. Warn, do not guess.
    """
    main_df = perturbation_main_df(
        treatments_lower_bound_duration=["8; 24"] * 2,
        treatments_upper_bound_duration=["8; 24"] * 2,
        treatments_description=["stimulation; washout"] * 2,
    )

    make_flattener().create_sra_biosample_dataframe(main_df)

    assert "several treatments collapsed into one cell" in capsys.readouterr().out


def test_perturbation_does_not_warn_for_a_single_treatment(capsys):
    make_flattener().create_sra_biosample_dataframe(perturbation_main_df())

    assert "several treatments" not in capsys.readouterr().out


# experimental_perturbation_factors


def factors_main_df(rows=2, **columns):
    base = {
        "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * rows,
        "treatments_ontological_term_term_name": [["anti-CD2_HUMAN", "IL2_HUMAN"]] * rows,
    }
    base.update(columns)
    return main_frame(base)


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


# preservation_method


def test_preservation_method_reads_tissues():
    main_df = provider_main_df(tissues_preservation_method=["fresh", "fresh"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PRESERVATION_COLUMN] == "fresh"


@pytest.mark.parametrize("prefix", ["cell_lines", "organoids", "primary_cell_cultures"])
def test_preservation_method_ignores_other_sample_types(prefix):
    """The field is on tissues alone, so the column is scoped to that prefix."""
    main_df = provider_main_df(**{f"{prefix}_preservation_method": ["fresh", "fresh"]})

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert PRESERVATION_COLUMN not in sra_df.columns


def test_preservation_method_gap_is_not_applicable():
    main_df = provider_main_df(tissues_preservation_method=["fresh", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, PRESERVATION_COLUMN] == "fresh; not applicable"


def test_preservation_method_absent_when_no_sample_has_one():
    main_df = provider_main_df(tissues_preservation_method=[None, None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert PRESERVATION_COLUMN not in sra_df.columns


# genetic_perturbation_strategy


def strategy_main_df(values=("interference screen", "interference screen")):
    return main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * len(values),
            "genetic_modifications_strategy": list(values),
        }
    )


@pytest.mark.parametrize(
    ("stored", "expected"),
    [
        pytest.param("interference screen", "CRISPR interference screen", id="interference"),
        pytest.param("activation screen", "CRISPR activation screen", id="activation"),
        pytest.param("knockout mutation", "CRISPR knockout mutant", id="knockout-mutation"),
        pytest.param("knockout screen", "CRISPR knockout screen", id="knockout-screen"),
        # anything the map does not name passes through untouched
        pytest.param("base editing", "base editing", id="unmapped-passes-through"),
    ],
)
def test_genetic_strategy_is_translated(stored, expected):
    sra_df = make_flattener().create_sra_biosample_dataframe(strategy_main_df((stored, stored)))

    assert sra_df.loc[0, STRATEGY_COLUMN] == expected


def test_genetic_strategy_matches_the_map_biohub_uses():
    """Both sheets report this field, so they must agree on the wording."""
    for stored, expected in GENETIC_PERTURBATION_MAP.items():
        sra_df = make_flattener().create_sra_biosample_dataframe(strategy_main_df((stored, stored)))
        assert sra_df.loc[0, STRATEGY_COLUMN] == expected


def test_genetic_strategy_ignores_sample_prefixed_columns():
    """'strategy' is read off genetic_modifications, not off the sample."""
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
            "tissues_strategy": ["interference screen"] * 2,
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert STRATEGY_COLUMN not in sra_df.columns


# intended_cell_type / the suspension_* columns


def cell_type_main_df(prefix="cell_lines", values=(["HeLa"], ["HeLa"]), field="intended"):
    return main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * len(values),
            f"{prefix}_{field}_cell_types_term_name": list(values),
        }
    )


@pytest.mark.parametrize("prefix", ["cell_lines", "organoids"])
def test_intended_cell_type_reads_intended_cell_types(prefix):
    sra_df = make_flattener().create_sra_biosample_dataframe(cell_type_main_df(prefix))

    assert sra_df.loc[0, INTENDED_CELL_TYPE_COLUMN] == "HeLa"


@pytest.mark.parametrize("prefix", ["tissues", "primary_cell_cultures"])
def test_intended_cell_type_ignores_sample_types_without_intended_cell_types(prefix):
    """The field is only on cell lines and organoids."""
    sra_df = make_flattener().create_sra_biosample_dataframe(cell_type_main_df(prefix))

    assert INTENDED_CELL_TYPE_COLUMN not in sra_df.columns


def test_both_cell_type_columns_can_coexist():
    """A cell line run carries both, each with its own gap marker."""
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * 2,
            "cell_lines_intended_cell_types_term_name": [["HeLa"], None],
            "cell_lines_enriched_cell_types_term_name": [["T cell"], None],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, INTENDED_CELL_TYPE_COLUMN] == "HeLa; not applicable"
    assert sra_df.loc[0, ENRICHED_COLUMN] == "T cell; not provided"


def test_enriched_and_depleted_are_separate_columns():
    """One sample enriched for a type and depleted of another reports both."""
    main_df = main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "tissues_enriched_cell_types_term_name": [["T cell"]],
            "tissues_depleted_cell_types_term_name": [["B cell"]],
        }
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ENRICHED_COLUMN] == "T cell"
    assert sra_df.loc[0, DEPLETED_COLUMN] == "B cell"


# suspension_selection_kits / suspension_enrichment_factors


def selection_main_df(field, values, prefix="tissues"):
    return main_frame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"] * len(values),
            f"{prefix}_{field}": list(values),
        }
    )


@pytest.mark.parametrize("prefix", ["tissues", "cell_lines", "organoids", "primary_cell_cultures"])
def test_selection_kits_reads_any_sample_type(prefix):
    main_df = selection_main_df("selection_kits", (["EasySep CD4"],), prefix)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, KITS_COLUMN] == "EasySep CD4"


@pytest.mark.parametrize("prefix", ["tissues", "cell_lines", "organoids", "primary_cell_cultures"])
def test_enrichment_factors_read_selection_markers(prefix):
    """BIOHUB's mapping: the factors column is fed by selection_markers."""
    main_df = selection_main_df("selection_markers", (["lipid stain (top 30%)"],), prefix)

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ENRICHMENT_FACTORS_COLUMN] == "lipid stain (top 30%)"


def test_enrichment_factors_pool_across_samples_that_disagree():
    """One group whose samples were sorted on opposite tails of the same stain."""
    main_df = selection_main_df(
        "selection_markers", (["lipid stain (top 30%)"], ["lipid stain (bottom 30%)"])
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, ENRICHMENT_FACTORS_COLUMN] == (
        "lipid stain (bottom 30%); lipid stain (top 30%)"
    )


def test_kits_and_factors_are_independent_columns():
    """Distinct sources, so a run with only kits gets no factors column."""
    main_df = selection_main_df("selection_kits", (["EasySep CD4"],))

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert KITS_COLUMN in sra_df.columns
    assert ENRICHMENT_FACTORS_COLUMN not in sra_df.columns


# suspension_type


def test_suspension_type_uses_the_sample_field():
    main_df = provider_main_df(tissues_suspension_type=["cell", "cell"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SUSPENSION_COLUMN] == "cell"


def test_suspension_type_fills_a_gap_like_a_required_column():
    """Being optional decides only whether the column exists, not how gaps fill."""
    main_df = provider_main_df(tissues_suspension_type=["cell", None])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, SUSPENSION_COLUMN] == "cell; not provided"


def test_no_suspension_type_column_omits_it():
    sra_df = make_flattener().create_sra_biosample_dataframe(provider_main_df())

    assert SUSPENSION_COLUMN not in sra_df.columns


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
    """donor_main_df plus age bounds: 32-35 for one donor, 29-30 for the other."""
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


def test_age_bounds_include_the_gap_in_the_set():
    main_df = age_main_df()
    main_df["tissues_lower_bound_age"] = main_df["human_donors_cxg_donor_id"].map(
        {889023040: 32, 889081306: None}
    )

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "pooled: 32 years, not provided"


def test_lower_bound_alone_omits_only_the_upper_column():
    main_df = age_main_df().drop(columns=["tissues_upper_bound_age"])

    sra_df = make_flattener().create_sra_biosample_dataframe(main_df)

    assert sra_df.loc[0, AGE_LOWER_COLUMN] == "pooled: 29 years, 32 years"
    assert AGE_UPPER_COLUMN not in sra_df.columns


# PROP_MAP_SRA_BIOSAMPLE


def test_prop_map_sends_both_library_types_to_sample_name():
    """The rename relies on both library keys sharing one output name."""
    droplet_key = "droplet_based_libraries_CRO_group_identifier"
    assert PROP_MAP_SRA_BIOSAMPLE[droplet_key] == "sample_name"
    assert PROP_MAP_SRA_BIOSAMPLE["plate_based_libraries_CRO_group_identifier"] == "sample_name"
    # The per-sample alias is deliberately absent: the key is the CRO group.
    assert "sample_alias" not in PROP_MAP_SRA_BIOSAMPLE
    # The library alias is absent too, now the group identifier keys the sheet.
    assert "droplet_based_libraries_aliases" not in PROP_MAP_SRA_BIOSAMPLE
    # sample_name is the one output column without the SRA required marker.
    assert not PROP_MAP_SRA_BIOSAMPLE[droplet_key].startswith("*")
    assert ORGANISM_COLUMN.startswith("*")


def test_prop_map_holds_only_renames():
    """The derived columns are named in the flattener, not mapped from MAIN."""
    assert "isolate" not in PROP_MAP_SRA_BIOSAMPLE
    assert "age" not in PROP_MAP_SRA_BIOSAMPLE
    assert ISOLATE_COLUMN not in PROP_MAP_SRA_BIOSAMPLE.values()
    assert AGE_COLUMN not in PROP_MAP_SRA_BIOSAMPLE.values()
    assert BARCODE_COLUMN not in PROP_MAP_SRA_BIOSAMPLE.values()
