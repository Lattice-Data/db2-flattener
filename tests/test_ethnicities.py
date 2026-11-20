"""HumanDonor.ethnicities (always a list) on BIOHUB, GEO, SRA biosample and SAMPLES.

Several ethnicities from exactly one donor read 'a; b'. With more than one
donor (or no donor ids to count), the merged set keeps 'pooled: a, b'.
BIOHUB takes sorted ontology ids joined with ' || '.
"""

import pandas as pd
import pytest
from test_create_dataframe import make_flattener as make_biohub_flattener
from test_geo_dataframe import gex_row
from test_geo_dataframe import make_flattener as make_geo_flattener
from test_samples_dataframe import make_flattener as make_samples_flattener
from test_sra_biosample import ETHNICITY_COLUMN, donor_main_df
from test_sra_biosample import make_flattener as make_sra_flattener

D1 = "/human_donors/d1/"
D2 = "/human_donors/d2/"
BOTH = ["Asian", "European"]


# --- BIOHUB ---


def biohub_frame(ethnicity_ids, organisms):
    n = len(ethnicity_ids)
    return pd.DataFrame(
        {
            "tissues_@type": [["Tissue", "Biosample", "Item"]] * n,
            "tissues_sample_terms_term_id": ["UBERON:0002048"] * n,
            "human_donors_cxg_donor_id": ["donor1"] * n,
            "human_donors_sex": ["female"] * n,
            "human_donors_taxa_term_id": organisms,
            "human_donors_ethnicities_term_id": ethnicity_ids,
        }
    )


def test_biohub_multiple_ethnicity_ids_are_sorted_and_joined():
    f = make_biohub_flattener()
    main_df = biohub_frame([["HANCESTRO:0008", "HANCESTRO:0005"]], ["NCBITaxon:9606"])

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["self_reported_ethnicity_ontology_term_id"]) == [
        "HANCESTRO:0005 || HANCESTRO:0008"
    ]


def test_biohub_single_ethnicity_id_is_unchanged():
    f = make_biohub_flattener()
    main_df = biohub_frame(["HANCESTRO:0005"], ["NCBITaxon:9606"])

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["self_reported_ethnicity_ontology_term_id"]) == ["HANCESTRO:0005"]


@pytest.mark.parametrize("empty", [[], (), None, ""])
def test_biohub_empty_ethnicities_default_to_unknown_for_human(empty):
    f = make_biohub_flattener()
    main_df = biohub_frame(["HANCESTRO:0005", empty], ["NCBITaxon:9606", "NCBITaxon:9606"])
    # a second row so the column exists with a real value; the empty one must default
    main_df["raw_file_samples"] = ["s1", "s2"]

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["self_reported_ethnicity_ontology_term_id"]) == [
        "HANCESTRO:0005",
        "unknown",
    ]


def test_biohub_empty_ethnicities_default_to_na_for_non_human():
    f = make_biohub_flattener()
    main_df = biohub_frame(["HANCESTRO:0005", []], ["NCBITaxon:9606", "NCBITaxon:10090"])
    main_df["raw_file_samples"] = ["s1", "s2"]

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["self_reported_ethnicity_ontology_term_id"]) == ["HANCESTRO:0005", "na"]


def test_biohub_ignores_the_removed_ethnicity_column():
    """The old single-valued column is gone from the schema; it is not read."""
    f = make_biohub_flattener()
    main_df = biohub_frame(["HANCESTRO:0005"], ["NCBITaxon:9606"]).drop(
        columns="human_donors_ethnicities_term_id"
    )
    main_df["human_donors_ethnicity_term_id"] = ["HANCESTRO:0005"]

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["self_reported_ethnicity_ontology_term_id"]) == ["unknown"]


# --- GEO ---


def geo_cell(*rows):
    main_df = pd.DataFrame(
        [
            gex_row(
                **{
                    "tissues_@id": f"/tissues/t{i}/",
                    "human_donors_ethnicities_term_name": terms,
                    **({} if donor is None else {"human_donors_@id": donor}),
                }
            )
            for i, (terms, donor) in enumerate(rows, start=1)
        ]
    ).dropna(axis=1, how="all")
    return list(make_geo_flattener().create_geo_dataframe(main_df)["donor_ethnicity"])


def test_geo_one_donor_with_several_ethnicities_has_no_pooled_prefix():
    assert geo_cell((BOTH, D1)) == ["Asian; European"]


def test_geo_one_donor_on_several_sample_rows_stays_one_donor():
    assert geo_cell((BOTH, D1), (BOTH, D1)) == ["Asian; European"]


def test_geo_two_donors_pool():
    assert geo_cell((["European"], D1), (["Asian"], D2)) == ["pooled: Asian, European"]


def test_geo_two_donors_one_with_none_adds_the_gap():
    assert geo_cell((BOTH, D1), (None, D2)) == ["pooled: Asian, European, not provided"]


def test_geo_one_donor_with_a_single_ethnicity_is_bare():
    assert geo_cell((["European"], D1)) == ["European"]


def test_geo_without_donor_ids_the_count_is_unknown_so_it_pools():
    assert geo_cell((BOTH, None)) == ["pooled: Asian, European"]


# --- SRA biosample ---


def sra_cell(ids, terms):
    main_df = donor_main_df(
        human_donors_ethnicities_term_name=terms,
        **({} if ids is None else {"human_donors_@id": ids}),
    )
    return make_sra_flattener().create_sra_biosample_dataframe(main_df).loc[0, ETHNICITY_COLUMN]


def test_sra_one_donor_with_several_ethnicities_has_no_pooled_prefix():
    assert sra_cell([D1] * 4, [BOTH] * 4) == "Asian; European"


def test_sra_two_donors_pool():
    assert sra_cell([D1, D2, D1, D2], [["European"], ["Asian"]] * 2) == "pooled: Asian, European"


def test_sra_two_donors_one_with_none_adds_the_gap():
    terms = [BOTH, None, BOTH, None]
    assert sra_cell([D1, D2, D1, D2], terms) == "pooled: Asian, European, not provided"


def test_sra_donors_joined_in_one_cell_count_as_two():
    """create_dataframe collapses a multi-donor sample into one '; '-joined @id."""
    assert sra_cell([f"{D1}; {D2}"] * 4, [BOTH] * 4) == "pooled: Asian, European"


def test_sra_without_donor_ids_the_count_is_unknown_so_it_pools():
    assert sra_cell(None, [BOTH] * 4) == "pooled: Asian, European"


# --- SAMPLES ---


def samples_cells(*rows):
    sample_df = pd.DataFrame(
        {
            "sample_alias": [f"s{i}" for i in range(len(rows))],
            "tissues_@id": [f"/tissues/t{i}/" for i in range(len(rows))],
            "human_donors_ethnicities_term_name": [terms for terms, _ in rows],
            **(
                {"human_donors_@id": [donor for _, donor in rows]}
                if all(donor is not None for _, donor in rows)
                else {}
            ),
        }
    )
    result = make_samples_flattener().create_samples_dataframe(sample_df)
    return list(result["donor_ethnicity"]), result


def test_samples_one_donor_with_several_ethnicities_has_no_pooled_prefix():
    cells, result = samples_cells((BOTH, D1))

    assert cells == ["Asian; European"]
    assert "human_donors_@id" not in result.columns
    assert "human_donors_ethnicities_term_name" not in result.columns


def test_samples_a_sample_pooling_two_donors_is_pooled():
    cells, _ = samples_cells((BOTH, f"{D1}; {D2}"))

    assert cells == ["pooled: Asian, European"]


def test_samples_per_row_donor_counts():
    cells, _ = samples_cells((BOTH, D1), (["European"], D2), (None, D2), ([], D1))

    assert cells == ["Asian; European", "European", "not provided", "not provided"]


def test_samples_without_donor_ids_the_count_is_unknown_so_it_pools():
    cells, _ = samples_cells((BOTH, None))

    assert cells == ["pooled: Asian, European"]
