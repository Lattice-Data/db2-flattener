from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import (
    GENETIC_PERTURBATION_MAP,
    GEO_EXPERIMENTAL_CONDITION_COLS,
    PROP_MAP_SAMPLES,
    Configs,
)


def make_flattener():
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = None
    flattener.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    flattener.gatherer = SimpleNamespace()
    return flattener


def test_prop_map_samples_names():
    assert PROP_MAP_SAMPLES["raw_matrix_file_alias"] == "processed data file"
    assert PROP_MAP_SAMPLES["sample_alias"] == "pre_pooled_sample"
    assert PROP_MAP_SAMPLES["human_donors_cxg_donor_id"] == "donor_id"
    assert PROP_MAP_SAMPLES["non_human_donors_cxg_donor_id"] == "donor_id"
    assert PROP_MAP_SAMPLES["human_donors_sex"] == "donor_sex"
    assert PROP_MAP_SAMPLES["non_human_donors_sex"] == "donor_sex"
    assert PROP_MAP_SAMPLES["human_donors_taxa"] == "organism"
    assert PROP_MAP_SAMPLES["non_human_donors_taxa"] == "organism"
    assert "human_donors_ethnicity_term_name" not in PROP_MAP_SAMPLES
    assert PROP_MAP_SAMPLES["cell_lines_intended_cell_types_term_name"] == "intended_cell_types"
    assert PROP_MAP_SAMPLES["cell_lines_sample_terms_term_name"] == "**cell_line"
    assert PROP_MAP_SAMPLES["tissues_sample_terms_term_name"] == "tissue"
    assert PROP_MAP_SAMPLES["tissues_enriched_cell_types_term_name"] == "enriched_cell_type"
    assert PROP_MAP_SAMPLES["tissues_multiplexing_barcodes"] == "sample_probe_barcode"
    assert PROP_MAP_SAMPLES["tissues_selection_kits"] == "selection_kits"
    assert PROP_MAP_SAMPLES["tissues_selection_markers"] == "selection_markers"
    assert PROP_MAP_SAMPLES["tissues_selection_methods"] == "selection_methods"
    assert PROP_MAP_SAMPLES["tissues_developmental_stages_term_name"] == "donor_dev_stage"
    assert PROP_MAP_SAMPLES["tissues_sources_title"] == "source"
    assert PROP_MAP_SAMPLES["treatments_ontological_term_term_name"] == "treatment"
    assert PROP_MAP_SAMPLES["treatments_description"] == "treatment_description"
    assert PROP_MAP_SAMPLES["genetic_modifications_strategy"] == "genetic_modifications_strategy"
    assert "treatments_lower_bound_duration" not in PROP_MAP_SAMPLES
    assert "treatments_upper_bound_duration" not in PROP_MAP_SAMPLES
    assert "treatments_duration_units" not in PROP_MAP_SAMPLES
    for col in GEO_EXPERIMENTAL_CONDITION_COLS:
        assert col not in PROP_MAP_SAMPLES


def test_create_samples_dataframe_renames_and_drops_unmapped():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1"],
            "sample_alias": ["s1"],
            "human_donors_cxg_donor_id": ["H1"],
            "human_donors_sex": ["female"],
            "human_donors_taxa": ["Homo sapiens"],
            "human_donors_ethnicity_term_name": ["European"],
            "cell_lines_intended_cell_types_term_name": ["hepatocyte"],
            "cell_lines_sample_terms_term_name": ["HeLa"],
            "tissues_sample_terms_term_name": ["liver"],
            "tissues_enriched_cell_types_term_name": ["hepatocyte"],
            "tissues_multiplexing_barcodes": [["BC001", "CR001"]],
            "tissues_selection_kits": ["EasySep"],
            "tissues_selection_markers": ["CD4"],
            "tissues_selection_methods": ["positive"],
            "tissues_developmental_stages_term_name": ["adult"],
            "tissues_sources_title": ["Vendor X"],
            "treatments_ontological_term_term_name": ["LPS"],
            "treatments_description": ["stimulation"],
            "genetic_modifications_strategy": ["knockout screen"],
            "raw_file_samples": ["s1"],
            "tissues_@id": ["/tissues/s1/"],
        }
    ).set_index("raw_matrix_file_alias")

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["processed data file"]) == ["rmf1"]
    assert list(result["pre_pooled_sample"]) == ["s1"]
    assert list(result["donor_id"]) == ["H1"]
    assert list(result["donor_sex"]) == ["female"]
    assert list(result["organism"]) == ["Homo sapiens"]
    assert list(result["donor_ethnicity"]) == ["European"]
    assert "human_donors_ethnicity_term_name" not in result.columns
    assert "tissues_@id" not in result.columns
    assert list(result["intended_cell_types"]) == ["hepatocyte"]
    assert list(result["**cell_line"]) == ["HeLa"]
    assert list(result["tissue"]) == ["liver"]
    assert list(result["enriched_cell_type"]) == ["hepatocyte"]
    assert list(result["sample_probe_barcode"]) == ["BC001|CR001"]
    assert list(result["selection_kits"]) == ["EasySep"]
    assert list(result["selection_markers"]) == ["CD4"]
    assert list(result["selection_methods"]) == ["positive"]
    assert list(result["donor_dev_stage"]) == ["adult"]
    assert list(result["source"]) == ["Vendor X"]
    assert list(result["treatment"]) == ["LPS"]
    assert list(result["treatment_description"]) == ["stimulation"]
    assert "treatment_duration" not in result.columns
    assert list(result["genetic_modifications_strategy"]) == ["CRISPR knockout screen"]
    assert "raw_file_samples" not in result.columns
    assert "tissues_@id" not in result.columns
    assert "sample_alias" not in result.columns


def test_create_samples_dataframe_maps_genetic_modifications_strategy():
    flattener = make_flattener()
    keys = list(GENETIC_PERTURBATION_MAP)
    sample_df = pd.DataFrame(
        {
            "sample_alias": [f"s{i}" for i in range(len(keys) + 1)],
            "genetic_modifications_strategy": [*keys, "custom edit"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["genetic_modifications_strategy"]) == [
        *GENETIC_PERTURBATION_MAP.values(),
        "custom edit",
    ]


def test_create_samples_dataframe_coalesces_human_and_non_human_donors():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1", "rmf2"],
            "sample_alias": ["human_s", "mouse_s"],
            "human_donors_cxg_donor_id": ["H1", None],
            "non_human_donors_cxg_donor_id": [None, "M1"],
            "human_donors_sex": ["female", None],
            "non_human_donors_sex": [None, "male"],
            "human_donors_taxa": ["Homo sapiens", None],
            "non_human_donors_taxa": [None, "Mus musculus"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result.columns).count("donor_id") == 1
    assert list(result["donor_id"]) == ["H1", "M1"]
    assert list(result["donor_sex"]) == ["female", "male"]
    assert list(result["organism"]) == ["Homo sapiens", "Mus musculus"]


def test_create_samples_dataframe_joins_description_and_duration():
    flattener = make_flattener()
    with_duration = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1", "s2"],
                "treatments_description": ["stimulation", None],
                "treatments_lower_bound_duration": [4, 2],
                "treatments_upper_bound_duration": [4, 8],
                "treatments_duration_units": ["hours", "hours"],
            }
        )
    )
    without_duration = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s3", "s4"],
                "treatments_description": ["washout", None],
            }
        )
    )

    assert list(with_duration["treatment_description"]) == [
        "stimulation; 4 hours",
        "2-8 hours",
    ]
    assert list(without_duration["treatment_description"]) == ["washout", "na"]
    assert "treatment_duration" not in with_duration.columns
    assert "treatments_lower_bound_duration" not in with_duration.columns
    assert "treatments_upper_bound_duration" not in with_duration.columns
    assert "treatments_duration_units" not in with_duration.columns


def test_create_samples_dataframe_fills_empty_treatment():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["s1", "s2", "s3"],
            "treatments_ontological_term_term_name": ["LPS", None, ""],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["treatment"]) == ["LPS", "no treatment", "no treatment"]


def test_create_samples_dataframe_strips_author_metadata_prefix():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["s1"],
            "tissues_author_metadata_mouse_litter_batch": ["A1"],
            "tissues_author_metadata_diet": ["fasted"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["mouse_litter_batch"]) == ["A1"]
    assert list(result["diet"]) == ["fasted"]
    assert "tissues_author_metadata_mouse_litter_batch" not in result.columns
    assert "tissues_author_metadata_diet" not in result.columns


def _assert_no_exp_source_cols(sample_df):
    for col in (*GEO_EXPERIMENTAL_CONDITION_COLS, "_exp_duration"):
        assert col not in sample_df.columns


def test_create_samples_dataframe_experimental_condition_equal_duration():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1"],
                "experimental_conditions_condition": ["treatment"],
                "experimental_conditions_text_value": ["LPS"],
                "experimental_conditions_lower_bound_duration": [4],
                "experimental_conditions_upper_bound_duration": [4],
                "experimental_conditions_duration_units": ["hours"],
            }
        )
    )

    assert list(result["experimental_condition"]) == ["treatment; LPS 4 hours"]
    _assert_no_exp_source_cols(result)


def test_create_samples_dataframe_experimental_condition_unequal_duration():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1"],
                "experimental_conditions_condition": ["treatment"],
                "experimental_conditions_text_value": ["LPS"],
                "experimental_conditions_lower_bound_duration": [2],
                "experimental_conditions_upper_bound_duration": [4],
                "experimental_conditions_duration_units": ["hours"],
            }
        )
    )

    assert list(result["experimental_condition"]) == ["treatment; LPS 2-4 hours"]
    _assert_no_exp_source_cols(result)


def test_create_samples_dataframe_experimental_condition_no_duration_columns():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1"],
                "experimental_conditions_condition": ["treatment"],
                "experimental_conditions_text_value": ["LPS"],
            }
        )
    )

    assert list(result["experimental_condition"]) == ["treatment; LPS"]
    _assert_no_exp_source_cols(result)


def test_create_samples_dataframe_experimental_condition_does_not_collapse_rows():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1", "s2"],
                "experimental_conditions_condition": ["treatment", "control"],
                "experimental_conditions_text_value": ["LPS", "PBS"],
            }
        )
    )

    assert list(result["pre_pooled_sample"]) == ["s1", "s2"]
    assert list(result["experimental_condition"]) == ["treatment; LPS", "control; PBS"]
    _assert_no_exp_source_cols(result)


def test_create_samples_dataframe_donor_ethnicity_sra_pooling():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1", "s2", "s3"],
                "tissues_@id": ["/tissues/t1/", "/tissues/t2/", "/tissues/t3/"],
                "human_donors_ethnicity_term_name": [
                    "European",
                    "Asian; European",
                    None,
                ],
            }
        )
    )

    assert list(result["donor_ethnicity"]) == [
        "European",
        "pooled: Asian, European",
        "not provided",
    ]


def test_create_samples_dataframe_donor_ethnicity_skips_non_tissue():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["tissue", "line", "org", "pcc"],
                "tissues_@id": ["/tissues/t1/", None, None, None],
                "cell_lines_@id": [None, "/cell_lines/c1/", None, None],
                "organoids_@id": [None, None, "/organoids/o1/", None],
                "primary_cell_cultures_@id": [None, None, None, "/primary_cell_cultures/p1/"],
                "human_donors_ethnicity_term_name": ["European"] * 4,
            }
        )
    )

    assert result["donor_ethnicity"].iloc[0] == "European"
    assert result["donor_ethnicity"].iloc[1:].isna().all()


def test_create_samples_dataframe_omits_donor_ethnicity_without_tissue():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["line"],
                "cell_lines_@id": ["/cell_lines/c1/"],
                "human_donors_ethnicity_term_name": ["European"],
            }
        )
    )

    assert "donor_ethnicity" not in result.columns
