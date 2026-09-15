from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import PROP_MAP_SAMPLES, Configs


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
    assert PROP_MAP_SAMPLES["human_donors_ethnicity_term_id"] == "self_reported_ethnicity"
    assert PROP_MAP_SAMPLES["tissues_sample_terms_term_name"] == "tissue"
    assert PROP_MAP_SAMPLES["tissues_enriched_cell_types_term_name"] == "enriched_cell_tye"
    assert PROP_MAP_SAMPLES["tissues_multiplexing_barcodes"] == "sample_probe_barcode"
    assert PROP_MAP_SAMPLES["tissues_selection_kits"] == "selection_kits"
    assert PROP_MAP_SAMPLES["tissues_developmental_stages_term_name"] == "donor_dev_stage"
    assert PROP_MAP_SAMPLES["treatments_ontological_term_term_name"] == "treatment"
    assert PROP_MAP_SAMPLES["treatments_description"] == "treatment_description"
    assert PROP_MAP_SAMPLES["genetic_modifications_strategy"] == "genetic_modifications_strategy"
    assert "treatments_lower_bound_duration" not in PROP_MAP_SAMPLES
    assert "treatments_upper_bound_duration" not in PROP_MAP_SAMPLES
    assert "treatments_duration_units" not in PROP_MAP_SAMPLES


def test_create_samples_dataframe_renames_and_drops_unmapped():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1"],
            "sample_alias": ["s1"],
            "human_donors_cxg_donor_id": ["H1"],
            "human_donors_sex": ["female"],
            "human_donors_taxa": ["Homo sapiens"],
            "human_donors_ethnicity_term_id": ["HANCESTRO:0005"],
            "tissues_sample_terms_term_name": ["liver"],
            "tissues_enriched_cell_types_term_name": ["hepatocyte"],
            "tissues_multiplexing_barcodes": [["BC001", "CR001"]],
            "tissues_selection_kits": ["EasySep"],
            "tissues_developmental_stages_term_name": ["adult"],
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
    assert list(result["self_reported_ethnicity"]) == ["HANCESTRO:0005"]
    assert list(result["tissue"]) == ["liver"]
    assert list(result["enriched_cell_tye"]) == ["hepatocyte"]
    assert list(result["sample_probe_barcode"]) == ["BC001|CR001"]
    assert list(result["selection_kits"]) == ["EasySep"]
    assert list(result["donor_dev_stage"]) == ["adult"]
    assert list(result["treatment"]) == ["LPS"]
    assert list(result["treatment_description"]) == ["stimulation"]
    assert list(result["genetic_modifications_strategy"]) == ["knockout screen"]
    assert "raw_file_samples" not in result.columns
    assert "tissues_@id" not in result.columns
    assert "sample_alias" not in result.columns


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


def test_create_samples_dataframe_combines_equal_and_unequal_duration():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["s1", "s2"],
            "treatments_lower_bound_duration": [4, 2],
            "treatments_upper_bound_duration": [4, 8],
            "treatments_duration_units": ["hours", "hours"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["treatment_duration"]) == ["4 hours", "2-8 hours"]
    assert "treatments_lower_bound_duration" not in result.columns
    assert "treatments_upper_bound_duration" not in result.columns
    assert "treatments_duration_units" not in result.columns


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
