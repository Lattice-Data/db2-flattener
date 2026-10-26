from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import (
    GENETIC_PERTURBATION_MAP,
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
    assert PROP_MAP_SAMPLES["human_donors_taxa_term_name"] == "organism"
    assert PROP_MAP_SAMPLES["non_human_donors_taxa_term_name"] == "organism"
    assert "human_donors_ethnicity_term_name" not in PROP_MAP_SAMPLES
    assert PROP_MAP_SAMPLES["cell_lines_intended_cell_types_term_name"] == "**cell type"
    assert PROP_MAP_SAMPLES["cell_lines_sample_terms_term_name"] == "**cell_line"
    assert PROP_MAP_SAMPLES["tissues_sample_terms_term_name"] == "tissue"
    assert PROP_MAP_SAMPLES["tissues_enriched_cell_types_term_name"] == "**cell type"
    assert PROP_MAP_SAMPLES["primary_cell_cultures_enriched_cell_types_term_name"] == "**cell type"
    assert PROP_MAP_SAMPLES["tissues_multiplexing_barcodes"] == "sample_probe_barcode"
    assert PROP_MAP_SAMPLES["cell_lines_multiplexing_barcodes"] == "sample_probe_barcode"
    assert PROP_MAP_SAMPLES["tissues_RT_indexes"] == "RT_index"
    assert PROP_MAP_SAMPLES["cell_lines_RT_indexes"] == "RT_index"
    assert PROP_MAP_SAMPLES["organoids_RT_indexes"] == "RT_index"
    assert PROP_MAP_SAMPLES["primary_cell_cultures_RT_indexes"] == "RT_index"
    assert PROP_MAP_SAMPLES["tissues_selection_kits"] == "selection_kits"
    assert PROP_MAP_SAMPLES["tissues_selection_markers"] == "selection_markers"
    assert PROP_MAP_SAMPLES["tissues_selection_methods"] == "selection_methods"
    assert PROP_MAP_SAMPLES["tissues_developmental_stages_term_name"] == "donor_dev_stage"
    assert PROP_MAP_SAMPLES["tissues_sources_title"] == "source"
    assert PROP_MAP_SAMPLES["cell_lines_sources_title"] == "source"
    assert PROP_MAP_SAMPLES["experimental_conditions_summary"] == "experimental_condition_summary"
    assert (
        PROP_MAP_SAMPLES["experimental_conditions_controlled_terms_term_name"]
        == "experimental_condition_term"
    )
    assert PROP_MAP_SAMPLES["genetic_modifications_aliases"] == "genetic_modification_alias"
    assert PROP_MAP_SAMPLES["genetic_modifications_strategy"] == "genetic_modifications_strategy"
    assert not any(key.startswith("treatments_") for key in PROP_MAP_SAMPLES)


def test_create_samples_dataframe_renames_and_drops_unmapped():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1"],
            "sample_alias": ["s1"],
            "human_donors_cxg_donor_id": ["H1"],
            "human_donors_sex": ["female"],
            "human_donors_taxa_term_name": ["Homo sapiens"],
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
            "experimental_conditions_summary": ["LPS for 4 hours"],
            "genetic_modifications_aliases": ["lab:gm-screen-a"],
            "genetic_modifications_strategy": ["knockout screen"],
            "raw_file_samples": ["s1"],
            "tissues_@id": ["/tissues/s1/"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["processed data file"]) == ["rmf1"]
    assert list(result["pre_pooled_sample"]) == ["s1"]
    assert list(result["donor_id"]) == ["H1"]
    assert list(result["donor_sex"]) == ["female"]
    assert list(result["organism"]) == ["Homo sapiens"]
    assert list(result["donor_ethnicity"]) == ["European"]
    assert "human_donors_ethnicity_term_name" not in result.columns
    assert "tissues_@id" not in result.columns
    assert list(result["**cell type"]) == ["hepatocyte"]
    assert "enriched_cell_type" not in result.columns
    assert list(result["**cell_line"]) == ["HeLa"]
    assert list(result["tissue"]) == ["liver"]
    assert list(result["sample_probe_barcode"]) == ["BC001|CR001"]
    assert list(result["selection_kits"]) == ["EasySep"]
    assert list(result["selection_markers"]) == ["CD4"]
    assert list(result["selection_methods"]) == ["positive"]
    assert list(result["donor_dev_stage"]) == ["adult"]
    assert list(result["source"]) == ["Vendor X"]
    assert list(result["experimental_condition_summary"]) == ["LPS for 4 hours"]
    assert list(result["genetic_modification_alias"]) == ["lab:gm-screen-a"]
    assert "genetic_modifications_aliases" not in result.columns
    assert list(result["genetic_modifications_strategy"]) == ["CRISPR knockout screen"]
    assert "raw_file_samples" not in result.columns
    assert "tissues_@id" not in result.columns
    assert "sample_alias" not in result.columns


def test_create_samples_dataframe_maps_pcc_enriched_to_cell_type():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["p1"],
            "primary_cell_cultures_enriched_cell_types_term_name": ["T cell"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["**cell type"]) == ["T cell"]
    assert "primary_cell_cultures_enriched_cell_types_term_name" not in result.columns
    assert "enriched_cell_type" not in result.columns


def test_create_samples_dataframe_maps_cell_line_source_and_barcodes():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["cl1"],
            "cell_lines_sources_title": ["ATCC"],
            "cell_lines_multiplexing_barcodes": [["BC010", "CR010"]],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["source"]) == ["ATCC"]
    assert list(result["sample_probe_barcode"]) == ["BC010|CR010"]
    assert "cell_lines_sources_title" not in result.columns
    assert "cell_lines_multiplexing_barcodes" not in result.columns


def test_create_samples_dataframe_plate_library_names_barcode_column_hash_index():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["s1"],
            "tissues_multiplexing_barcodes": [["BC001", "CR001"]],
            "plate_based_libraries_@id": ["/plate_based_libraries/p1/"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["hash_index"]) == ["BC001|CR001"]
    assert "sample_probe_barcode" not in result.columns
    assert "plate_based_libraries_@id" not in result.columns


def test_attach_plate_library_id_does_not_duplicate_sample_rows():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1", "rmf1"],
            "plate_based_libraries_@id": [
                "/plate_based_libraries/p1/",
                "/plate_based_libraries/p2/",
            ],
        }
    )
    sample_df = pd.DataFrame(
        {
            "raw_matrix_file_alias": ["rmf1"],
            "sample_alias": ["s1"],
            "tissues_multiplexing_barcodes": [["BC001"]],
        }
    )

    attached = flattener._attach_plate_library_id(main_df, sample_df)
    result = flattener.create_samples_dataframe(attached)

    assert len(attached) == 1
    assert list(result["hash_index"]) == ["BC001"]
    assert "sample_probe_barcode" not in result.columns
    assert "plate_based_libraries_@id" not in result.columns


def test_create_samples_dataframe_plate_library_maps_rt_index():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["s1"],
            "tissues_RT_indexes": [["SCALEQUANT-A1", "SCALEQUANT-B1"]],
            "plate_based_libraries_@id": ["/plate_based_libraries/p1/"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["RT_index"]) == ["SCALEQUANT-A1,SCALEQUANT-B1"]
    assert "tissues_RT_indexes" not in result.columns


def test_create_samples_dataframe_plate_cell_line_maps_rt_index():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["cl1"],
            "cell_lines_RT_indexes": [["SCALEQUANT-C1"]],
            "plate_based_libraries_@id": ["/plate_based_libraries/p1/"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result["RT_index"]) == ["SCALEQUANT-C1"]
    assert "cell_lines_RT_indexes" not in result.columns


def test_create_samples_dataframe_droplet_omits_rt_index():
    flattener = make_flattener()
    sample_df = pd.DataFrame(
        {
            "sample_alias": ["s1"],
            "tissues_RT_indexes": [["SCALEQUANT-A1", "SCALEQUANT-B1"]],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert "RT_index" not in result.columns
    assert "tissues_RT_indexes" not in result.columns


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
            "human_donors_taxa_term_name": ["Homo sapiens", None],
            "non_human_donors_taxa_term_name": [None, "Mus musculus"],
        }
    )

    result = flattener.create_samples_dataframe(sample_df)

    assert list(result.columns).count("donor_id") == 1
    assert list(result["donor_id"]) == ["H1", "M1"]
    assert list(result["donor_sex"]) == ["female", "male"]
    assert list(result["organism"]) == ["Homo sapiens", "Mus musculus"]


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


def test_create_samples_dataframe_drops_treatment_and_built_condition_columns():
    flattener = make_flattener()
    result = flattener.create_samples_dataframe(
        pd.DataFrame(
            {
                "sample_alias": ["s1", "s2"],
                "treatments_ontological_term_term_name": ["LPS", None],
                "treatments_description": ["stimulation", None],
                "treatments_lower_bound_duration": [4, None],
                "experimental_conditions_summary": ["LPS for 4 hours", "control"],
                "experimental_conditions_controlled_terms_term_name": ["lipopolysaccharide", None],
                "experimental_conditions_controlled_terms_term_id": ["CHEBI:16412", None],
                "experimental_conditions_condition": ["treatment", "control"],
                "experimental_conditions_text_value": ["LPS", "PBS"],
            }
        )
    )

    assert list(result["pre_pooled_sample"]) == ["s1", "s2"]
    assert list(result["experimental_condition_summary"]) == ["LPS for 4 hours", "control"]
    assert result["experimental_condition_term"].iloc[0] == "lipopolysaccharide"
    assert pd.isna(result["experimental_condition_term"].iloc[1])
    for col in ("treatment", "treatment_description", "experimental_condition"):
        assert col not in result.columns
    assert not any(col.startswith(("treatments_", "experimental_conditions_")) for col in result)


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
