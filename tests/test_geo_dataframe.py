from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import (
    GEO_SUSPENSION_TYPE_COLS,
    GEO_TITLE_TREATMENT_COLS,
    GEO_TREATMENT_COLS,
    PROP_MAP_GEO,
    Configs,
)


def make_flattener():
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = None
    flattener.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    flattener.gatherer = SimpleNamespace()
    return flattener


def gex_row(**overrides):
    row = {
        "droplet_based_libraries_CRO_group_identifier": "libA",
        "droplet_based_libraries_feature_types": "Gene Expression",
        "droplet_based_libraries_library_construction_technology_term_name": "10x 3' v3",
        "droplet_based_libraries_library_cardinality": "paired",
        "raw_matrix_file_alias": "libA.h5",
        "raw_file_samples": "sample1",
        "tissues_sample_terms_term_name": "liver",
        "tissues_enriched_cell_types_term_name": "hepatocyte",
        "tissues_developmental_stages_term_name": "adult",
        "human_donors_cxg_donor_id": None,
        "human_donors_sex": None,
        "human_donors_taxa": None,
        "non_human_donors_cxg_donor_id": None,
        "non_human_donors_sex": None,
        "non_human_donors_taxa": None,
    }
    row.update(overrides)
    return row


def test_prop_map_geo_keeps_library_strategy():
    assert (
        PROP_MAP_GEO["droplet_based_libraries_library_construction_technology_term_name"]
        == "library_protocol"
    )
    assert "*library strategy" not in PROP_MAP_GEO.values()
    assert PROP_MAP_GEO["human_donors_cxg_donor_id"] == "donor_ids"
    assert PROP_MAP_GEO["non_human_donors_cxg_donor_id"] == "donor_ids"
    assert PROP_MAP_GEO["human_donors_sex"] == "donor_sex"
    assert PROP_MAP_GEO["non_human_donors_sex"] == "donor_sex"
    assert PROP_MAP_GEO["human_donors_taxa"] == "*organism"
    assert PROP_MAP_GEO["raw_file_samples"] == "samples"
    assert PROP_MAP_GEO["tissues_developmental_stages_term_name"] == "donor_dev_stage"
    assert PROP_MAP_GEO["sequence_file_sets_sequencing_platform"] == "*instrument model"
    assert PROP_MAP_GEO["genetic_modifications_strategy"] == "genetic_modifications_strategy"


def test_create_geo_dataframe_adds_new_columns():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                human_donors_cxg_donor_id="H1",
                human_donors_sex="female",
                human_donors_taxa="Homo sapiens",
                sequence_file_sets_sequencing_platform="Illumina NovaSeq 6000",
            )
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert "library_protocol" in geo_df.columns
    assert "*library strategy" not in geo_df.columns
    assert list(geo_df["donor_ids"]) == ["H1"]
    assert list(geo_df["donor_sex"]) == ["female"]
    assert list(geo_df["*organism"]) == ["Homo sapiens"]
    assert list(geo_df["samples"]) == ["sample1"]
    assert list(geo_df["donor_dev_stage"]) == ["adult"]
    assert list(geo_df["**tissue"]) == ["liver"]
    assert list(geo_df["**cell_type"]) == ["hepatocyte"]
    assert list(geo_df["single or paired-end"]) == ["paired"]
    assert list(geo_df["*instrument model"]) == ["Illumina NovaSeq 6000"]


def test_create_geo_dataframe_collapses_human_and_non_human_donors():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libMouse",
                non_human_donors_cxg_donor_id="M1",
                non_human_donors_sex="male",
                non_human_donors_taxa="Mus musculus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libBoth",
                human_donors_cxg_donor_id="H1",
                human_donors_sex="female",
                human_donors_taxa="Homo sapiens",
                non_human_donors_cxg_donor_id="M2",
                non_human_donors_sex="male",
                non_human_donors_taxa="Mus musculus",
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert list(geo_df.columns).count("donor_ids") == 1
    assert list(geo_df.columns).count("donor_sex") == 1
    assert list(geo_df.columns).count("*organism") == 1

    mouse = geo_df.loc["libMouse"]
    assert mouse["donor_ids"] == "M1"
    assert mouse["donor_sex"] == "male"
    assert mouse["*organism"] == "Mus musculus"

    both = geo_df.loc["libBoth"]
    assert both["donor_ids"] == ["H1", "M2"]
    assert both["donor_sex"] == ["female", "male"]
    assert both["*organism"] == ["Homo sapiens", "Mus musculus"]


def test_create_geo_dataframe_expands_raw_file_to_rightmost_columns():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(raw_matrix_file_alias="file1.h5"),
            gex_row(raw_matrix_file_alias="file2.h5"),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)
    raw_values = geo_df.loc[:, geo_df.columns == "raw_file"].iloc[0].tolist()

    assert list(geo_df.columns[-2:]) == ["raw_file", "raw_file"]
    assert list(geo_df.columns).count("raw_file") == 2
    assert raw_values == ["file1.h5", "file2.h5"]
    assert all(isinstance(v, str) for v in raw_values)


def test_create_geo_dataframe_maps_dual_cardinality_to_paired_end():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(droplet_based_libraries_library_cardinality="dual")]).dropna(
        axis=1, how="all"
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["single or paired-end"]) == ["paired-end"]


def _assert_no_exp_source_cols(geo_df):
    for col in (
        "experimental_conditions_condition",
        "experimental_conditions_text_value",
        "experimental_conditions_lower_bound_duration",
        "experimental_conditions_upper_bound_duration",
        "experimental_conditions_duration_units",
        "_exp_duration",
    ):
        assert col not in geo_df.columns


def test_create_geo_dataframe_experimental_condition_equal_duration():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                experimental_conditions_condition="treatment",
                experimental_conditions_text_value="LPS",
                experimental_conditions_lower_bound_duration=4,
                experimental_conditions_upper_bound_duration=4,
                experimental_conditions_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["experimental_condition"]) == ["treatment; LPS 4 hours"]
    _assert_no_exp_source_cols(geo_df)
    assert list(geo_df.columns)[-1:] == ["raw_file"]


def test_create_geo_dataframe_experimental_condition_unequal_duration():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                experimental_conditions_condition="treatment",
                experimental_conditions_text_value="LPS",
                experimental_conditions_lower_bound_duration=2,
                experimental_conditions_upper_bound_duration=4,
                experimental_conditions_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["experimental_condition"]) == ["treatment; LPS 2-4 hours"]
    _assert_no_exp_source_cols(geo_df)
    assert list(geo_df.columns)[-1:] == ["raw_file"]


def test_create_geo_dataframe_experimental_condition_no_duration_columns():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                experimental_conditions_condition="treatment",
                experimental_conditions_text_value="LPS",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["experimental_condition"]) == ["treatment; LPS"]
    _assert_no_exp_source_cols(geo_df)
    assert list(geo_df.columns)[-1:] == ["raw_file"]


def _assert_no_library_strategy_source_cols(geo_df):
    for col in (
        "droplet_based_libraries_feature_types",
        "plate_based_libraries_feature_types",
        *GEO_SUSPENSION_TYPE_COLS,
    ):
        assert col not in geo_df.columns


def test_create_geo_dataframe_library_strategy_scrna_seq():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(tissues_suspension_type="cell")]).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["library_strategy"]) == ["scRNA-seq"]
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]
    _assert_no_library_strategy_source_cols(geo_df)


def test_create_geo_dataframe_library_strategy_snrna_seq():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(tissues_suspension_type="nucleus")]).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["library_strategy"]) == ["snRNA-seq"]
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]
    _assert_no_library_strategy_source_cols(geo_df)


def test_create_geo_dataframe_library_strategy_scatac_seq_keeps_atac_row():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_feature_types="ATAC",
                tissues_suspension_type="cell",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert len(geo_df) == 1
    assert list(geo_df["library_strategy"]) == ["scATAC-seq"]
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]
    _assert_no_library_strategy_source_cols(geo_df)


def test_create_geo_dataframe_library_strategy_collapses_suspension_sources():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                organoids_suspension_type="cell",
                cell_lines_suspension_type="cell",
                primary_cell_cultures_suspension_type="cell",
                tissues_suspension_type="cell",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["library_strategy"]) == ["scRNA-seq"]
    _assert_no_library_strategy_source_cols(geo_df)


def _assert_no_treatment_source_cols(geo_df):
    for col in (*GEO_TREATMENT_COLS, *GEO_TITLE_TREATMENT_COLS, "_title_treatment"):
        assert col not in geo_df.columns


def test_create_geo_dataframe_treatment_all_fields():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                treatments_ontological_term_term_name="lipopolysaccharide",
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration=4,
                treatments_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["treatment"]) == ["lipopolysaccharide; LPS stimulation 4 hours"]
    _assert_no_treatment_source_cols(geo_df)
    assert list(geo_df.columns)[-1:] == ["raw_file"]


def test_create_geo_dataframe_treatment_no_duration():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                treatments_ontological_term_term_name="lipopolysaccharide",
                treatments_description="LPS stimulation",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["treatment"]) == ["lipopolysaccharide; LPS stimulation"]
    _assert_no_treatment_source_cols(geo_df)


def test_create_geo_dataframe_treatment_term_only():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(treatments_ontological_term_term_name="lipopolysaccharide")]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["treatment"]) == ["lipopolysaccharide"]
    _assert_no_treatment_source_cols(geo_df)


def test_create_geo_dataframe_treatment_description_and_duration_only():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration=4,
                treatments_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["treatment"]) == ["LPS stimulation 4 hours"]
    _assert_no_treatment_source_cols(geo_df)


def test_create_geo_dataframe_molecule_genomic_dna_for_atac():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(droplet_based_libraries_feature_types="ATAC")]).dropna(
        axis=1, how="all"
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*molecule"]) == ["genomic DNA"]


def test_create_geo_dataframe_molecule_total_rna_for_flex_gex():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_library_construction_technology_term_name=(
                    "10x gene expression flex v1"
                )
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*molecule"]) == ["total RNA"]
    assert list(geo_df["library_protocol"]) == ["10x gene expression flex v1"]


def test_create_geo_dataframe_molecule_polya_rna_for_other_gex():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row()]).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*molecule"]) == ["polyA RNA"]
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]


def test_create_geo_dataframe_includes_stripped_author_metadata():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_author_metadata_mouse_litter_batch="A1",
                tissues_author_metadata_diet="fasted",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["mouse_litter_batch"]) == ["A1"]
    assert list(geo_df["diet"]) == ["fasted"]
    assert "tissues_author_metadata_mouse_litter_batch" not in geo_df.columns
    assert "tissues_author_metadata_diet" not in geo_df.columns
    assert list(geo_df.columns)[-1:] == ["raw_file"]


def test_create_geo_dataframe_maps_genetic_modifications_strategy():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(genetic_modifications_strategy="knockout screen")]).dropna(
        axis=1, how="all"
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["genetic_modifications_strategy"]) == ["CRISPR knockout screen"]
    assert "genetic_perturbation_strategy" not in geo_df.columns


def test_create_geo_dataframe_title_full_set():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_suspension_type="cell",
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration=4,
                treatments_upper_bound_duration=4,
                treatments_duration_units="hours",
                treatments_schema_version="19",
                treatments_ontological_term_term_name="lipopolysaccharide",
                genetic_modifications_strategy="knockout screen",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["title"]) == [
        "libA scRNA-seq; LPS stimulation 4 hours; CRISPR knockout screen"
    ]
    assert all("lipopolysaccharide" not in title for title in geo_df["title"])
    assert list(geo_df["library_strategy"]) == ["scRNA-seq"]
    assert list(geo_df["genetic_modifications_strategy"]) == ["CRISPR knockout screen"]
    assert "treatments_schema_version" not in geo_df.columns
    assert "treatments_ontological_term_term_name" not in geo_df.columns
    _assert_no_treatment_source_cols(geo_df)
    assert list(geo_df.columns)[-1:] == ["raw_file"]


def test_create_geo_dataframe_title_unequal_duration():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_suspension_type="cell",
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration=2,
                treatments_upper_bound_duration=4,
                treatments_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["title"]) == ["libA scRNA-seq; LPS stimulation 2-4 hours"]
    _assert_no_treatment_source_cols(geo_df)


def test_create_geo_dataframe_title_skips_empty_fields():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(tissues_suspension_type="cell")]).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["title"]) == ["libA scRNA-seq"]
    assert all("no treatment" not in title for title in geo_df["title"])
    if "treatment" in geo_df.columns:
        assert list(geo_df["treatment"]) != ["no treatment"]


def test_create_geo_dataframe_title_pooled_when_samples_list():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_suspension_type="cell",
                raw_file_samples=["sample1", "sample2"],
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration=4,
                treatments_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["title"]) == ["libA scRNA-seq; pooled; LPS stimulation 4 hours"]


def test_create_geo_dataframe_title_pooled_when_joined_samples():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(tissues_suspension_type="cell", raw_file_samples="sample1; sample2")]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["title"]) == ["libA scRNA-seq; pooled"]
    assert all("no treatment" not in title for title in geo_df["title"])
    if "treatment" in geo_df.columns:
        assert list(geo_df["treatment"]) != ["no treatment"]


def test_create_geo_dataframe_treatment_stays_empty_when_no_row_has_treatment():
    flattener = make_flattener()
    empty_treatment = {
        "tissues_suspension_type": "cell",
        "treatments_description": None,
        "treatments_ontological_term_term_name": None,
        "treatments_lower_bound_duration": None,
        "treatments_duration_units": None,
    }
    main_df = pd.DataFrame(
        [
            gex_row(**empty_treatment),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libB",
                raw_matrix_file_alias="libB.h5",
                **empty_treatment,
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert "treatment" in geo_df.columns
    assert all(val is pd.NA or pd.isna(val) or val != "no treatment" for val in geo_df["treatment"])
    assert list(geo_df["title"]) == ["libA scRNA-seq", "libB scRNA-seq"]
    assert all("no treatment" not in title for title in geo_df["title"])


def test_create_geo_dataframe_treatment_no_treatment_when_other_library_has_treatment():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_suspension_type="cell",
                treatments_ontological_term_term_name="lipopolysaccharide",
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration="4",
                treatments_duration_units="hours",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libB",
                raw_matrix_file_alias="libB.h5",
                tissues_suspension_type="cell",
                treatments_description=None,
                treatments_ontological_term_term_name=None,
                treatments_lower_bound_duration=None,
                treatments_duration_units=None,
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)
    geo_df = geo_df.set_index("*library name")

    assert geo_df.loc["libA", "treatment"] == "lipopolysaccharide; LPS stimulation 4 hours"
    assert geo_df.loc["libA", "title"] == "libA scRNA-seq; LPS stimulation 4 hours"
    assert geo_df.loc["libB", "treatment"] == "no treatment"
    assert geo_df.loc["libB", "title"] == "libB scRNA-seq; no treatment"


def test_create_geo_dataframe_title_lists_mixed_treatments_same_library():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_suspension_type="cell",
                treatments_ontological_term_term_name="lipopolysaccharide",
                treatments_description="LPS stimulation",
                treatments_lower_bound_duration="4",
                treatments_duration_units="hours",
            ),
            gex_row(
                tissues_suspension_type="cell",
                raw_matrix_file_alias="libA-2.h5",
                treatments_description=None,
                treatments_ontological_term_term_name=None,
                treatments_lower_bound_duration=None,
                treatments_duration_units=None,
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert len(geo_df) == 1
    assert list(geo_df["treatment"]) == [
        ["lipopolysaccharide; LPS stimulation 4 hours", "no treatment"]
    ]
    assert list(geo_df["title"]) == [
        "libA scRNA-seq; ['LPS stimulation 4 hours', 'no treatment']"
    ]
    assert all("lipopolysaccharide" not in title for title in geo_df["title"])
