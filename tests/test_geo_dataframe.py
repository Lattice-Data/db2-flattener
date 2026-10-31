from types import SimpleNamespace

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import (
    GEO_SUSPENSION_TYPE_COLS,
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
        "sequence_file_sets_run_cardinality": "paired-end",
        "raw_matrix_file_alias": "libA.h5",
        "raw_file_samples": "sample1",
        "tissues_sample_terms_term_name": "liver",
        "tissues_enriched_cell_types_term_name": "hepatocyte",
        "tissues_developmental_stages_term_name": "adult",
        "human_donors_cxg_donor_id": None,
        "human_donors_sex": None,
        "human_donors_taxa_term_name": None,
        "non_human_donors_cxg_donor_id": None,
        "non_human_donors_sex": None,
        "non_human_donors_taxa_term_name": None,
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
    assert PROP_MAP_GEO["human_donors_taxa_term_name"] == "*organism"
    assert PROP_MAP_GEO["raw_file_samples"] == "samples"
    assert PROP_MAP_GEO["tissues_developmental_stages_term_name"] == "donor_dev_stage"
    assert PROP_MAP_GEO["sequence_file_sets_sequencing_platform"] == "*instrument model"
    assert PROP_MAP_GEO["sequence_file_sets_run_cardinality"] == "*single or paired-end"
    assert PROP_MAP_GEO["genetic_modifications_aliases"] == "genetic_modification_alias"
    assert PROP_MAP_GEO["genetic_modifications_strategy"] == "genetic_modifications_strategy"
    assert PROP_MAP_GEO["tissues_selection_kits"] == "selection_kits"
    assert PROP_MAP_GEO["tissues_selection_markers"] == "selection_markers"
    assert PROP_MAP_GEO["tissues_selection_methods"] == "selection_methods"
    assert PROP_MAP_GEO["droplet_based_libraries_dbxrefs"] == "*SRA Experiment or Run"
    assert PROP_MAP_GEO["tissues_sources_title"] == "source"
    assert PROP_MAP_GEO["cell_lines_intended_cell_types_term_name"] == "**cell type"
    assert PROP_MAP_GEO["tissues_enriched_cell_types_term_name"] == "**cell type"
    assert PROP_MAP_GEO["primary_cell_cultures_enriched_cell_types_term_name"] == "**cell type"


def test_create_geo_dataframe_adds_new_columns():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                human_donors_cxg_donor_id="H1",
                human_donors_sex="female",
                human_donors_taxa_term_name="Homo sapiens",
                sequence_file_sets_sequencing_platform="Illumina NovaSeq 6000",
            )
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert "library_protocol" in geo_df.columns
    assert "*library strategy" in geo_df.columns
    assert "*title" in geo_df.columns
    assert list(geo_df["donor_ids"]) == ["H1"]
    assert list(geo_df["donor_sex"]) == ["female"]
    assert list(geo_df["*organism"]) == ["Homo sapiens"]
    assert list(geo_df["samples"]) == ["sample1"]
    assert list(geo_df["donor_dev_stage"]) == ["adult"]
    assert list(geo_df["**tissue"]) == ["liver"]
    assert list(geo_df["**cell type"]) == ["hepatocyte"]
    assert list(geo_df["*single or paired-end"]) == ["paired-end"]
    assert list(geo_df["*instrument model"]) == ["Illumina NovaSeq 6000"]
    assert "selection_kits" not in geo_df.columns
    assert "selection_markers" not in geo_df.columns
    assert "selection_methods" not in geo_df.columns


def test_create_geo_dataframe_donor_ethnicity_for_tissue():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                **{
                    "tissues_@id": "/tissues/t1/",
                    "human_donors_ethnicity_term_name": "European",
                }
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["donor_ethnicity"]) == ["European"]


def test_create_geo_dataframe_donor_ethnicity_pools_and_gaps():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                **{
                    "tissues_@id": "/tissues/t1/",
                    "human_donors_ethnicity_term_name": "European",
                }
            ),
            gex_row(
                **{
                    "tissues_@id": "/tissues/t2/",
                    "human_donors_ethnicity_term_name": None,
                }
            ),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["donor_ethnicity"]) == ["pooled: European, not provided"]


def test_create_geo_dataframe_omits_donor_ethnicity_for_cell_line():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                **{
                    "cell_lines_@id": "/cell_lines/c1/",
                    "human_donors_ethnicity_term_name": "European",
                }
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert "donor_ethnicity" not in geo_df.columns


def test_create_geo_dataframe_includes_intended_cell_types():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(cell_lines_intended_cell_types_term_name="hepatocyte")])

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["**cell type"]) == ["hepatocyte"]
    assert "intended_cell_types" not in geo_df.columns
    assert "cell_lines_intended_cell_types_term_name" not in geo_df.columns


def test_create_geo_dataframe_includes_tissue_selection_columns():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                tissues_selection_kits="EasySep",
                tissues_selection_markers="CD4",
                tissues_selection_methods="positive",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["selection_kits"]) == ["EasySep"]
    assert list(geo_df["selection_markers"]) == ["CD4"]
    assert list(geo_df["selection_methods"]) == ["positive"]
    assert "tissues_selection_kits" not in geo_df.columns
    assert "tissues_selection_markers" not in geo_df.columns
    assert "tissues_selection_methods" not in geo_df.columns
    assert list(geo_df.columns)[-1:] == ["processed data file"]


def test_create_geo_dataframe_collapses_human_and_non_human_donors():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libMouse",
                non_human_donors_cxg_donor_id="M1",
                non_human_donors_sex="male",
                non_human_donors_taxa_term_name="Mus musculus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libBoth",
                human_donors_cxg_donor_id="H1",
                human_donors_sex="female",
                human_donors_taxa_term_name="Homo sapiens",
                non_human_donors_cxg_donor_id="M2",
                non_human_donors_sex="male",
                non_human_donors_taxa_term_name="Mus musculus",
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
    assert both["donor_sex"] == "pooled male and female"
    assert both["*organism"] == ["Homo sapiens", "Mus musculus"]


def test_create_geo_dataframe_pools_mixed_donor_sex():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libMaleFemale",
                human_donors_sex="male",
                non_human_donors_sex="female",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libSameLib",
                human_donors_sex="female",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libSameLib",
                human_donors_sex="male",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libExtra",
                human_donors_sex="female",
                non_human_donors_sex="unknown",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libJoined",
                human_donors_sex="female; male",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libJoinedPlusMale",
                human_donors_sex="female; male",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libJoinedPlusMale",
                human_donors_sex="male",
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["libMaleFemale"]["donor_sex"] == "pooled male and female"
    assert geo_df.loc["libSameLib"]["donor_sex"] == "pooled male and female"
    assert geo_df.loc["libExtra"]["donor_sex"] == ["female", "unknown"]
    assert geo_df.loc["libJoined"]["donor_sex"] == "pooled male and female"
    assert geo_df.loc["libJoinedPlusMale"]["donor_sex"] == "pooled male and female"


def test_create_geo_dataframe_expands_raw_file_to_rightmost_columns():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(raw_matrix_file_alias="file1.h5"),
            gex_row(raw_matrix_file_alias="file2.h5"),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)
    raw_values = geo_df.loc[:, geo_df.columns == "processed data file"].iloc[0].tolist()

    assert list(geo_df.columns[-2:]) == ["processed data file", "processed data file"]
    assert list(geo_df.columns).count("processed data file") == 2
    assert raw_values == ["file1.h5", "file2.h5"]
    assert all(isinstance(v, str) for v in raw_values)


def test_create_geo_dataframe_maps_run_cardinality_single_end():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(sequence_file_sets_run_cardinality="single-end")]).dropna(
        axis=1, how="all"
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*single or paired-end"]) == ["single"]


def test_create_geo_dataframe_maps_run_cardinality_paired_end_with_index():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(sequence_file_sets_run_cardinality="paired-end-with-index")]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*single or paired-end"]) == ["paired-end"]


def test_create_geo_dataframe_maps_ultima_instrument_model():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(sequence_file_sets_sequencing_platform="Ultima Genomics UG 100")]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*instrument model"]) == ["UG 100"]


def test_create_geo_dataframe_extracts_sra_accession_from_dbxrefs():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(droplet_based_libraries_dbxrefs=["GEO:GSE1", "SRA:SRX123"])]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*SRA Experiment or Run"]) == ["SRX123"]
    assert "droplet_based_libraries_dbxrefs" not in geo_df.columns


def test_create_geo_dataframe_maps_tissues_sources_title_to_source():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(tissues_sources_title="Vendor X")]).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["source"]) == ["Vendor X"]
    assert "tissues_sources_title" not in geo_df.columns


def test_create_geo_dataframe_empty_when_dbxrefs_have_no_sra():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(droplet_based_libraries_dbxrefs=["GEO:GSE1"])]).dropna(
        axis=1, how="all"
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert len(geo_df) == 1
    assert pd.isna(geo_df["*SRA Experiment or Run"].iloc[0])
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]
    assert "droplet_based_libraries_dbxrefs" not in geo_df.columns


def test_create_geo_dataframe_extracts_biosample_from_dbxrefs():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_dbxrefs=[
                    "GEO:GSE1",
                    "SRA:SRX123",
                    "Biomaterial:SAMN53299868",
                ]
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*BioSample"]) == ["SAMN53299868"]
    # BioSample is read off the same cell, so the SRA accession must survive it.
    assert list(geo_df["*SRA Experiment or Run"]) == ["SRX123"]
    assert list(geo_df.columns).index("*BioSample") == (
        list(geo_df.columns).index("*SRA Experiment or Run") + 1
    )


def test_create_geo_dataframe_keeps_all_biosamples_from_dbxrefs():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_dbxrefs=[
                    "Biomaterial:SAMN1",
                    "Biomaterial:SAMN2",
                    "Biomaterial:SAMN1",
                ]
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert geo_df["*BioSample"].iloc[0] == ["SAMN1", "SAMN2"]


def test_create_geo_dataframe_empty_when_dbxrefs_have_no_biosample():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(droplet_based_libraries_dbxrefs=["GEO:GSE1", "SRA:SRX123"])]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert len(geo_df) == 1
    assert pd.isna(geo_df["*BioSample"].iloc[0])


def test_create_geo_dataframe_experimental_condition_summary():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                experimental_conditions_summary="LPS for 4 hours",
                experimental_conditions_condition="treatment",
                experimental_conditions_text_value="LPS",
                experimental_conditions_lower_bound_duration=4,
                experimental_conditions_upper_bound_duration=4,
                experimental_conditions_duration_units="hours",
            )
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["experimental_condition_summary"]) == ["LPS for 4 hours"]
    assert "experimental_condition" not in geo_df.columns
    for col in (
        "experimental_conditions_summary",
        "experimental_conditions_condition",
        "experimental_conditions_text_value",
        "experimental_conditions_lower_bound_duration",
        "experimental_conditions_upper_bound_duration",
        "experimental_conditions_duration_units",
    ):
        assert col not in geo_df.columns


def test_create_geo_dataframe_ignores_treatment_columns():
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
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert "treatment" not in geo_df.columns
    assert not any(col.startswith("treatments_") for col in geo_df.columns)


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

    assert list(geo_df["*library strategy"]) == ["scRNA-seq"]
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]
    _assert_no_library_strategy_source_cols(geo_df)


def test_create_geo_dataframe_library_strategy_snrna_seq():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(tissues_suspension_type="nucleus")]).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*library strategy"]) == ["snRNA-seq"]
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
    assert list(geo_df["*library strategy"]) == ["scATAC-seq"]
    assert list(geo_df["library_protocol"]) == ["10x 3' v3"]
    _assert_no_library_strategy_source_cols(geo_df)


def test_create_geo_dataframe_library_strategy_mixed_suspension_is_rna_seq():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libCellNucleus",
                tissues_suspension_type="cell; nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libNucleusCell",
                tissues_suspension_type="nucleus; cell",
            ),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["libCellNucleus"]["*library strategy"] == "RNA-Seq"
    assert geo_df.loc["libNucleusCell"]["*library strategy"] == "RNA-Seq"
    _assert_no_library_strategy_source_cols(geo_df)


def test_create_geo_dataframe_library_strategy_pools_cell_and_nucleus_rows():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libCellThenNucleus",
                tissues_suspension_type="cell",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libCellThenNucleus",
                tissues_suspension_type="nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libNucleusThenCell",
                tissues_suspension_type="nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libNucleusThenCell",
                tissues_suspension_type="cell",
            ),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["libCellThenNucleus"]["*library strategy"] == "RNA-Seq"
    assert geo_df.loc["libNucleusThenCell"]["*library strategy"] == "RNA-Seq"


def test_create_geo_dataframe_library_strategy_pools_rna_seq_with_single_cell_rows():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libMixed",
                tissues_suspension_type="cell; nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libMixed",
                tissues_suspension_type="cell",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libMixed",
                tissues_suspension_type="nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libWithAtac",
                tissues_suspension_type="cell; nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libWithAtac",
                tissues_suspension_type="cell",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libWithAtac",
                tissues_suspension_type="nucleus",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libWithAtac",
                droplet_based_libraries_feature_types="ATAC",
            ),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["libMixed"]["*library strategy"] == "RNA-Seq"
    atac_strategy = geo_df.loc["libWithAtac"]["*library strategy"]
    assert isinstance(atac_strategy, list)
    assert "scATAC-seq" in atac_strategy


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

    assert list(geo_df["*library strategy"]) == ["scRNA-seq"]
    _assert_no_library_strategy_source_cols(geo_df)


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


def test_create_geo_dataframe_molecule_from_plate_feature_types():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            {
                "plate_based_libraries_CRO_group_identifier": "plateGex",
                "plate_based_libraries_feature_types": "Gene Expression",
                "plate_based_libraries_library_construction_technology_term_name": "10x 3' v3",
            },
            {
                "plate_based_libraries_CRO_group_identifier": "plateFlex",
                "plate_based_libraries_feature_types": "Gene Expression",
                "plate_based_libraries_library_construction_technology_term_name": (
                    "10x gene expression flex v1"
                ),
            },
            {
                "plate_based_libraries_CRO_group_identifier": "plateAtac",
                "plate_based_libraries_feature_types": "ATAC",
                "plate_based_libraries_library_construction_technology_term_name": "10x ATAC",
            },
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["plateGex", "*molecule"] == "polyA RNA"
    assert geo_df.loc["plateFlex", "*molecule"] == "total RNA"
    assert geo_df.loc["plateAtac", "*molecule"] == "genomic DNA"


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
    assert list(geo_df.columns)[-1:] == ["processed data file"]


def test_create_geo_dataframe_maps_genetic_modifications_strategy():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(genetic_modifications_strategy="knockout screen"),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libJoined",
                genetic_modifications_strategy="activation screen; knockout screen",
            ),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["libA"]["genetic_modifications_strategy"] == "CRISPR knockout screen"
    assert (
        geo_df.loc["libJoined"]["genetic_modifications_strategy"]
        == "CRISPR activation screen; CRISPR knockout screen"
    )
    assert "genetic_perturbation_strategy" not in geo_df.columns


def test_create_geo_dataframe_maps_genetic_modification_alias():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(genetic_modifications_aliases="lab:gm-screen-a")]).dropna(
        axis=1, how="all"
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["genetic_modification_alias"]) == ["gm-screen-a"]
    assert "genetic_modifications_aliases" not in geo_df.columns


def test_create_geo_dataframe_title_includes_experimental_condition_summary():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                raw_file_samples=["sample1", "sample2"],
                genetic_modifications_strategy="knockout screen",
                experimental_conditions_summary="Liver cells stimulated with LPS",
                tissues_suspension_type="cell",
            )
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == [
        "libA scRNA-seq; Liver cells stimulated with LPS; pooled; CRISPR knockout screen"
    ]
    assert list(geo_df["experimental_condition_summary"]) == ["Liver cells stimulated with LPS"]
    assert "experimental_conditions_summary" not in geo_df.columns


def test_create_geo_dataframe_title_without_summary():
    flattener = make_flattener()
    main_df = pd.DataFrame([gex_row(tissues_suspension_type="cell")])

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == ["libA scRNA-seq"]


def test_create_geo_dataframe_title_pooled_and_strategy():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                raw_file_samples=["sample1", "sample2"],
                genetic_modifications_strategy="knockout screen",
                tissues_suspension_type="cell",
            )
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == ["libA scRNA-seq; pooled; CRISPR knockout screen"]


def test_create_geo_dataframe_title_maps_differing_strategies_in_library():
    """Per-sample matrix files carrying different strategies collapse to a list."""
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                raw_file_samples="sample1",
                genetic_modifications_strategy="knockout screen",
                tissues_suspension_type="cell",
            ),
            gex_row(
                raw_matrix_file_alias="libA-2.h5",
                raw_file_samples="sample2",
                genetic_modifications_strategy="activation screen",
                tissues_suspension_type="cell",
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == [
        "libA scRNA-seq; pooled; CRISPR activation screen; CRISPR knockout screen"
    ]


def test_create_geo_dataframe_title_pooled_from_joined_samples():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [gex_row(raw_file_samples="sample1; sample2", tissues_suspension_type="cell")]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == ["libA scRNA-seq; pooled"]


def test_create_geo_dataframe_title_joins_distinct_summaries_in_library():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(experimental_conditions_summary="zeta"),
            gex_row(raw_matrix_file_alias="libA-2.h5", experimental_conditions_summary="alpha"),
            gex_row(raw_matrix_file_alias="libA-3.h5", experimental_conditions_summary="zeta"),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == ["libA; alpha; zeta"]


def test_create_geo_dataframe_title_dedupes_within_joined_summaries():
    """A sample with two conditions arrives joined as 'a; b'."""
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(experimental_conditions_summary="a; b"),
            gex_row(raw_matrix_file_alias="libA-2.h5", experimental_conditions_summary="a"),
        ]
    ).dropna(axis=1, how="all")

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == ["libA; a; b"]
    assert list(geo_df["experimental_condition_summary"]) == ["a; b"]


def test_create_geo_dataframe_title_per_library():
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(experimental_conditions_summary="stimulated"),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libB",
                raw_matrix_file_alias="libB.h5",
                experimental_conditions_summary=None,
                tissues_suspension_type="cell",
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df).set_index("*library name")

    assert geo_df.loc["libA", "*title"] == "libA; stimulated"
    assert geo_df.loc["libB", "*title"] == "libB scRNA-seq"


def test_create_geo_dataframe_warns_on_title_over_120_characters(capsys):
    flattener = make_flattener()
    summary = "x" * 120
    main_df = pd.DataFrame([gex_row(experimental_conditions_summary=summary)])

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == [f"libA; {summary}"]
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "1 *title value(s) exceed 120 characters" in out


def test_create_geo_dataframe_no_warning_at_120_characters(capsys):
    flattener = make_flattener()
    summary = "x" * (120 - len("libA; "))
    main_df = pd.DataFrame([gex_row(experimental_conditions_summary=summary)])

    geo_df = flattener.create_geo_dataframe(main_df)

    assert len(geo_df["*title"].iloc[0]) == 120
    assert "exceed 120 characters" not in capsys.readouterr().out


def test_create_geo_dataframe_titles_stay_unique_for_a_shared_summary():
    """GEO requires unique sample titles; replicates often share a condition."""
    flattener = make_flattener()
    main_df = pd.DataFrame(
        [
            gex_row(
                experimental_conditions_summary="LPS for 4 hours",
                tissues_suspension_type="cell",
            ),
            gex_row(
                droplet_based_libraries_CRO_group_identifier="libB",
                raw_matrix_file_alias="libB.h5",
                experimental_conditions_summary="LPS for 4 hours",
                tissues_suspension_type="cell",
            ),
        ]
    )

    geo_df = flattener.create_geo_dataframe(main_df)

    assert list(geo_df["*title"]) == [
        "libA scRNA-seq; LPS for 4 hours",
        "libB scRNA-seq; LPS for 4 hours",
    ]
