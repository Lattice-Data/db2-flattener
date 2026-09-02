from dataclasses import dataclass
from typing import Any, TypeAlias

# Contains information about objects
# Their types, and what fields can be found in each

# typing for various objects used with JSON profile parsing
Hierarchy: TypeAlias = dict[str, dict[str, dict]]
JSONProfile: TypeAlias = dict[str, Any]
FieldTypes: TypeAlias = dict[str, dict[str, str]]
ObjectConfig: TypeAlias = dict[str, dict[str, Any]]


@dataclass
class Configs:
    """
    Use as a container to hold the parsed configs from profile schemas
    Structure:
    FIELD_TYPES: {
        {field}: {
            "type": {datatype value},
            "elements {optional}": {datatype of collection items},
        }
    }

    OBJECT_CONFIG: {
        {object url_prefix}: {
            "api_type": {object API Name},
            "fields": list[fields],
            "references": {
                "{field}": {object url_prefix}
            }
        }
    }

    """

    FIELD_TYPES: FieldTypes
    OBJECT_CONFIG: ObjectConfig


# Audit/provenance fields present on nearly every Lattice schema profile.
# Excluded from OBJECT_CONFIG so they don't get flattened into a column
# (and, in submitted_by's case, resolved as a reference) for every object type.
EXCLUDED_FIELDS = {"creation_timestamp", "submitted_by"}

# URL length limit for chunking (includes base URL overhead)
MAX_URL_LENGTH = 3800
# Base URL overhead for chunking calculations (base URL + field params + safety margin)
BASE_URL_OVERHEAD = 700

# Fields requested when gathering ControlledTerm objects.
# DB2_utils.split_term_cell derives the term id from '@id', so term_name is the
# only other field needed downstream. The profile has 18 fields; requesting the
# rest is dead weight on every request.
CONTROLLED_TERM_FIELDS = ["@id", "term_name"]

# Columns kept from a GeneticModification guide RNA tabular file, in output order.
# Only columns that exist in the source file are written.
GUIDE_METADATA_COLUMNS = [
    "guide_id",
    "guide_protospacer",
    "guide_role",
    "guide_PAM",
    "guide_target_gene_id",
    "guide_target_gene_name",
]

# Keys use _term_name suffix for columns produced by DB2_utils.split_controlled_term_columns
PROP_MAP_GEO = {
    "cell_lines_sample_terms_term_name": "**cell_line",
    "droplet_based_libraries_CRO_group_identifier": "*library name",
    "droplet_based_libraries_dbxrefs": "*SRA Experiment or Run",
    "droplet_based_libraries_library_cardinality": "single or paired-end",
    "droplet_based_libraries_library_construction_technology_term_name": "library_protocol",
    "genetic_modifications_strategy": "genetic_modifications_strategy",
    "human_donors_cxg_donor_id": "donor_ids",
    "human_donors_sex": "donor_sex",
    "human_donors_taxa": "*organism",
    "non_human_donors_cxg_donor_id": "donor_ids",
    "non_human_donors_sex": "donor_sex",
    "non_human_donors_taxa": "*organism",
    "organoids_sample_terms_term_name": "**tissue",
    "primary_cell_cultures_enriched_cell_types_term_name": "**cell_type",
    "raw_file_samples": "samples",
    "raw_matrix_file_alias": "processed data file",
    "sequence_file_sets_sequencing_platform": "*instrument model",
    "tissues_developmental_stages_term_name": "donor_dev_stage",
    "tissues_enriched_cell_types_term_name": "**cell_type",
    "tissues_sample_terms_term_name": "**tissue",
    "tissues_sources_title": "source",
    "tissues_selection_kits": "selection_kits",
    "tissues_selection_markers": "selection_markers",
    "tissues_selection_methods": "selection_methods",
}

PROP_MAP_BIOHUB = {
    "raw_file_samples": "sample_name",
    "non_human_donors_cxg_donor_id": "donor_id",
    "human_donors_cxg_donor_id": "donor_id",
    "non_human_donors_taxa": "organism",
    "human_donors_taxa": "organism",
    "non_human_donors_sex": "sex",
    "human_donors_sex": "sex",
    "human_donors_ethnicity_term_name": "self_reported_ethnicity",
    "human_donors_ethnicity_term_id": "self_reported_ethnicity_ontology_term_id",
    "droplet_based_libraries_library_construction_technology_term_name": "assay",
    "tissues_upper_bound_age": "tissues_upper_bound_age",
    "tissues_lower_bound_age": "tissues_lower_bound_age",
    "tissues_age_units": "tissues_age_units",
    "tissues_diseases_term_name": "disease",
    "tissues_enriched_cell_types_term_name": "suspension_enriched_cell_types",
    "tissues_sample_terms_term_name": "tissue",
    "tissues_developmental_stages_term_name": "development_stage",
    "tissues_multiplexing_barcodes": "sample_probe_barcode",
    "tissues_@type": "tissue_type",
    "tissues_selection_markers": "suspension_enrichment_factors",
    "tissues_selection_kits": "suspension_selection_kits",
    "tissues_suspension_type": "suspension_type",
    "tissues_preservation_method": "preservation_method",
    "treatments_ontological_term_term_id": "experimental_condition_ontology_term_id",
    "treatments_ontological_term_term_name": "experimental_condition",
    "experimental_conditions_condition": "experimental_condition",
    "experimental_conditions_text_value": "experimental_perturbation",
    "experimental_conditions_upper_bound_duration": "experimental_conditions_upper_bound_duration",
    "experimental_conditions_lower_bound_duration": "experimental_conditions_lower_bound_duration",
    "experimental_conditions_duration_units": "experimental_conditions_duration_units",
    "genetic_modifications_strategy": "genetic_perturbation_strategy",
    "sequence_file_sets_is_pilot_order": "is_pilot_data",
}

BIOHUB_SORT_ONTOLOGY_IDS = [
    "experimental_condition_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
]

TISSUE_TYPE_MAP = {
    "CellLine": "cell line",
    "Organoid": "organoid",
    "PrimaryCellCulture": "primary cell culture",
    "Tissue": "tissue",
}

GEO_LIBRARY_CARDINALITY_MAP = {
    "dual": "paired-end",
}

GEO_INSTRUMENT_MODEL_MAP = {
    "Ultima Genomics UG 100": "UG 100",
}

GEO_EXPERIMENTAL_CONDITION_COLS = [
    "experimental_conditions_condition",
    "experimental_conditions_text_value",
    "experimental_conditions_lower_bound_duration",
    "experimental_conditions_upper_bound_duration",
    "experimental_conditions_duration_units",
]

GEO_LIBRARY_STRATEGY_FEATURE_COL = "droplet_based_libraries_feature_types"
GEO_LIBRARY_STRATEGY_PLATE_FEATURE_COL = "plate_based_libraries_feature_types"
GEO_SUSPENSION_TYPE_COLS = [
    "organoids_suspension_type",
    "cell_lines_suspension_type",
    "primary_cell_cultures_suspension_type",
    "tissues_suspension_type",
]
GEO_LIBRARY_STRATEGY_SOURCE_COLS = [
    GEO_LIBRARY_STRATEGY_FEATURE_COL,
    GEO_LIBRARY_STRATEGY_PLATE_FEATURE_COL,
    *GEO_SUSPENSION_TYPE_COLS,
]
GEO_LIBRARY_STRATEGY_MAP = {
    ("Gene Expression", "nucleus"): "snRNA-seq",
    ("Gene Expression", "cell"): "scRNA-seq",
    ("ATAC", None): "scATAC-seq",
}

GEO_FLEX_LIBRARY_PROTOCOLS = {
    "10x gene expression flex v1",
    "10x gene expression flex",
    "10x Flex Apex",
    "10x GEM-X Flex v1",
    "10x Next GEM Flex v1",
}

GEO_TREATMENT_COLS = [
    "treatments_description",
    "treatments_lower_bound_duration",
    "treatments_duration_units",
    "treatments_ontological_term_term_name",
]
GEO_TITLE_TREATMENT_COLS = [
    "treatments_description",
    "treatments_lower_bound_duration",
    "treatments_upper_bound_duration",
    "treatments_duration_units",
]

GENETIC_PERTURBATION_MAP = {
    "activation screen": "CRISPR activation screen",
    "interference screen": "CRISPR interference screen",
    "knockout mutation": "CRISPR knockout mutant",
    "knockout screen": "CRISPR knockout screen",
}

REFORMAT_LIST = [
    "sample_probe_barcode",
    "suspension_enrichment_factors",
    "suspension_selection_kits",
]
