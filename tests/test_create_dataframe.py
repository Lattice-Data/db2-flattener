import pandas as pd
import pytest

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.schema.constants import Configs

MIN_CONFIGS = Configs(
    FIELD_TYPES={},
    OBJECT_CONFIG={
        "droplet_based_libraries": {
            "api_type": "DropletBasedLibrary",
            "fields": ["@id", "feature_types", "CRO_group_identifier", "aliases"],
            "references": {},
        },
        "tissues": {
            "api_type": "Tissue",
            "fields": ["@id", "aliases", "author_metadata", "sources"],
            "references": {},
        },
    },
)


def make_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = MIN_CONFIGS
    f.gatherer = None
    return f


def _tissue(sample_id="/tissues/s1/", alias="lab:sample1", author_metadata=None, sources=None):
    return {
        "@id": sample_id,
        "@type": ["Tissue"],
        "aliases": [alias],
        "author_metadata": author_metadata,
        "sources": sources,
    }


def _lib(lib_id, feature_types, cro="RSJS_fast_1", alias=None):
    return {
        "@id": lib_id,
        "@type": ["DropletBasedLibrary"],
        "feature_types": feature_types,
        "CRO_group_identifier": cro,
        "aliases": [alias or f"lab:{feature_types[0].replace(' ', '_')}"],
    }


def _rmf(rmf_id="/raw_matrix_files/r1/", alias="lab:rmf1", sample_id="/tissues/s1/"):
    return {
        "@id": rmf_id,
        "aliases": [alias],
        "samples": [sample_id],
    }


def _complete_data(libs_with_rmfs):
    """libs_with_rmfs: list of (library, [raw_matrix_files], [samples])"""
    libraries = {}
    for i, (lib, rmfs, samples) in enumerate(libs_with_rmfs):
        libraries[f"uuid-{i}"] = {
            "library": lib,
            "raw_matrix_files": rmfs,
            "samples": samples,
        }
    return {"libraries": libraries, "resolved_objects": {"ControlledTerm": {}}}


# --- _row_is_geo_library ---


@pytest.mark.parametrize(
    "row,expected",
    [
        ({"droplet_based_libraries_feature_types": ["Gene Expression"]}, True),
        ({"droplet_based_libraries_feature_types": ["ATAC"]}, True),
        ({"droplet_based_libraries_feature_types": ["CRISPR Guide Capture"]}, False),
        ({"droplet_based_libraries_feature_types": "Gene Expression"}, True),
        ({"plate_based_libraries_feature_types": ["Gene Expression"]}, True),
        # missing FT: plate assumed GEX
        ({"plate_based_libraries_@id": "/plate_based_libraries/p1/"}, True),
        # missing FT: droplet assumed non-GEX
        ({"droplet_based_libraries_@id": "/droplet_based_libraries/d1/"}, False),
    ],
)
def test_row_is_geo_library(row, expected):
    f = make_flattener()
    assert f._row_is_geo_library(pd.Series(row)) is expected


# --- create_dataframe: one row per (RMF, library) ---


def test_create_dataframe_one_row_per_rmf_library():
    f = make_flattener()
    sample = _tissue()
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"], alias="lab:gex")
    cri = _lib("/droplet_based_libraries/cri/", ["CRISPR Guide Capture"], alias="lab:cri")

    main_df, sample_df = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
                (cri, [rmf], [sample]),
            ]
        )
    )

    assert len(main_df) == 2
    assert set(main_df["droplet_based_libraries_@id"]) == {
        "/droplet_based_libraries/gex/",
        "/droplet_based_libraries/cri/",
    }
    assert set(main_df["raw_matrix_file_alias"]) == {"rmf1"}
    ftypes = set()
    for ft in main_df["droplet_based_libraries_feature_types"]:
        ftypes.update(ft if isinstance(ft, list) else [ft])
    assert "Gene Expression" in ftypes
    assert "CRISPR Guide Capture" in ftypes


def test_create_dataframe_dedupes_same_library_twice_on_rmf():
    """Same lib attached via two library buckets should not double MAIN rows."""
    f = make_flattener()
    sample = _tissue()
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
                (gex, [rmf], [sample]),
            ]
        )
    )

    assert len(main_df) == 1
    assert main_df.iloc[0]["droplet_based_libraries_@id"] == "/droplet_based_libraries/gex/"


def test_sample_df_not_inflated_when_rmf_has_two_libraries():
    f = make_flattener()
    sample = _tissue()
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"], alias="lab:gex")
    cri = _lib("/droplet_based_libraries/cri/", ["CRISPR Guide Capture"], alias="lab:cri")

    main_df, sample_df = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
                (cri, [rmf], [sample]),
            ]
        )
    )

    assert len(main_df) == 2
    assert sample_df is not None
    assert sample_df["raw_matrix_file_alias"].nunique() == 1
    assert len(sample_df) == 1


# --- author_metadata explode ---


def test_author_metadata_dict_exploded_to_columns():
    f = make_flattener()
    sample = _tissue(author_metadata={"mouse litter batch": "A1", "diet": "fasted"})
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
            ]
        )
    )

    # The field name is retained in the column: create_biohub_dataframe selects
    # these columns with re.search('_author_metadata_', k), and
    # strip_author_metadata_column_prefix only drops the prefix later, for BIOHUB
    assert "tissues_author_metadata_mouse_litter_batch" in main_df.columns
    assert "tissues_author_metadata_diet" in main_df.columns
    assert main_df.iloc[0]["tissues_author_metadata_mouse_litter_batch"] == "A1"
    assert main_df.iloc[0]["tissues_author_metadata_diet"] == "fasted"
    assert "tissues_author_metadata" not in main_df.columns


def test_author_metadata_non_dict_kept_as_normal_field():
    f = make_flattener()
    sample = _tissue(author_metadata=None)
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
            ]
        )
    )

    assert "tissues_author_metadata" in main_df.columns
    assert "tissues_mouse_litter_batch" not in main_df.columns


def test_sources_title_extracted_from_embedded_dict():
    f = make_flattener()
    sample = _tissue(sources={"title": "Vendor X"})
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
            ]
        )
    )

    assert main_df.iloc[0]["tissues_sources_title"] == "Vendor X"
    assert main_df.iloc[0]["tissues_sources"] == {"title": "Vendor X"}


def test_sources_title_extracted_from_list_of_dicts():
    f = make_flattener()
    sample = _tissue(sources=[{"title": "Vendor X"}])
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
            ]
        )
    )

    assert main_df.iloc[0]["tissues_sources_title"] == "Vendor X"


def test_sources_without_title_does_not_add_title_column():
    f = make_flattener()
    sample = _tissue(sources={"@id": "/sources/s1/"})
    rmf = _rmf()
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
            ]
        )
    )

    assert "tissues_sources_title" not in main_df.columns
    assert main_df.iloc[0]["tissues_sources"] == {"@id": "/sources/s1/"}


# --- _join_unique: keep boolean False ---


def test_join_unique_keeps_false():
    f = make_flattener()
    assert f._join_unique([False]) == "False"
    assert f._join_unique([True]) == "True"
    assert f._join_unique([False, True]) == "False; True"


def test_join_unique_still_drops_empty_values():
    f = make_flattener()
    assert f._join_unique([None, "", "  "]) is None
    assert f._join_unique(["a", None, "", "b"]) == "a; b"


def test_create_dataframe_keeps_is_pilot_order_false():
    f = make_flattener()
    f.configs = Configs(
        FIELD_TYPES={},
        OBJECT_CONFIG={
            **MIN_CONFIGS.OBJECT_CONFIG,
            "sequence_file_sets": {
                "api_type": "SequenceFileSet",
                "fields": ["is_pilot_order"],
                "references": {},
            },
        },
    )
    sample = _tissue()
    rmf = _rmf()
    rmf["sequence_file_sets"] = [{"is_pilot_order": False}]
    gex = _lib("/droplet_based_libraries/gex/", ["Gene Expression"])

    main_df, _ = f.create_dataframe(
        _complete_data(
            [
                (gex, [rmf], [sample]),
            ]
        )
    )

    assert main_df.iloc[0]["sequence_file_sets_is_pilot_order"] == "False"


def test_biohub_tissue_type_from_tissues_cell_lines_or_both():
    f = make_flattener()
    main_df = pd.DataFrame(
        {
            "tissues_@type": [
                ["Tissue", "Biosample", "Item"],
                None,
                ["Tissue", "Biosample", "Item"],
            ],
            "cell_lines_@type": [
                None,
                ["CellLine", "Biosample", "Item"],
                ["CellLine", "Biosample", "Item"],
            ],
            "tissues_sample_terms_term_name": ["lung", None, "lung"],
            "cell_lines_sample_terms_term_id": [None, "CL:0000000", None],
            "tissues_developmental_stages_term_name": ["adult", None, "adult"],
            "tissues_multiplexing_barcodes": ["BC001", None, "BC003"],
            "cell_lines_multiplexing_barcodes": [None, "BC005", None],
            "tissues_suspension_type": ["cell", None, "cell"],
            "cell_lines_suspension_type": [None, "cell", None],
            "tissues_preservation_method": ["fresh", None, "fresh"],
            "cell_lines_preservation_method": [None, "frozen", None],
            "human_donors_cxg_donor_id": ["donor1", "donor2", "donor3"],
            "human_donors_sex": ["female", "male", "female"],
            "human_donors_ethnicity_term_name": ["European", "Asian", "European"],
            "human_donors_taxa": ["Homo sapiens", "Homo sapiens", "Homo sapiens"],
        }
    )

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["tissue_type"]) == ["tissue", "cell line", "tissue"]
    assert list(biohub_df["tissue"]) == ["lung", "CL:0000000", "lung"]
    assert list(biohub_df["development_stage"]) == ["adult", "na", "adult"]
    assert list(biohub_df["donor_id"]) == ["donor1", "na", "donor3"]
    assert list(biohub_df["sex"]) == ["female", "na", "female"]
    assert list(biohub_df["self_reported_ethnicity"]) == ["European", "na", "European"]
    assert list(biohub_df["sample_probe_barcode"]) == ["BC001", "BC005", "BC003"]
    assert list(biohub_df["suspension_type"]) == ["cell", "cell", "cell"]
    assert list(biohub_df["preservation_method"]) == ["fresh", "frozen", "fresh"]


def test_biohub_empty_development_stage_defaults_to_unknown():
    f = make_flattener()
    main_df = pd.DataFrame(
        {
            "tissues_@type": [
                ["Tissue", "Biosample", "Item"],
                ["Tissue", "Biosample", "Item"],
                ["Tissue", "Biosample", "Item"],
                ["Tissue", "Biosample", "Item"],
                None,
            ],
            "cell_lines_@type": [
                None,
                None,
                None,
                None,
                ["CellLine", "Biosample", "Item"],
            ],
            "tissues_sample_terms_term_name": ["lung", "lung", "lung", "lung", None],
            "cell_lines_sample_terms_term_id": [None, None, None, None, "CL:0000000"],
            "tissues_developmental_stages_term_name": ["", None, pd.NA, "adult", None],
            "human_donors_cxg_donor_id": ["donor1", "donor2", "donor3", "donor4", "donor5"],
            "human_donors_sex": ["female", "female", "female", "female", "male"],
            "human_donors_ethnicity_term_name": [
                "European",
                "European",
                "European",
                "European",
                "Asian",
            ],
            "human_donors_taxa": [
                "Homo sapiens",
                "Homo sapiens",
                "Homo sapiens",
                "Homo sapiens",
                "Homo sapiens",
            ],
        }
    )

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["tissue_type"]) == [
        "tissue",
        "tissue",
        "tissue",
        "tissue",
        "cell line",
    ]
    assert list(biohub_df["development_stage"]) == [
        "unknown",
        "unknown",
        "unknown",
        "adult",
        "na",
    ]


def test_biohub_unspecified_sex_defaults_to_unknown():
    f = make_flattener()
    main_df = pd.DataFrame(
        {
            "tissues_@type": [
                ["Tissue", "Biosample", "Item"],
                ["Tissue", "Biosample", "Item"],
                None,
            ],
            "cell_lines_@type": [
                None,
                None,
                ["CellLine", "Biosample", "Item"],
            ],
            "tissues_sample_terms_term_name": ["lung", "lung", None],
            "cell_lines_sample_terms_term_id": [None, None, "CL:0000000"],
            "tissues_developmental_stages_term_name": ["adult", "adult", None],
            "human_donors_cxg_donor_id": ["donor1", "donor2", "donor3"],
            "human_donors_sex": ["unspecified", "female", "unspecified"],
            "human_donors_ethnicity_term_name": ["European", "European", "Asian"],
            "human_donors_taxa": ["Homo sapiens", "Homo sapiens", "Homo sapiens"],
        }
    )

    biohub_df = f.create_biohub_dataframe(main_df)

    assert list(biohub_df["tissue_type"]) == ["tissue", "tissue", "cell line"]
    assert list(biohub_df["sex"]) == ["unknown", "female", "na"]


# --- multiplexed raw matrix files: one file carrying several samples ---

MULTIPLEX_CONFIGS = Configs(
    FIELD_TYPES={},
    OBJECT_CONFIG={
        "droplet_based_libraries": {
            "api_type": "DropletBasedLibrary",
            "fields": ["@id", "feature_types", "CRO_group_identifier", "aliases"],
            "references": {},
        },
        "tissues": {
            "api_type": "Tissue",
            "fields": ["@id", "aliases", "suspension_type", "selection_kits"],
            "references": {},
        },
    },
)


def _multiplex_flattener():
    f = DB2Flattener.__new__(DB2Flattener)
    f.connection = None
    f.configs = MULTIPLEX_CONFIGS
    f.gatherer = None
    return f


def _mx_tissue(sample_id, alias, suspension_type="cell", selection_kits=None):
    return {
        "@id": sample_id,
        "@type": ["Tissue"],
        "aliases": [alias],
        "suspension_type": suspension_type,
        "selection_kits": selection_kits,
    }


def _mx_run(tissues):
    """One library, one raw matrix file, carrying every tissue given."""
    rmf = {
        "@id": "/raw_matrix_files/r1/",
        "aliases": ["lab:rmf1"],
        "samples": [t["@id"] for t in tissues],
    }
    lib = _lib("/droplet_based_libraries/d1/", ["Gene Expression"])
    return _multiplex_flattener().create_dataframe(_complete_data([(lib, [rmf], tissues)]))


MX_TWO = [
    _mx_tissue("/tissues/s1/", "lab:H1", "cell", ["EasySep A"]),
    _mx_tissue("/tissues/s2/", "lab:H2", "nucleus", ["EasySep B"]),
]


def test_multiplexed_file_gives_samples_a_row_each():
    """The point of the fix: neither sample is lost to the joined-alias merge."""
    _, sample_df = _mx_run(MX_TWO)

    assert len(sample_df) == 2
    assert list(sample_df["sample_alias"]) == ["H1", "H2"]
    assert list(sample_df["raw_matrix_file_alias"]) == ["rmf1", "rmf1"]
    assert list(sample_df["tissues_suspension_type"]) == ["cell", "nucleus"]


def test_multiplexed_file_keeps_main_one_row_per_library():
    """MAIN must not be multiplied by the file's sample count."""
    main_df, _ = _mx_run(MX_TWO)

    assert len(main_df) == 1
    assert main_df.loc[0, "raw_file_samples"] == "H1; H2"


def test_multiplexed_file_joins_differing_sample_values_in_main():
    main_df, _ = _mx_run(MX_TWO)

    assert main_df.loc[0, "sample_alias"] == "H1; H2"
    assert main_df.loc[0, "tissues_suspension_type"] == "cell; nucleus"
    assert main_df.loc[0, "tissues_selection_kits"] == "EasySep A; EasySep B"


def test_multiplexed_file_does_not_repeat_a_shared_value():
    """Both samples agreeing is one value, not 'cell; cell'."""
    both_cell = [
        _mx_tissue("/tissues/s1/", "lab:H1", "cell"),
        _mx_tissue("/tissues/s2/", "lab:H2", "cell"),
    ]
    main_df, sample_df = _mx_run(both_cell)

    assert len(sample_df) == 2
    assert main_df.loc[0, "tissues_suspension_type"] == "cell"


def test_single_sample_file_passes_its_values_through_unchanged():
    """A list-valued field stays a list rather than being joined into a string."""
    main_df, sample_df = _mx_run([MX_TWO[0]])

    assert len(sample_df) == 1
    assert main_df.loc[0, "tissues_selection_kits"] == ["EasySep A"]
    assert main_df.loc[0, "tissues_suspension_type"] == "cell"


# --- _combine_sample_values ---


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        pytest.param(["cell"], "cell", id="one-value-passes-through"),
        pytest.param(["cell", "cell"], "cell", id="identical-values-collapse"),
        pytest.param(["cell", "nucleus"], "cell; nucleus", id="distinct-values-join"),
        pytest.param([["A"]], ["A"], id="lone-list-is-not-joined"),
        pytest.param([None, "cell"], "cell", id="empty-values-ignored"),
        pytest.param([None, None], None, id="all-empty"),
        pytest.param([], None, id="no-values"),
    ],
)
def test_combine_sample_values(values, expected):
    assert make_flattener()._combine_sample_values(values) == expected


def test_combine_sample_values_keeps_controlled_terms_as_dicts():
    """_join_unique would str() these, and the term columns could not be split."""
    blood = {"@id": "/controlled_terms/UBERON:1/", "term_name": "blood"}
    liver = {"@id": "/controlled_terms/UBERON:2/", "term_name": "liver"}

    assert make_flattener()._combine_sample_values([blood, liver]) == [blood, liver]


def test_combine_sample_values_keeps_dicts_that_are_not_terms():
    """An embedded 'sources' has no '@id', so it must not be deduped away."""
    one = {"title": "Vendor X"}
    two = {"title": "Vendor Y"}

    assert make_flattener()._combine_sample_values([one, two]) == [one, two]
