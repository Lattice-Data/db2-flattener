from types import SimpleNamespace

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.schema.constants import PROP_MAP_BIOHUB, PROP_MAP_GEO, Configs

GM_ID = "/genetic_modifications/gm1/"
TF_ID = "/tabular_files/tf1/"
SAMPLE_ID = "/tissues/s1/"
LIB_ID = "/droplet_based_libraries/lib1/"
RMF_ID = "/raw_matrix_files/rmf1/"
SEQ_ID = "/sequence_files/sf1/"
SFS_ID = "/sequence_file_sets/sfs1/"
LIB_UUID = "lib1"

GUIDE_FILE = {
    "@id": TF_ID,
    "s3_uri": "s3://bucket/guide.tsv",
    "file_format": "tsv",
    "aliases": ["lab:guide"],
}

CONFIGS = Configs(
    FIELD_TYPES={
        "genetic_modification": {"type": "string"},
        "guide_rna_files": {"type": "array", "elements": "string"},
    },
    OBJECT_CONFIG={
        "tissues": {
            "api_type": "Tissue",
            "fields": ["@id", "genetic_modification"],
            "references": {"genetic_modification": "genetic_modifications"},
        },
        "genetic_modifications": {
            "api_type": "GeneticModification",
            "fields": ["@id", "guide_rna_files"],
            "references": {"guide_rna_files": "tabular_files"},
        },
        "tabular_files": {
            "api_type": "TabularFile",
            "fields": ["@id", "aliases", "file_format", "s3_uri"],
            "references": {},
        },
        "droplet_based_libraries": {
            "api_type": "DropletBasedLibrary",
            "fields": ["@id"],
            "references": {},
        },
        "sequence_files": {"api_type": "SequenceFile", "fields": ["@id"], "references": {}},
        "sequence_file_sets": {"api_type": "SequenceFileSet", "fields": ["@id"], "references": {}},
        "raw_matrix_files": {"api_type": "RawMatrixFile", "fields": ["@id"], "references": {}},
    },
)


def make_flattener():
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = None
    flattener.configs = CONFIGS
    flattener.gatherer = SimpleNamespace(resolved_objects={})
    return flattener


def test_create_dataframe_flattens_tabular_file_fields():
    flattener = make_flattener()
    complete_data = {
        "libraries": {
            LIB_UUID: {
                "library": {"@id": LIB_ID, "uuid": LIB_UUID},
                "samples": [],
                "raw_matrix_files": [
                    {
                        "@id": RMF_ID,
                        "aliases": ["lab:rmf1.h5"],
                        "samples": [],
                        "sequence_files": [],
                        "sequence_file_sets": [],
                        "tabular_files": [GUIDE_FILE],
                    }
                ],
            }
        },
        "resolved_objects": {},
    }

    main_df, sample_df = flattener.create_dataframe(complete_data)

    assert list(main_df["tabular_files_s3_uri"]) == [GUIDE_FILE["s3_uri"]]
    assert list(main_df["tabular_files_file_format"]) == ["tsv"]
    assert list(main_df["tabular_files_@id"]) == [TF_ID]
    assert list(main_df["tabular_files_aliases"]) == ["lab:guide"]
    assert sample_df is None
    assert "tabular_files_s3_uri" not in PROP_MAP_BIOHUB
    assert "tabular_files_s3_uri" not in PROP_MAP_GEO
    assert "tabular_files_file_format" not in PROP_MAP_BIOHUB
    assert "tabular_files_file_format" not in PROP_MAP_GEO


def test_create_dataframe_joins_multiple_tabular_files():
    flattener = make_flattener()
    other = {
        "@id": "/tabular_files/tf2/",
        "s3_uri": "s3://bucket/other.tsv",
        "file_format": "tsv",
    }
    complete_data = {
        "libraries": {
            LIB_UUID: {
                "library": {"@id": LIB_ID, "uuid": LIB_UUID},
                "samples": [],
                "raw_matrix_files": [
                    {
                        "@id": RMF_ID,
                        "aliases": ["lab:rmf1.h5"],
                        "samples": [],
                        "tabular_files": [GUIDE_FILE, other],
                    }
                ],
            }
        },
        "resolved_objects": {},
    }

    main_df, _sample_df = flattener.create_dataframe(complete_data)

    assert main_df["tabular_files_s3_uri"].iloc[0] == (
        "s3://bucket/guide.tsv; s3://bucket/other.tsv"
    )
    assert main_df["tabular_files_@id"].iloc[0] == f"{TF_ID}; /tabular_files/tf2/"


def make_gatherer(fetched_by_type):
    gatherer = DB2Gatherer.__new__(DB2Gatherer)
    gatherer.connection = None
    gatherer.configs = CONFIGS
    gatherer.resolved_objects = {}
    calls = []

    def fake_fetch(obj_type, object_ids, filter_field="uuid", fields=None):
        calls.append((obj_type, list(object_ids)))
        return list(fetched_by_type.get(obj_type, []))

    gatherer.chunk_and_fetch = fake_fetch
    gatherer.fetch_calls = calls
    return gatherer


def test_resolve_references_fetches_tabular_file_from_guide_rna_files():
    gm = {"@id": GM_ID, "guide_rna_files": [TF_ID]}
    gatherer = make_gatherer({"GeneticModification": [gm], "TabularFile": [GUIDE_FILE]})
    sample = {"@id": SAMPLE_ID, "genetic_modification": GM_ID}

    gatherer.resolve_references_for_samples({SAMPLE_ID: sample})

    assert gatherer.resolved_objects["GeneticModification"][GM_ID] == gm
    assert gatherer.resolved_objects["TabularFile"][TF_ID] == GUIDE_FILE
    assert ("TabularFile", ["tf1"]) in gatherer.fetch_calls


def test_add_references_attaches_tabular_files_to_library():
    gm = {"@id": GM_ID, "guide_rna_files": [TF_ID]}
    gatherer = make_gatherer({})
    gatherer.resolved_objects = {
        "GeneticModification": {GM_ID: gm},
        "TabularFile": {TF_ID: GUIDE_FILE},
    }
    library_data = {}
    sample = {"@id": SAMPLE_ID, "genetic_modification": GM_ID}

    gatherer.add_references_to_library(library_data, [sample])

    assert library_data["genetic_modifications"] == [gm]
    assert library_data["tabular_files"] == [GUIDE_FILE]


def test_map_raw_matrix_files_copies_tabular_files_onto_raw_file():
    gatherer = make_gatherer({})
    seq_file = {"@id": SEQ_ID, "sequence_file_sets": [SFS_ID]}
    file_set = {"@id": SFS_ID, "library": LIB_ID}
    gatherer.resolved_objects = {
        "SequenceFile": {SEQ_ID: seq_file},
        "SequenceFileSet": {SFS_ID: file_set},
    }
    libraries_data = {
        LIB_UUID: {
            "library": {"@id": LIB_ID, "uuid": LIB_UUID},
            "samples": [],
            "raw_matrix_files": [],
            "tabular_files": [GUIDE_FILE],
        }
    }
    raw_file = {"@id": RMF_ID, "derived_from": [SEQ_ID]}

    gatherer._map_raw_matrix_files_to_libraries([raw_file], libraries_data)

    attached = libraries_data[LIB_UUID]["raw_matrix_files"]
    assert len(attached) == 1
    assert attached[0]["tabular_files"] == [GUIDE_FILE]
