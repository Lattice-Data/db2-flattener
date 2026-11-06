"""SRA file/biosample chunking by sample_name under a byte budget."""

import pandas as pd

from db2_flattener.flatten.flattener import SRA_CHUNK_MAX_BYTES, DB2Flattener


def make_flattener():
    return DB2Flattener.__new__(DB2Flattener)


def test_chunk_keeps_sample_name_rows_together():
    sra_files = pd.DataFrame(
        {
            "sample_name": ["A", "A", "B"],
            "library_ID": ["A_GEX", "A_CRI", "B_GEX"],
            "file_size": [10, 20, 5],
        }
    )
    chunks = make_flattener()._chunk_sra_by_sample_name(sra_files, max_bytes=100)
    assert len(chunks) == 1
    assert list(chunks[0]["library_ID"]) == ["A_GEX", "A_CRI", "B_GEX"]


def test_chunk_splits_when_next_sample_would_exceed():
    sra_files = pd.DataFrame(
        {
            "sample_name": ["A", "A", "B", "B"],
            "library_ID": ["A_GEX", "A_CRI", "B_GEX", "B_CRI"],
            "file_size": [40, 40, 30, 30],
        }
    )
    chunks = make_flattener()._chunk_sra_by_sample_name(sra_files, max_bytes=100)
    assert len(chunks) == 2
    assert list(chunks[0]["sample_name"]) == ["A", "A"]
    assert list(chunks[1]["sample_name"]) == ["B", "B"]


def test_oversized_sample_alone_warns(capsys):
    sra_files = pd.DataFrame(
        {
            "sample_name": ["HUGE", "SMALL"],
            "library_ID": ["H", "S"],
            "file_size": [150, 10],
        }
    )
    chunks = make_flattener()._chunk_sra_by_sample_name(sra_files, max_bytes=100)
    assert len(chunks) == 2
    assert list(chunks[0]["sample_name"]) == ["HUGE"]
    assert list(chunks[1]["sample_name"]) == ["SMALL"]
    assert "writing alone" in capsys.readouterr().out


def test_write_sra_chunked_outputs_pairs_biosample(tmp_path, monkeypatch):
    flattener = make_flattener()
    monkeypatch.setattr(
        "db2_flattener.flatten.flattener.SRA_CHUNK_MAX_BYTES",
        100,
    )
    # Re-bind the method default by wrapping chunk helper used inside write
    real_chunk = flattener._chunk_sra_by_sample_name

    def chunk_with_budget(sra_files_df, max_bytes=100):
        return real_chunk(sra_files_df, max_bytes=max_bytes)

    flattener._chunk_sra_by_sample_name = chunk_with_budget

    sra_dir = tmp_path / "SRA"
    sra_dir.mkdir()
    sra_files = pd.DataFrame(
        {
            "sample_name": ["A", "A", "B"],
            "library_ID": ["A_GEX", "A_CRI", "B_GEX"],
            "file_size": [60, 60, 10],
        }
    )
    biosample = pd.DataFrame(
        {
            "sample_name": ["A", "B"],
            "organism": ["human", "mouse"],
        }
    )
    flattener._write_sra_chunked_outputs(sra_dir, sra_files, biosample)

    assert (sra_dir / "sra_file_set1.csv").is_file()
    assert (sra_dir / "sra_biosample_set1.csv").is_file()
    assert (sra_dir / "sra_file_set2.csv").is_file()
    assert (sra_dir / "sra_biosample_set2.csv").is_file()
    assert (sra_dir / "sra_sequence_fileset_set1.csv").is_file()
    assert (sra_dir / "sra_sequence_fileset_set2.csv").is_file()

    set1_files = pd.read_csv(sra_dir / "sra_file_set1.csv")
    set1_bio = pd.read_csv(sra_dir / "sra_biosample_set1.csv")
    set2_files = pd.read_csv(sra_dir / "sra_file_set2.csv")
    set2_bio = pd.read_csv(sra_dir / "sra_biosample_set2.csv")

    assert "file_size" not in set1_files.columns
    assert list(set1_files["sample_name"]) == ["A", "A"]
    assert list(set1_bio["sample_name"]) == ["A"]
    assert list(set2_files["sample_name"]) == ["B"]
    assert list(set2_bio["sample_name"]) == ["B"]


def test_write_sra_drops_empty_filename_columns(tmp_path):
    flattener = make_flattener()
    sra_dir = tmp_path / "SRA"
    sra_dir.mkdir()
    sra_files = pd.DataFrame(
        {
            "sample_name": ["A", "B"],
            "library_ID": ["A_GEX", "B_GEX"],
            "filename": ["a_R1.fastq.gz", "b_R1.fastq.gz"],
            "filename2": [pd.NA, pd.NA],
            "filename3": [None, ""],
            "file_size": [10, 10],
        }
    )
    biosample = pd.DataFrame({"sample_name": ["A", "B"]})
    flattener._write_sra_chunked_outputs(sra_dir, sra_files, biosample)

    written = pd.read_csv(sra_dir / "sra_file_set1.csv")
    assert "filename" in written.columns
    assert "filename2" not in written.columns
    assert "filename3" not in written.columns
    assert list(written["filename"]) == ["a_R1.fastq.gz", "b_R1.fastq.gz"]


def test_crispr_row_carries_sequence_file_set_file_size():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["lab:LIB_A_CRI"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/cri/"],
            "droplet_based_libraries_feature_types": [["CRISPR Guide Capture"]],
        }
    )
    sequence_file_set_df = pd.DataFrame(
        [
            {
                "sequence_file_set_@id": "/sequence_file_sets/a/",
                "library_@id": "/droplet_based_libraries/cri/",
                "filename": "setA_R1.fastq.gz",
                "file_size": 999,
            }
        ]
    )
    sra_df = make_flattener().create_sra_files_dataframe(
        main_df, sequence_file_set_df=sequence_file_set_df
    )
    assert list(sra_df["file_size"]) == [999]


def test_non_crispr_sums_sequence_file_set_file_sizes():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["lab:LIB_A_GEX"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
            "droplet_based_libraries_feature_types": [["Gene Expression"]],
            "sequence_files_sra_file_names": ["R1.fastq.gz; R2.fastq.gz"],
        }
    )
    sequence_file_set_df = pd.DataFrame(
        [
            {
                "sequence_file_set_@id": "/sequence_file_sets/a/",
                "library_@id": "/droplet_based_libraries/gex/",
                "filename": "R1.fastq.gz",
                "file_size": 100,
            },
            {
                "sequence_file_set_@id": "/sequence_file_sets/b/",
                "library_@id": "/droplet_based_libraries/gex/",
                "filename": "R2.fastq.gz",
                "file_size": 50,
            },
        ]
    )
    sra_df = make_flattener().create_sra_files_dataframe(
        main_df, sequence_file_set_df=sequence_file_set_df
    )
    assert list(sra_df["file_size"]) == [150]


def test_non_crispr_counts_equal_sequence_file_set_sizes():
    """Equal set totals must both count (MAIN _join_unique would have undercounted)."""
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["lab:LIB_A_GEX"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
            "droplet_based_libraries_feature_types": [["Gene Expression"]],
            "sequence_files_sra_file_names": ["R1.fastq.gz; R2.fastq.gz"],
        }
    )
    sequence_file_set_df = pd.DataFrame(
        [
            {
                "sequence_file_set_@id": "/sequence_file_sets/a/",
                "library_@id": "/droplet_based_libraries/gex/",
                "filename": "R1.fastq.gz",
                "file_size": 100,
            },
            {
                "sequence_file_set_@id": "/sequence_file_sets/b/",
                "library_@id": "/droplet_based_libraries/gex/",
                "filename": "R2.fastq.gz",
                "file_size": 100,
            },
        ]
    )
    sra_df = make_flattener().create_sra_files_dataframe(
        main_df, sequence_file_set_df=sequence_file_set_df
    )
    assert list(sra_df["file_size"]) == [200]


def test_sra_chunk_max_bytes_is_four_point_eight_tb():
    assert SRA_CHUNK_MAX_BYTES == 4_800_000_000_000


GEX_MAIN = {
    "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
    "droplet_based_libraries_aliases": [["lab:LIB_A_GEX"]],
    "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
    "droplet_based_libraries_feature_types": [["Gene Expression"]],
    "sequence_files_sra_file_names": ["R1.fastq.gz; R2.fastq.gz"],
}
CRI_MAIN = {
    "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
    "droplet_based_libraries_aliases": [["lab:LIB_A_CRI"]],
    "droplet_based_libraries_@id": ["/droplet_based_libraries/cri/"],
    "droplet_based_libraries_feature_types": [["CRISPR Guide Capture"]],
}


def _sets(library_atid, *uuids):
    return pd.DataFrame(
        [
            {
                "sequence_file_set_@id": f"/sequence_file_sets/{u}/",
                "sequence_file_set_uuid": u,
                "library_@id": library_atid,
                "filename": f"{u}_R1.fastq.gz",
                "file_size": 10,
            }
            for u in uuids
        ]
    )


def test_non_crispr_library_lists_every_sequence_file_set_uuid():
    sra_df = make_flattener().create_sra_files_dataframe(
        pd.DataFrame(GEX_MAIN),
        sequence_file_set_df=_sets("/droplet_based_libraries/gex/", "uuid-b", "uuid-a"),
    )
    assert list(sra_df["_sequence_file_set_uuids"]) == [["uuid-a", "uuid-b"]]

    rows = make_flattener()._sra_sequence_file_set_rows(sra_df)
    assert list(rows.columns) == ["library_ID", "sequence_file_set_uuid"]
    assert rows.values.tolist() == [["LIB_A_GEX", "uuid-a"], ["LIB_A_GEX", "uuid-b"]]


def test_crispr_suffixed_libraries_each_get_their_own_uuid():
    sra_df = make_flattener().create_sra_files_dataframe(
        pd.DataFrame(CRI_MAIN),
        sequence_file_set_df=_sets("/droplet_based_libraries/cri/", "uuid-a", "uuid-b"),
    )
    assert "sequence_file_set_uuid" not in sra_df.columns

    rows = make_flattener()._sra_sequence_file_set_rows(sra_df)
    assert rows.values.tolist() == [["LIB_A_CRI.1", "uuid-a"], ["LIB_A_CRI.2", "uuid-b"]]


def test_library_without_sequence_file_set_warns_and_is_left_out(capsys):
    sra_df = make_flattener().create_sra_files_dataframe(pd.DataFrame(GEX_MAIN))
    assert list(sra_df["_sequence_file_set_uuids"]) == [[]]

    rows = make_flattener()._sra_sequence_file_set_rows(sra_df)
    assert rows.empty
    assert list(rows.columns) == ["library_ID", "sequence_file_set_uuid"]
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "LIB_A_GEX" in out


def test_sequence_file_set_uuid_falls_back_to_atid():
    complete_data = {
        "libraries": {
            "lib-uuid": {
                "raw_matrix_files": [
                    {
                        "sequence_files": [],
                        "sequence_file_sets": [
                            {
                                "@id": "/sequence_file_sets/from-atid/",
                                "library": "/droplet_based_libraries/gex/",
                            },
                            {
                                "@id": "/sequence_file_sets/other/",
                                "uuid": "real-uuid",
                                "library": "/droplet_based_libraries/gex/",
                            },
                        ],
                    }
                ],
            }
        },
        "resolved_objects": {"SequenceFile": {}},
    }
    sfs_df = make_flattener().create_sequence_file_set_dataframe(complete_data)
    by_id = sfs_df.set_index("sequence_file_set_@id")["sequence_file_set_uuid"]
    assert by_id["/sequence_file_sets/from-atid/"] == "from-atid"
    assert by_id["/sequence_file_sets/other/"] == "real-uuid"


def test_write_sra_sequence_fileset_per_chunk(tmp_path):
    flattener = make_flattener()
    real_chunk = flattener._chunk_sra_by_sample_name
    flattener._chunk_sra_by_sample_name = lambda df: real_chunk(df, max_bytes=100)

    sra_dir = tmp_path / "SRA"
    sra_dir.mkdir()
    sra_files = pd.DataFrame(
        {
            "sample_name": ["A", "A", "B"],
            "library_ID": ["A_GEX", "A_CRI.1", "B_GEX"],
            "file_size": [60, 30, 50],
            "_sequence_file_set_uuids": [["u1", "u2"], ["u3"], ["u4"]],
        }
    )
    biosample = pd.DataFrame({"sample_name": ["A", "B"]})
    flattener._write_sra_chunked_outputs(sra_dir, sra_files, biosample)

    set1 = pd.read_csv(sra_dir / "sra_sequence_fileset_set1.csv")
    set2 = pd.read_csv(sra_dir / "sra_sequence_fileset_set2.csv")
    assert set1.values.tolist() == [["A_GEX", "u1"], ["A_GEX", "u2"], ["A_CRI.1", "u3"]]
    assert set2.values.tolist() == [["B_GEX", "u4"]]

    for n in (1, 2):
        files = pd.read_csv(sra_dir / f"sra_file_set{n}.csv")
        assert "_sequence_file_set_uuids" not in files.columns
        assert "sequence_file_set_uuid" not in files.columns
