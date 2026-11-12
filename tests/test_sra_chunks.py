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
    assert (sra_dir / "sra_sequence_fileset_set1.tsv").is_file()
    assert (sra_dir / "sra_sequence_fileset_set2.tsv").is_file()

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

    set1 = pd.read_csv(sra_dir / "sra_sequence_fileset_set1.tsv", sep="\t")
    set2 = pd.read_csv(sra_dir / "sra_sequence_fileset_set2.tsv", sep="\t")
    assert set1.values.tolist() == [["A_GEX", "u1"], ["A_GEX", "u2"], ["A_CRI.1", "u3"]]
    assert set2.values.tolist() == [["B_GEX", "u4"]]

    for n in (1, 2):
        files = pd.read_csv(sra_dir / f"sra_file_set{n}.csv")
        assert "_sequence_file_set_uuids" not in files.columns
        assert "sequence_file_set_uuid" not in files.columns


def _seq_file(name, size):
    return {
        "@id": f"/sequence_files/{name}/",
        "sra_file_names": [f"{name}.fastq.gz"],
        "file_size": size,
    }


def _seq_set(set_uuid, library_atid, *reads):
    file_set = {
        "@id": f"/sequence_file_sets/{set_uuid}/",
        "uuid": set_uuid,
        "library": library_atid,
    }
    for slot, name in zip(("read1", "read2"), reads, strict=False):
        file_set[slot] = f"/sequence_files/{name}/"
    return file_set


def test_sequence_fileset_sheet_from_gathered_data_to_tsv(tmp_path):
    """Gathered sets → SRA file frame → chunked CSVs, with no stubbed steps."""
    gex, cri = "/droplet_based_libraries/gex/", "/droplet_based_libraries/cri/"
    complete_data = {
        "libraries": {
            "gex-uuid": {
                "raw_matrix_files": [
                    {
                        "sequence_files": [_seq_file(n, 10) for n in ("g1", "g2", "g3")],
                        "sequence_file_sets": [
                            _seq_set("gex-set-b", gex, "g3"),
                            _seq_set("gex-set-a", gex, "g1", "g2"),
                        ],
                    }
                ]
            },
            "cri-uuid": {
                "raw_matrix_files": [
                    {
                        "sequence_files": [_seq_file(n, 5) for n in ("c1", "c2")],
                        "sequence_file_sets": [
                            _seq_set("cri-set-a", cri, "c1"),
                            _seq_set("cri-set-b", cri, "c2"),
                        ],
                    }
                ]
            },
        },
        "resolved_objects": {"SequenceFile": {}},
    }
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A", "LIB_A"],
            "droplet_based_libraries_@id": [gex, cri],
            "droplet_based_libraries_feature_types": [
                ["Gene Expression"],
                ["CRISPR Guide Capture"],
            ],
            "sequence_files_sra_file_names": ["g1.fastq.gz; g2.fastq.gz; g3.fastq.gz", None],
        }
    )

    flattener = make_flattener()
    sets_df = flattener.create_sequence_file_set_dataframe(complete_data)
    sra_files = flattener.create_sra_files_dataframe(main_df, sequence_file_set_df=sets_df)
    sra_dir = tmp_path / "SRA"
    sra_dir.mkdir()
    flattener._write_sra_chunked_outputs(
        sra_dir, sra_files, pd.DataFrame({"sample_name": ["LIB_A"]})
    )

    fileset = pd.read_csv(sra_dir / "sra_sequence_fileset_set1.tsv", sep="\t")
    assert list(fileset.columns) == ["library_ID", "sequence_file_set_uuid"]
    assert sorted(fileset.values.tolist()) == [
        ["LIB_A_CRI.1", "cri-set-a"],
        ["LIB_A_CRI.2", "cri-set-b"],
        ["LIB_A_GEX", "gex-set-a"],
        ["LIB_A_GEX", "gex-set-b"],
    ]
    files = pd.read_csv(sra_dir / "sra_file_set1.csv")
    assert set(fileset["library_ID"]) == set(files["library_ID"])
    assert not (sra_dir / "sra_sequence_fileset_set2.tsv").exists()


def _slot_check_data(sequence_files, file_set):
    return {
        "libraries": {
            "lib-uuid": {
                "raw_matrix_files": [
                    {"sequence_files": sequence_files, "sequence_file_sets": [file_set]}
                ]
            }
        },
        "resolved_objects": {"SequenceFile": {}},
    }


def _linked_file(name, set_id):
    return {**_seq_file(name, 1), "sequence_file_sets": [set_id]}


SLOT_SET = "/sequence_file_sets/s1/"


def test_slot_cross_check_is_quiet_when_slots_match_files(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1", "r2")
    data = _slot_check_data([_linked_file("r1", SLOT_SET), _linked_file("r2", SLOT_SET)], file_set)

    make_flattener().create_sequence_file_set_dataframe(data)

    assert "====WARNING====" not in capsys.readouterr().out


def test_slot_file_outside_derived_from_warns(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1", "r2")
    file_set["index1"] = {"@id": "/sequence_files/i1/", "aliases": ["lab:sample_I1.fastq.gz"]}
    data = _slot_check_data([_linked_file("r1", SLOT_SET), _linked_file("r2", SLOT_SET)], file_set)

    sfs_df = make_flattener().create_sequence_file_set_dataframe(data)

    out = capsys.readouterr().out
    assert "====WARNING==== 1 SequenceFileSet slot file(s) are not in any raw matrix" in out
    assert f"{SLOT_SET} index1: /sequence_files/i1/ (lab:sample_I1.fastq.gz)" in out
    assert "no slot pointing back" not in out
    # Filenames are unchanged: the unresolved I1 is still left out.
    assert sfs_df.loc[0, "filename"] == "r1.fastq.gz"
    assert "sample_I1.fastq.gz" not in sfs_df.iloc[0].tolist()


def test_file_listing_a_set_without_a_slot_warns(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1")
    data = _slot_check_data([_linked_file("r1", SLOT_SET), _linked_file("r2", SLOT_SET)], file_set)

    make_flattener().create_sequence_file_set_dataframe(data)

    out = capsys.readouterr().out
    assert "====WARNING==== 1 SequenceFile(s) list a SequenceFileSet" in out
    assert f"{SLOT_SET}: /sequence_files/r2/" in out
    assert "slot file(s)" not in out


def test_slot_mismatch_seen_on_two_raw_files_warns_once(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1")
    file_set["index1"] = "/sequence_files/i1/"
    raw_file = {"sequence_files": [_linked_file("r1", SLOT_SET)], "sequence_file_sets": [file_set]}
    data = {
        "libraries": {"lib-uuid": {"raw_matrix_files": [raw_file, dict(raw_file)]}},
        "resolved_objects": {"SequenceFile": {}},
    }

    make_flattener().create_sequence_file_set_dataframe(data)

    out = capsys.readouterr().out
    assert "====WARNING==== 1 SequenceFileSet slot file(s)" in out
    assert out.count(f"{SLOT_SET} index1: /sequence_files/i1/") == 1


def _gathered(raw_files, extra=()):
    """complete_data with resolved SequenceFiles as the gatherer leaves them."""
    resolved = {sf["@id"]: sf for rf in raw_files for sf in rf["sequence_files"]}
    resolved.update({sf["@id"]: sf for sf in extra})
    return {
        "libraries": {"lib-uuid": {"raw_matrix_files": raw_files}},
        "resolved_objects": {"SequenceFile": resolved},
    }


def test_set_split_across_raw_files_of_one_library_is_quiet(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1", "r2")
    data = _gathered(
        [
            {"sequence_files": [_linked_file("r1", SLOT_SET)], "sequence_file_sets": [file_set]},
            {"sequence_files": [_linked_file("r2", SLOT_SET)], "sequence_file_sets": [file_set]},
        ]
    )

    sfs_df = make_flattener().create_sequence_file_set_dataframe(data)

    assert "====WARNING====" not in capsys.readouterr().out
    assert sfs_df.loc[0, ["filename", "filename2"]].tolist() == ["r1.fastq.gz", "r2.fastq.gz"]


def test_gathered_slot_file_without_back_link_names_the_link(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1", "r2")
    r2 = {**_seq_file("r2", 1), "aliases": ["lab:r2"], "sequence_file_sets": []}
    data = _gathered(
        [{"sequence_files": [_linked_file("r1", SLOT_SET)], "sequence_file_sets": [file_set]}],
        extra=[r2],
    )

    make_flattener().create_sequence_file_set_dataframe(data)

    out = capsys.readouterr().out
    assert "1 SequenceFileSet slot file(s) are in derived_from but do not list that set" in out
    assert f"{SLOT_SET} read2: /sequence_files/r2/ (lab:r2)" in out
    assert "not in any raw matrix" not in out
    assert "reach no gathered library" in out
    assert "/sequence_files/r2/ (lab:r2): sequence_file_sets none" in out


def test_derived_from_file_reaching_no_library_warns(capsys):
    file_set = _seq_set("s1", "/droplet_based_libraries/gex/", "r1")
    orphan = {**_seq_file("orphan", 1), "sequence_file_sets": ["/sequence_file_sets/gone/"]}
    data = _gathered(
        [{"sequence_files": [_linked_file("r1", SLOT_SET)], "sequence_file_sets": [file_set]}],
        extra=[orphan],
    )

    make_flattener().create_sequence_file_set_dataframe(data)

    out = capsys.readouterr().out
    assert "====WARNING==== 1 derived_from SequenceFile(s) reach no gathered library" in out
    assert "/sequence_files/orphan/: sequence_file_sets /sequence_file_sets/gone/" in out
    assert "slot file(s)" not in out
