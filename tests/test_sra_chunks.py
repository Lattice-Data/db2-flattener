"""SRA file/biosample chunking by sample_name under a byte budget."""

from pathlib import Path

import pandas as pd

from db2_flattener.flatten.flattener import DB2Flattener, SRA_CHUNK_MAX_BYTES


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

    set1_files = pd.read_csv(sra_dir / "sra_file_set1.csv")
    set1_bio = pd.read_csv(sra_dir / "sra_biosample_set1.csv")
    set2_files = pd.read_csv(sra_dir / "sra_file_set2.csv")
    set2_bio = pd.read_csv(sra_dir / "sra_biosample_set2.csv")

    assert "file_size" not in set1_files.columns
    assert list(set1_files["sample_name"]) == ["A", "A"]
    assert list(set1_bio["sample_name"]) == ["A"]
    assert list(set2_files["sample_name"]) == ["B"]
    assert list(set2_bio["sample_name"]) == ["B"]


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


def test_non_crispr_sums_sequence_files_file_size():
    main_df = pd.DataFrame(
        {
            "droplet_based_libraries_CRO_group_identifier": ["LIB_A"],
            "droplet_based_libraries_aliases": [["lab:LIB_A_GEX"]],
            "droplet_based_libraries_@id": ["/droplet_based_libraries/gex/"],
            "droplet_based_libraries_feature_types": [["Gene Expression"]],
            "sequence_files_sra_file_names": ["R1.fastq.gz; R2.fastq.gz"],
            "sequence_files_file_size": ["100; 50"],
        }
    )
    sra_df = make_flattener().create_sra_files_dataframe(main_df)
    assert list(sra_df["file_size"]) == [150]


def test_sra_chunk_max_bytes_is_five_tb():
    assert SRA_CHUNK_MAX_BYTES == 5 * 1000**4
