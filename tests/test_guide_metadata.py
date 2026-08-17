from types import SimpleNamespace
from urllib.parse import urljoin

import pandas as pd
import pytest

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.gather.lattice import _separator_for_file, download_file, read_tabular_file
from db2_flattener.schema.constants import Configs

GUIDE_FILE = {
    "@id": "/tabular_files/e0763e74-9046-4547-931f-f8788ec203fc/",
    "s3_uri": "s3://submissions-lattice/BCP/tabular_files/guide_template.tsv",
    "file_format": "tsv",
}

GUIDE_ROWS = pd.DataFrame(
    {
        "guide_id": ["g1", "g2"],
        "guide_protospacer": ["ACGT", "TGCA"],
        "guide_role": ["targeting", "non-targeting"],
        "guide_PAM": ["NGG", "NGG"],
        "guide_target_gene_id": ["ENSG1", ""],
        "guide_target_gene_name": ["GENE1", ""],
        "extra_col": ["drop-me", "drop-me"],
    }
)


def make_flattener(monkeypatch, read_df=None):
    flattener = DB2Flattener.__new__(DB2Flattener)
    flattener.connection = SimpleNamespace(server="https://example.com/", auth=("k", "s"))
    flattener.configs = Configs(FIELD_TYPES={}, OBJECT_CONFIG={})
    calls = []

    def fake_read(file_info, connection):
        calls.append(file_info)
        return read_df.copy() if read_df is not None else pd.DataFrame()

    monkeypatch.setattr("db2_flattener.flatten.flattener.DB2lattice.read_tabular_file", fake_read)
    flattener.gatherer = SimpleNamespace()
    return flattener, calls


def complete_data_with_gms(*gms):
    return {
        "resolved_objects": {
            "GeneticModification": {gm["@id"]: gm for gm in gms},
        }
    }


def guide_df_from_data(flattener, data):
    return flattener.create_guide_metadata_dataframe(flattener._resolve_guide_rna_file(data))


def test_no_genetic_modifications_returns_none(monkeypatch):
    flattener, calls = make_flattener(monkeypatch)
    assert flattener._resolve_guide_rna_file({"libraries": {}}) is None
    assert flattener._resolve_guide_rna_file({"resolved_objects": {}}) is None
    assert flattener.create_guide_metadata_dataframe(None) is None
    assert calls == []


def test_gm_without_guide_rna_files_returns_none(monkeypatch):
    flattener, calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "strategy": "interference screen"}
    )
    assert guide_df_from_data(flattener, data) is None
    assert calls == []


def test_one_unique_file_subsets_columns(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {
            "@id": "/genetic_modifications/aaa/",
            "guide_rna_files": GUIDE_FILE,
        }
    )
    result = guide_df_from_data(flattener, data)
    assert list(result.columns) == [
        "guide_id",
        "guide_protospacer",
        "guide_role",
        "guide_PAM",
        "guide_target_gene_id",
        "guide_target_gene_name",
    ]
    assert list(result["guide_id"]) == ["g1", "g2"]
    assert len(calls) == 1


def test_same_file_on_two_gms_is_read_once(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": [GUIDE_FILE]},
        {"@id": "/genetic_modifications/bbb/", "guide_rna_files": [GUIDE_FILE]},
    )
    result = guide_df_from_data(flattener, data)
    assert result is not None
    assert len(calls) == 1


def test_two_distinct_files_warns_and_skips(monkeypatch, capsys):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    other = {
        "@id": "/tabular_files/22222222-2222-2222-2222-222222222222/",
        "s3_uri": "s3://submissions-lattice/BCP/tabular_files/other.tsv",
        "file_format": "tsv",
    }
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": GUIDE_FILE},
        {"@id": "/genetic_modifications/bbb/", "guide_rna_files": other},
    )
    assert guide_df_from_data(flattener, data) is None
    assert calls == []
    captured = capsys.readouterr()
    assert "found 2 unique guide RNA TabularFiles" in captured.out
    assert GUIDE_FILE["@id"] in captured.out
    assert other["@id"] in captured.out


def test_guide_rna_files_as_id_strings_use_gathered_tabular_file(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {
            "@id": "/genetic_modifications/aaa/",
            "guide_rna_files": [GUIDE_FILE["@id"]],
        }
    )
    data["resolved_objects"]["TabularFile"] = {GUIDE_FILE["@id"]: GUIDE_FILE}
    result = guide_df_from_data(flattener, data)
    assert result is not None
    assert calls == [GUIDE_FILE]


def test_missing_expected_columns_returns_none(monkeypatch, capsys):
    flattener, _calls = make_flattener(monkeypatch, read_df=pd.DataFrame({"unrelated": [1]}))
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": GUIDE_FILE}
    )
    assert guide_df_from_data(flattener, data) is None
    assert "none of the expected GUIDE_METADATA columns" in capsys.readouterr().out


def test_keeps_only_present_guide_columns(monkeypatch):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame({"guide_id": ["g1"], "guide_role": ["targeting"]}),
    )
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": GUIDE_FILE}
    )
    result = guide_df_from_data(flattener, data)
    assert list(result.columns) == ["guide_id", "guide_role"]


@pytest.mark.parametrize(
    "file_info, expected",
    [
        ({"file_format": "tsv"}, "\t"),
        ({"file_format": "csv"}, ","),
        ({"s3_uri": "s3://bucket/file.tsv"}, "\t"),
        ({"s3_uri": "s3://bucket/file.csv"}, ","),
        ({}, None),
        ({"@id": "/tabular_files/abc/"}, None),
    ],
)
def test_separator_for_file(file_info, expected):
    assert _separator_for_file(file_info) == expected


def test_download_file_uses_lattice_download_url(monkeypatch):
    captured = {}

    class FakeResponse:
        content = b"guide_id\ng1\n"

        def raise_for_status(self):
            return None

    def fake_get(url, auth=None):
        captured["url"] = url
        captured["auth"] = auth
        return FakeResponse()

    monkeypatch.setattr("db2_flattener.gather.lattice.requests.get", fake_get)
    connection = SimpleNamespace(server="https://example.com/", auth=("k", "s"))
    content = download_file("/tabular_files/abc/", connection)
    assert content == b"guide_id\ng1\n"
    assert captured["url"] == urljoin("https://example.com/", "/tabular_files/abc/") + "@@download"
    assert captured["auth"] == ("k", "s")


def test_read_tabular_file_uses_s3_uri(monkeypatch):
    def fake_read_csv(path, sep=None):
        assert path == GUIDE_FILE["s3_uri"]
        assert sep == "\t"
        return GUIDE_ROWS.copy()

    monkeypatch.setattr("db2_flattener.gather.lattice.pd.read_csv", fake_read_csv)
    result = read_tabular_file(GUIDE_FILE, connection=None)
    assert list(result["guide_id"]) == ["g1", "g2"]


def test_read_tabular_file_falls_back_to_lattice(monkeypatch, capsys):
    real_read_csv = pd.read_csv

    def fake_read_csv(path, sep=None):
        if isinstance(path, str) and path.startswith("s3://"):
            raise PermissionError("AccessDenied")
        return real_read_csv(path, sep=sep)

    monkeypatch.setattr("db2_flattener.gather.lattice.pd.read_csv", fake_read_csv)
    monkeypatch.setattr(
        "db2_flattener.gather.lattice.download_file",
        lambda object_id, connection: b"guide_id\tguide_role\ng1\ttargeting\n",
    )
    result = read_tabular_file(GUIDE_FILE, connection=SimpleNamespace())
    assert list(result.columns) == ["guide_id", "guide_role"]
    assert list(result["guide_id"]) == ["g1"]
    assert "falling back to Lattice @@download" in capsys.readouterr().out


def test_read_tabular_file_sniffs_tsv_without_format_hints(monkeypatch):
    monkeypatch.setattr(
        "db2_flattener.gather.lattice.download_file",
        lambda object_id, connection: b"guide_id\tguide_role\ng1\ttargeting\n",
    )
    result = read_tabular_file({"@id": "/tabular_files/abc/"}, connection=SimpleNamespace())
    assert list(result.columns) == ["guide_id", "guide_role"]
    assert list(result["guide_id"]) == ["g1"]


def test_resolve_guide_rna_file_returns_complete_embed(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": GUIDE_FILE}
    )
    assert flattener._resolve_guide_rna_file(data) == GUIDE_FILE


def test_resolve_guide_rna_file_id_only_without_gathered_object(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": [GUIDE_FILE["@id"]]}
    )
    assert flattener._resolve_guide_rna_file(data) == {"@id": GUIDE_FILE["@id"]}


def test_resolve_guide_rna_file_prefers_gathered_tabular_file(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": "/genetic_modifications/aaa/", "guide_rna_files": [GUIDE_FILE["@id"]]}
    )
    data["resolved_objects"]["TabularFile"] = {GUIDE_FILE["@id"]: GUIDE_FILE}
    assert flattener._resolve_guide_rna_file(data) == GUIDE_FILE


def test_resolve_guide_rna_file_none_without_warning(monkeypatch, capsys):
    flattener, _calls = make_flattener(monkeypatch)
    assert flattener._resolve_guide_rna_file({"libraries": {}}) is None
    assert "Warning" not in capsys.readouterr().out


def test_create_guide_metadata_reuses_resolved_file_info(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    result = flattener.create_guide_metadata_dataframe(GUIDE_FILE)
    assert result is not None
    assert calls == [GUIDE_FILE]


def test_read_tabular_file_requires_id_when_s3_fails(monkeypatch):
    def fake_read_csv(path, sep=None):
        raise PermissionError("AccessDenied")

    monkeypatch.setattr("db2_flattener.gather.lattice.pd.read_csv", fake_read_csv)
    with pytest.raises(ValueError, match="neither a readable s3_uri nor an @id"):
        read_tabular_file({"s3_uri": "s3://bucket/missing.tsv"}, connection=None)
