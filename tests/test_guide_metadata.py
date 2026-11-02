import warnings
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

GUIDE_FILE_OTHER = {
    "@id": "/tabular_files/22222222-2222-2222-2222-222222222222/",
    "s3_uri": "s3://submissions-lattice/BCP/tabular_files/other.tsv",
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

GM_AAA = "/genetic_modifications/aaa/"
GM_BBB = "/genetic_modifications/bbb/"
LIB_A = "/droplet_based_libraries/lib-a/"
LIB_B = "/droplet_based_libraries/lib-b/"


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


def complete_data_with_gms(*gms, libraries=None):
    normalized = []
    for gm in gms:
        gm = dict(gm)
        if "aliases" not in gm:
            slug = gm["@id"].rstrip("/").rsplit("/", 1)[-1]
            gm["aliases"] = [f"lab:gm-{slug}"]
        normalized.append(gm)
    data = {
        "resolved_objects": {
            "GeneticModification": {gm["@id"]: gm for gm in normalized},
        },
        "libraries": libraries or {},
    }
    return data


def library_entry(lib_atid, *gm_atids):
    return {
        "library": {"@id": lib_atid},
        "genetic_modifications": [{"@id": gm_id} for gm_id in gm_atids],
    }


def guide_df_from_first_set(flattener, data):
    sets = flattener._biohub_guide_sets(data)
    if not sets or not sets[0]["file_info"]:
        return None
    return flattener.create_guide_metadata_dataframe(sets[0]["file_info"])


def test_no_genetic_modifications_returns_empty_sets(monkeypatch):
    flattener, calls = make_flattener(monkeypatch)
    assert flattener._biohub_guide_sets({"libraries": {}}) == []
    assert flattener._biohub_guide_sets({"resolved_objects": {}}) == []
    assert flattener.create_guide_metadata_dataframe(None) is None
    assert calls == []


def test_gm_without_guide_rna_files_skips_guide(monkeypatch, capsys):
    flattener, calls = make_flattener(monkeypatch)
    data = complete_data_with_gms({"@id": GM_AAA, "strategy": "interference screen"})
    sets = flattener._biohub_guide_sets(data)
    assert len(sets) == 1
    assert sets[0]["file_info"] is None
    assert sets[0]["set_number"] == 1
    assert guide_df_from_first_set(flattener, data) is None
    assert calls == []
    assert "has no guide_rna_files" in capsys.readouterr().out


def test_one_unique_file_subsets_columns(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {
            "@id": GM_AAA,
            "guide_rna_files": GUIDE_FILE,
        }
    )
    result = guide_df_from_first_set(flattener, data)
    assert list(result.columns) == [
        "id",
        "protospacer_sequence",
        "role",
        "protospacer_adjacent_motif",
        "intended_features_feature_id",
        "intended_features_gene_name",
    ]
    assert list(result["id"]) == ["g1", "g2"]
    assert len(calls) == 1


def test_same_file_on_two_gms_warns_and_skips_second_guide(monkeypatch, capsys):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {"@id": GM_AAA, "guide_rna_files": [GUIDE_FILE]},
        {"@id": GM_BBB, "guide_rna_files": [GUIDE_FILE]},
    )
    sets = flattener._biohub_guide_sets(data)
    assert len(sets) == 2
    assert sets[0]["gm_id"] == GM_AAA
    assert sets[0]["file_info"] == GUIDE_FILE
    assert sets[1]["gm_id"] == GM_BBB
    assert sets[1]["file_info"] is None
    assert "shared by" in capsys.readouterr().out
    assert calls == []


def test_two_distinct_files_yield_two_sets(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {"@id": GM_AAA, "guide_rna_files": GUIDE_FILE},
        {"@id": GM_BBB, "guide_rna_files": GUIDE_FILE_OTHER},
    )
    sets = flattener._biohub_guide_sets(data)
    assert len(sets) == 2
    assert sets[0]["set_number"] == 1
    assert sets[0]["file_info"] == GUIDE_FILE
    assert sets[1]["set_number"] == 2
    assert sets[1]["file_info"] == GUIDE_FILE_OTHER
    assert calls == []


def test_gm_with_multiple_guide_files_skips_guide(monkeypatch, capsys):
    flattener, _calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {
            "@id": GM_AAA,
            "guide_rna_files": [GUIDE_FILE, GUIDE_FILE_OTHER],
        }
    )
    sets = flattener._biohub_guide_sets(data)
    assert len(sets) == 1
    assert sets[0]["file_info"] is None
    assert "expected one" in capsys.readouterr().out


def test_guide_rna_files_as_id_strings_use_gathered_tabular_file(monkeypatch):
    flattener, calls = make_flattener(monkeypatch, read_df=GUIDE_ROWS)
    data = complete_data_with_gms(
        {
            "@id": GM_AAA,
            "guide_rna_files": [GUIDE_FILE["@id"]],
        }
    )
    data["resolved_objects"]["TabularFile"] = {GUIDE_FILE["@id"]: GUIDE_FILE}
    result = guide_df_from_first_set(flattener, data)
    assert result is not None
    assert calls == [GUIDE_FILE]


def test_missing_expected_columns_returns_none(monkeypatch, capsys):
    flattener, _calls = make_flattener(monkeypatch, read_df=pd.DataFrame({"unrelated": [1]}))
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    assert guide_df_from_first_set(flattener, data) is None
    assert "none of the expected GUIDE_METADATA columns" in capsys.readouterr().out


def test_keeps_only_present_guide_columns(monkeypatch):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame({"guide_id": ["g1"], "guide_role": ["targeting"]}),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result.columns) == ["id", "role"]


def test_control_guide_with_gene_values_blanks_and_warns(monkeypatch, capsys):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["c1"],
                "guide_role": ["Control"],
                "guide_target_gene_id": ["ENSG1"],
                "guide_target_gene_name": ["GENE1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["intended_features_feature_id"]) == [""]
    assert list(result["intended_features_gene_name"]) == [""]
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "intended_features_feature_id" in out
    assert "ENSG1" in out
    assert "intended_features_gene_name" in out
    assert "GENE1" in out
    assert "id=c1" in out


@pytest.mark.parametrize("placeholder", ["NA", "na", "N/A", "n/a", "", None])
def test_control_guide_na_placeholders_blank_without_warning(monkeypatch, capsys, placeholder):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["c1"],
                "guide_role": ["Control"],
                "guide_target_gene_id": [placeholder],
                "guide_target_gene_name": [placeholder],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["intended_features_feature_id"]) == [""]
    assert list(result["intended_features_gene_name"]) == [""]
    assert "====WARNING====" not in capsys.readouterr().out


def test_control_guide_float_nan_columns_blank_without_futurewarning(monkeypatch, capsys):
    """read_csv all-NaN intended-feature columns are float64; blanking must not warn."""
    read_df = pd.DataFrame(
        {
            "guide_id": ["c1"],
            "guide_role": ["Control"],
            "guide_target_gene_id": [float("nan")],
            "guide_target_gene_name": [float("nan")],
        }
    )
    assert read_df["guide_target_gene_id"].dtype == float
    flattener, _calls = make_flattener(monkeypatch, read_df=read_df)
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        result = guide_df_from_first_set(flattener, data)
    assert list(result["intended_features_feature_id"]) == [""]
    assert list(result["intended_features_gene_name"]) == [""]
    assert "====WARNING====" not in capsys.readouterr().out


def test_targeting_guide_intended_features_unchanged(monkeypatch, capsys):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["t1"],
                "guide_role": ["Targeting a Gene"],
                "guide_target_gene_id": ["ENSG1"],
                "guide_target_gene_name": ["GENE1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["intended_features_feature_id"]) == ["ENSG1"]
    assert list(result["intended_features_gene_name"]) == ["GENE1"]
    assert "====WARNING====" not in capsys.readouterr().out


def test_mixed_roles_only_control_intended_features_cleared(monkeypatch, capsys):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["t1", "c1"],
                "guide_role": ["Targeting a Gene", "Control"],
                "guide_target_gene_id": ["ENSG1", "ENSG2"],
                "guide_target_gene_name": ["GENE1", "GENE2"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["intended_features_feature_id"]) == ["ENSG1", ""]
    assert list(result["intended_features_gene_name"]) == ["GENE1", ""]
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "id=c1" in out
    assert "ENSG2" in out


def test_control_blanking_handles_missing_intended_feature_column(monkeypatch, capsys):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["c1"],
                "guide_role": ["Control"],
                "guide_target_gene_id": ["ENSG1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result.columns) == ["id", "role", "intended_features_feature_id"]
    assert list(result["intended_features_feature_id"]) == [""]
    out = capsys.readouterr().out
    assert "intended_features_feature_id" in out
    assert "intended_features_gene_name" not in out


def test_missing_role_skips_control_blanking(monkeypatch, capsys):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["g1"],
                "guide_target_gene_id": ["ENSG1"],
                "guide_target_gene_name": ["GENE1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["intended_features_feature_id"]) == ["ENSG1"]
    assert list(result["intended_features_gene_name"]) == ["GENE1"]
    assert "====WARNING====" not in capsys.readouterr().out


@pytest.mark.parametrize("alias", ["control", "CONTROL"])
def test_control_role_aliases_normalized_and_blanked(monkeypatch, capsys, alias):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["c1"],
                "guide_role": [alias],
                "guide_target_gene_id": ["ENSG1"],
                "guide_target_gene_name": ["GENE1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["role"]) == ["Control"]
    assert list(result["intended_features_feature_id"]) == [""]
    assert list(result["intended_features_gene_name"]) == [""]
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "Control guide(s) had non-empty" in out


@pytest.mark.parametrize(
    "role",
    ["Targeting a Gene", "Targeting a Region", "Control"],
)
def test_allowed_guide_roles_unchanged_without_role_warning(monkeypatch, capsys, role):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["g1"],
                "guide_role": [role],
                "guide_target_gene_id": [""],
                "guide_target_gene_name": [""],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["role"]) == [role]
    assert "====WARNING====" not in capsys.readouterr().out


def test_invalid_guide_role_warns_and_leaves_value(monkeypatch, capsys):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["g1"],
                "guide_role": ["non-targeting"],
                "guide_target_gene_id": ["ENSG1"],
                "guide_target_gene_name": ["GENE1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    assert list(result["role"]) == ["non-targeting"]
    assert list(result["intended_features_feature_id"]) == ["ENSG1"]
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "id=g1" in out
    assert "non-targeting" in out


@pytest.mark.parametrize("empty_role", ["", None])
def test_empty_guide_role_warns_and_leaves_value(monkeypatch, capsys, empty_role):
    flattener, _calls = make_flattener(
        monkeypatch,
        read_df=pd.DataFrame(
            {
                "guide_id": ["g1"],
                "guide_role": [empty_role],
                "guide_target_gene_id": ["ENSG1"],
                "guide_target_gene_name": ["GENE1"],
            }
        ),
    )
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": GUIDE_FILE})
    result = guide_df_from_first_set(flattener, data)
    # None becomes NaN in the DataFrame; leave whatever is there after read
    left = result["role"].iloc[0]
    assert left == empty_role or (empty_role is None and pd.isna(left))
    out = capsys.readouterr().out
    assert "====WARNING====" in out
    assert "id=g1" in out


def test_biohub_sets_map_libraries_to_gms(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": GM_AAA, "guide_rna_files": GUIDE_FILE},
        {"@id": GM_BBB, "guide_rna_files": GUIDE_FILE_OTHER},
        libraries={
            "uuid-a": library_entry(LIB_A, GM_AAA),
            "uuid-b": library_entry(LIB_B, GM_BBB),
        },
    )
    sets = flattener._biohub_guide_sets(data)
    assert sets[0]["library_ids"] == {LIB_A}
    assert sets[1]["library_ids"] == {LIB_B}


def test_biohub_sets_rejects_library_linked_to_multiple_gms(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": GM_AAA, "guide_rna_files": GUIDE_FILE},
        {"@id": GM_BBB, "guide_rna_files": GUIDE_FILE_OTHER},
        libraries={"uuid-a": library_entry(LIB_A, GM_AAA, GM_BBB)},
    )
    with pytest.raises(ValueError, match=rf"{LIB_A}.*{GM_AAA}.*{GM_BBB}"):
        flattener._biohub_guide_sets(data)


def test_biohub_sets_sorted_by_gm_id(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {"@id": GM_BBB, "guide_rna_files": GUIDE_FILE_OTHER},
        {"@id": GM_AAA, "guide_rna_files": GUIDE_FILE},
    )
    sets = flattener._biohub_guide_sets(data)
    assert [s["gm_id"] for s in sets] == [GM_AAA, GM_BBB]
    assert [s["set_number"] for s in sets] == [1, 2]


def test_biohub_sets_include_cleaned_alias(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {
            "@id": GM_AAA,
            "aliases": ["lab:gm-screen-a"],
            "guide_rna_files": GUIDE_FILE,
        }
    )
    sets = flattener._biohub_guide_sets(data)
    assert sets[0]["alias"] == "gm-screen-a"


def test_biohub_sets_sanitize_alias_for_filename(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {
            "@id": GM_AAA,
            "aliases": ["lab:foo/bar baz"],
            "guide_rna_files": GUIDE_FILE,
        }
    )
    sets = flattener._biohub_guide_sets(data)
    assert sets[0]["alias"] == "foo_bar_baz"


def test_biohub_sets_reject_multiple_aliases(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms(
        {
            "@id": GM_AAA,
            "aliases": ["lab:gm-a", "lab:gm-b"],
            "guide_rna_files": GUIDE_FILE,
        }
    )
    with pytest.raises(ValueError, match="exactly one alias"):
        flattener._biohub_guide_sets(data)


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


def test_biohub_set_prefers_gathered_tabular_file(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": [GUIDE_FILE["@id"]]})
    data["resolved_objects"]["TabularFile"] = {GUIDE_FILE["@id"]: GUIDE_FILE}
    sets = flattener._biohub_guide_sets(data)
    assert sets[0]["file_info"] == GUIDE_FILE


def test_biohub_set_id_only_without_gathered_object(monkeypatch):
    flattener, _calls = make_flattener(monkeypatch)
    data = complete_data_with_gms({"@id": GM_AAA, "guide_rna_files": [GUIDE_FILE["@id"]]})
    sets = flattener._biohub_guide_sets(data)
    assert sets[0]["file_info"] == {"@id": GUIDE_FILE["@id"]}


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
