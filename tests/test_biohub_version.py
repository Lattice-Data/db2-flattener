import pytest

from db2_flattener.flatten.biohub_version import (
    latest_biohub_versions_from_s3,
    latest_versions_from_metadata_names,
    next_biohub_version,
)


def test_latest_versions_tracks_kind_and_set_independently():
    names = [
        "bucket/order/metadata/samples_template_set1_v2.csv",
        "samples_template_set1_v1.csv",
        "samples_template_set2_v1.csv",
        "guide_template_set1_v1.csv",
        "guide_template_set2_v4.csv",
    ]
    assert latest_versions_from_metadata_names(names) == {
        ("samples", 1): 2,
        ("samples", 2): 1,
        ("guide", 1): 1,
        ("guide", 2): 4,
    }


def test_latest_versions_ignores_names_without_set():
    names = [
        "samples_template_v9.csv",
        "guide_template_v9.csv",
        "samples_template_set1_v2.csv",
    ]
    assert latest_versions_from_metadata_names(names) == {("samples", 1): 2}


@pytest.mark.parametrize(
    "name",
    [
        "readme.txt",
        "samples_template_set1.csv",
        "samples_template_set1_v2.tsv",
        "samples_template_set_v2.csv",
        "samples_template_set1_v.csv",
        "old_samples_template_set1_v2.csv",
        "samples_template_set1_v2.csv.bak",
        "guide_metadata_set1_v2.csv",
    ],
)
def test_latest_versions_ignores_non_matching(name):
    assert latest_versions_from_metadata_names([name]) == {}


def test_latest_versions_case_insensitive():
    names = ["Samples_Template_Set1_V3.CSV", "GUIDE_TEMPLATE_SET2_V2.csv"]
    assert latest_versions_from_metadata_names(names) == {
        ("samples", 1): 3,
        ("guide", 2): 2,
    }


def test_latest_versions_does_not_fill_gaps():
    names = ["samples_template_set1_v1.csv", "samples_template_set1_v3.csv"]
    assert latest_versions_from_metadata_names(names) == {("samples", 1): 3}


def test_next_version_absent_is_v1():
    assert next_biohub_version({}, "samples", 1) == "v1"
    assert next_biohub_version(None, "guide", 3) == "v1"


def test_next_version_independent_per_kind_and_set():
    latest = {("samples", 1): 2, ("guide", 1): 1, ("samples", 2): 5}
    assert next_biohub_version(latest, "samples", 1) == "v3"
    assert next_biohub_version(latest, "guide", 1) == "v2"
    assert next_biohub_version(latest, "samples", 2) == "v6"
    assert next_biohub_version(latest, "guide", 2) == "v1"


def test_next_version_rejects_unknown_kind():
    with pytest.raises(ValueError, match="unknown BIOHUB template kind"):
        next_biohub_version({}, "guide_metadata", 1)


class FakeFS:
    def __init__(self, *, exists_paths=None, ls_map=None, exists_errors=None, ls_errors=None):
        self.exists_paths = set(exists_paths or [])
        self.ls_map = dict(ls_map or {})
        self.exists_errors = dict(exists_errors or {})
        self.ls_errors = dict(ls_errors or {})

    def exists(self, path):
        path = path.rstrip("/")
        if path in self.exists_errors:
            raise self.exists_errors[path]
        return path in self.exists_paths or any(
            p == path or p.startswith(path + "/") for p in self.exists_paths
        )

    def ls(self, path, detail=False):
        path = path.rstrip("/")
        if path in self.ls_errors:
            raise self.ls_errors[path]
        return list(self.ls_map.get(path, []))


ORDER = "s3://bucket/orders/AN0001"
METADATA = f"{ORDER}/metadata"


def test_s3_missing_metadata_returns_empty():
    fs = FakeFS(exists_paths={ORDER})
    assert latest_biohub_versions_from_s3(ORDER, fs=fs) == {}


def test_s3_empty_metadata_returns_empty():
    fs = FakeFS(exists_paths={ORDER, METADATA}, ls_map={METADATA: []})
    assert latest_biohub_versions_from_s3(ORDER, fs=fs) == {}


def test_s3_no_matching_names_returns_empty():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_map={METADATA: [f"{METADATA}/notes.txt", f"{METADATA}/samples_template_v2.csv"]},
    )
    assert latest_biohub_versions_from_s3(ORDER, fs=fs) == {}


def test_s3_samples_and_guide_resolve_independently():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_map={
            METADATA: [
                f"{METADATA}/samples_template_set1_v1.csv",
                f"{METADATA}/samples_template_set1_v2.csv",
                f"{METADATA}/guide_template_set1_v1.csv",
            ]
        },
    )
    latest = latest_biohub_versions_from_s3(ORDER, fs=fs)
    assert next_biohub_version(latest, "samples", 1) == "v3"
    assert next_biohub_version(latest, "guide", 1) == "v2"


def test_s3_trailing_slash_on_order_dir():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_map={METADATA: [f"{METADATA}/samples_template_set1_v1.csv"]},
    )
    assert latest_biohub_versions_from_s3(ORDER + "/", fs=fs) == {("samples", 1): 1}


def test_s3_order_dir_missing_raises():
    fs = FakeFS(exists_paths=set())
    with pytest.raises(FileNotFoundError, match="not found"):
        latest_biohub_versions_from_s3(ORDER, fs=fs)


def test_s3_order_dir_access_error_raises():
    fs = FakeFS(exists_errors={ORDER: PermissionError("denied")})
    with pytest.raises(RuntimeError, match="Failed to access"):
        latest_biohub_versions_from_s3(ORDER, fs=fs)


def test_s3_metadata_list_error_raises():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_errors={METADATA: PermissionError("denied")},
    )
    with pytest.raises(RuntimeError, match="Failed to list"):
        latest_biohub_versions_from_s3(ORDER, fs=fs)


def test_s3_rejects_non_s3_uri():
    with pytest.raises(ValueError, match="s3://"):
        latest_biohub_versions_from_s3("/local/path", fs=FakeFS())
