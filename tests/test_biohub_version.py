import pytest

from db2_flattener.flatten.biohub_version import (
    next_biohub_version_from_versions,
    normalize_biohub_version,
    resolve_biohub_version,
    resolve_biohub_version_from_s3,
    versions_from_metadata_names,
)


@pytest.mark.parametrize(
    "value, expected",
    [
        ("1", "v1"),
        ("v1", "v1"),
        ("V3", "v3"),
        (" 02 ", "v2"),
    ],
)
def test_normalize_biohub_version(value, expected):
    assert normalize_biohub_version(value) == expected


@pytest.mark.parametrize("value", ["", "v", "v0", "0", "x1", "v1.2"])
def test_normalize_biohub_version_rejects_invalid(value):
    with pytest.raises(ValueError):
        normalize_biohub_version(value)


def test_versions_from_metadata_names():
    names = [
        "bucket/order/metadata/samples_template_v2.csv",
        "guide_template_set1_v2.csv",
        "readme.txt",
        "samples_template_set2_v2.csv",
    ]
    assert versions_from_metadata_names(names) == {2}


def test_next_version_empty_is_v1():
    assert next_biohub_version_from_versions(set()) == "v1"


def test_next_version_increments():
    assert next_biohub_version_from_versions({3}) == "v4"


def test_next_version_mixed_raises():
    with pytest.raises(ValueError, match="mixed BIOHUB template versions"):
        next_biohub_version_from_versions({1, 2})


def test_resolve_explicit_version():
    assert resolve_biohub_version(biohub_version="v2") == "v2"
    assert resolve_biohub_version(biohub_version="2") == "v2"


def test_resolve_requires_exactly_one():
    with pytest.raises(ValueError, match="required"):
        resolve_biohub_version()
    with pytest.raises(ValueError, match="only one"):
        resolve_biohub_version(s3_order_dir="s3://b/o", biohub_version="v1")


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


def test_s3_missing_metadata_returns_v1():
    fs = FakeFS(exists_paths={ORDER})
    assert resolve_biohub_version_from_s3(ORDER, fs=fs) == "v1"


def test_s3_empty_metadata_returns_v1():
    fs = FakeFS(exists_paths={ORDER, METADATA}, ls_map={METADATA: []})
    assert resolve_biohub_version_from_s3(ORDER, fs=fs) == "v1"


def test_s3_no_matching_names_returns_v1():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_map={METADATA: [f"{METADATA}/notes.txt"]},
    )
    assert resolve_biohub_version_from_s3(ORDER, fs=fs) == "v1"


def test_s3_single_version_increments():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_map={
            METADATA: [
                f"{METADATA}/samples_template_set1_v2.csv",
                f"{METADATA}/guide_template_set1_v2.csv",
            ]
        },
    )
    assert resolve_biohub_version_from_s3(ORDER, fs=fs) == "v3"


def test_s3_mixed_versions_raise():
    fs = FakeFS(
        exists_paths={ORDER, METADATA},
        ls_map={
            METADATA: [
                f"{METADATA}/samples_template_v1.csv",
                f"{METADATA}/guide_template_v2.csv",
            ]
        },
    )
    with pytest.raises(ValueError, match="mixed"):
        resolve_biohub_version_from_s3(ORDER, fs=fs)


def test_s3_order_dir_missing_raises():
    fs = FakeFS(exists_paths=set())
    with pytest.raises(FileNotFoundError, match="not found"):
        resolve_biohub_version_from_s3(ORDER, fs=fs)


def test_s3_order_dir_access_error_raises():
    fs = FakeFS(exists_errors={ORDER: PermissionError("denied")})
    with pytest.raises(RuntimeError, match="Failed to access"):
        resolve_biohub_version_from_s3(ORDER, fs=fs)


def test_s3_rejects_non_s3_uri():
    with pytest.raises(ValueError, match="s3://"):
        resolve_biohub_version_from_s3("/local/path", fs=FakeFS())
