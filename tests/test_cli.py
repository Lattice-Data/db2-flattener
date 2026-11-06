import pytest

from db2_flattener import __version__, cli


class FakeConnection:
    def __init__(self, mode):
        pass


class FakeFlattener:
    return_value = "ABC-123_20260814_095000"
    error = None
    last_biohub_versions = None
    last_audit_mode = None

    def __init__(self, connection, configs):
        pass

    def flatten_matrix_file_set(self, uuid, *, biohub_versions=None, audit_mode="strict"):
        FakeFlattener.last_biohub_versions = biohub_versions
        FakeFlattener.last_audit_mode = audit_mode
        if self.error is not None:
            raise self.error
        return self.return_value


ORDER = "s3://bucket/orders/AN1"


def _stub_cli(monkeypatch, flattener_cls=FakeFlattener, argv=None, latest_versions=None):
    monkeypatch.setattr(
        cli,
        "latest_biohub_versions_from_s3",
        lambda s3_order_dir: dict(latest_versions or {}),
    )
    monkeypatch.setattr(cli, "Connection", FakeConnection)
    monkeypatch.setattr(cli, "load_and_return_constant_dicts", lambda mode: ({}, {}))
    monkeypatch.setattr(cli, "DB2Flattener", flattener_cls)
    FakeFlattener.last_biohub_versions = None
    FakeFlattener.last_audit_mode = None
    monkeypatch.setattr(
        "sys.argv",
        argv or ["db2-flattener", "-u", "abc-uuid", "-s", ORDER],
    )


def test_cli_version_exits_without_uuid(monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["db2-flattener", "--version"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_cli_success_prints_path(monkeypatch, capsys):
    _stub_cli(monkeypatch)
    cli.main()
    out = capsys.readouterr().out
    assert f"db2-flattener {__version__}" in out
    assert "No BIOHUB templates in metadata/; all templates start at v1" in out
    assert "Success! Output directory created: ABC-123_20260814_095000" in out
    assert FakeFlattener.last_biohub_versions == {}
    assert FakeFlattener.last_audit_mode == "strict"


def test_cli_s3_order_dir_passes_latest_versions(monkeypatch, capsys):
    latest = {("samples", 1): 2, ("guide", 1): 1}
    _stub_cli(monkeypatch, latest_versions=latest)
    cli.main()
    out = capsys.readouterr().out
    assert FakeFlattener.last_biohub_versions == latest
    assert (
        "Latest BIOHUB templates in metadata/: "
        "guide_template_set1_v1, samples_template_set1_v2" in out
    )


def test_cli_s3_error_exits_1(monkeypatch, capsys):
    def boom(s3_order_dir):
        raise FileNotFoundError(f"S3 order directory not found: {s3_order_dir}")

    _stub_cli(monkeypatch)
    monkeypatch.setattr(cli, "latest_biohub_versions_from_s3", boom)
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 1
    assert "S3 order directory not found" in capsys.readouterr().out


def test_cli_requires_s3_order_dir(monkeypatch):
    monkeypatch.setattr("sys.argv", ["db2-flattener", "-u", "abc-uuid"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2


@pytest.mark.parametrize("flag", ["-b", "--biohub-version"])
def test_cli_rejects_removed_biohub_version_flag(monkeypatch, flag):
    monkeypatch.setattr(
        "sys.argv",
        ["db2-flattener", "-u", "abc-uuid", "-s", ORDER, flag, "v1"],
    )
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2


def test_cli_force_all_audit_sets_mode(monkeypatch):
    _stub_cli(
        monkeypatch,
        argv=["db2-flattener", "-u", "abc-uuid", "-s", ORDER, "--force-all-audit"],
    )
    cli.main()
    assert FakeFlattener.last_audit_mode == "all_ok"


def test_cli_force_only_warning_sets_mode(monkeypatch):
    _stub_cli(
        monkeypatch,
        argv=["db2-flattener", "-u", "abc-uuid", "-s", ORDER, "--force-only-warning"],
    )
    cli.main()
    assert FakeFlattener.last_audit_mode == "warnings_ok"


def test_cli_rejects_both_audit_force_flags(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        [
            "db2-flattener",
            "-u",
            "abc-uuid",
            "-s",
            ORDER,
            "--force-all-audit",
            "--force-only-warning",
        ],
    )
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2


def test_cli_no_data_exits_1(monkeypatch, capsys):
    class NoDataFlattener(FakeFlattener):
        return_value = None

    _stub_cli(monkeypatch, NoDataFlattener)
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 1
    assert "Failed to process MatrixFileSet abc-uuid" in capsys.readouterr().out


def test_cli_exception_exits_1(monkeypatch, capsys):
    class ErrorFlattener(FakeFlattener):
        error = RuntimeError("boom")

    _stub_cli(monkeypatch, ErrorFlattener)
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 1
    assert "Error: boom" in capsys.readouterr().out
