import pytest

from db2_flattener import __version__, cli


class FakeConnection:
    def __init__(self, mode):
        pass


class FakeFlattener:
    return_value = "ABC-123_20260814_095000"
    error = None
    last_biohub_version = None
    last_audit_mode = None

    def __init__(self, connection, configs):
        pass

    def flatten_matrix_file_set(self, uuid, *, biohub_version, audit_mode="strict"):
        FakeFlattener.last_biohub_version = biohub_version
        FakeFlattener.last_audit_mode = audit_mode
        if self.error is not None:
            raise self.error
        return self.return_value


def _stub_cli(monkeypatch, flattener_cls=FakeFlattener, argv=None):
    monkeypatch.setattr(cli, "Connection", FakeConnection)
    monkeypatch.setattr(cli, "load_and_return_constant_dicts", lambda mode: ({}, {}))
    monkeypatch.setattr(cli, "DB2Flattener", flattener_cls)
    FakeFlattener.last_biohub_version = None
    FakeFlattener.last_audit_mode = None
    monkeypatch.setattr(
        "sys.argv",
        argv or ["db2-flattener", "-u", "abc-uuid", "-b", "v1"],
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
    assert "BIOHUB template version: v1" in out
    assert "Success! Output directory created: ABC-123_20260814_095000" in out
    assert FakeFlattener.last_biohub_version == "v1"
    assert FakeFlattener.last_audit_mode == "strict"


def test_cli_biohub_version_normalizes(monkeypatch, capsys):
    _stub_cli(monkeypatch, argv=["db2-flattener", "-u", "abc-uuid", "-b", "3"])
    cli.main()
    assert FakeFlattener.last_biohub_version == "v3"
    assert "BIOHUB template version: v3" in capsys.readouterr().out


def test_cli_s3_order_dir_resolves_version(monkeypatch, capsys):
    monkeypatch.setattr(
        cli,
        "resolve_biohub_version",
        lambda **kwargs: "v4",
    )
    _stub_cli(
        monkeypatch,
        argv=["db2-flattener", "-u", "abc-uuid", "-s", "s3://bucket/orders/AN1"],
    )
    cli.main()
    assert FakeFlattener.last_biohub_version == "v4"


def test_cli_requires_version_flag(monkeypatch):
    monkeypatch.setattr("sys.argv", ["db2-flattener", "-u", "abc-uuid"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2


def test_cli_rejects_both_version_flags(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        ["db2-flattener", "-u", "abc-uuid", "-b", "v1", "-s", "s3://bucket/o"],
    )
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 2


def test_cli_force_all_audit_sets_mode(monkeypatch):
    _stub_cli(
        monkeypatch,
        argv=["db2-flattener", "-u", "abc-uuid", "-b", "v1", "--force_ALL_AUDIT"],
    )
    cli.main()
    assert FakeFlattener.last_audit_mode == "all_ok"


def test_cli_force_only_warning_sets_mode(monkeypatch):
    _stub_cli(
        monkeypatch,
        argv=["db2-flattener", "-u", "abc-uuid", "-b", "v1", "--force_ONLY_WARNING"],
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
            "-b",
            "v1",
            "--force_ALL_AUDIT",
            "--force_ONLY_WARNING",
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
