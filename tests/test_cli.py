import pytest

from db2_flattener import cli


class FakeConnection:
    def __init__(self, mode):
        pass


class FakeFlattener:
    return_value = "myrun_MAIN.csv"
    error = None

    def __init__(self, connection, configs):
        pass

    def flatten_matrix_file_set(self, uuid, output):
        if self.error is not None:
            raise self.error
        return self.return_value


def _stub_cli(monkeypatch, flattener_cls=FakeFlattener):
    monkeypatch.setattr(cli, "Connection", FakeConnection)
    monkeypatch.setattr(cli, "load_and_return_constant_dicts", lambda mode: ({}, {}))
    monkeypatch.setattr(cli, "DB2Flattener", flattener_cls)
    monkeypatch.setattr("sys.argv", ["db2-flattener", "-u", "abc-uuid"])


def test_cli_success_prints_path(monkeypatch, capsys):
    _stub_cli(monkeypatch)
    cli.main()
    assert "Success! CSV file created: myrun_MAIN.csv" in capsys.readouterr().out


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
