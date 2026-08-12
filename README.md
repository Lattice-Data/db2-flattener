# db2-flattener

[![CI](https://github.com/Lattice-Data/db2-flattener/actions/workflows/ci.yml/badge.svg)](https://github.com/Lattice-Data/db2-flattener/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/Lattice-Data/db2-flattener/badge.svg?branch=main)](https://coveralls.io/github/Lattice-Data/db2-flattener?branch=main)

Flattener utilities for DB2.

Requires Python 3.10 through 3.14.

## Install (dev)

```bash
pip install -e ".[dev]"
```

## Test

```bash
pytest
```

With coverage:

```bash
pytest --cov=db2_flattener --cov-report=term-missing
```

CI runs the suite on Python 3.10 and 3.14 — the floor declared by
`requires-python` and the newest release. This is a zero-dependency pure-Python
library, so the versions in between have nothing to break that those two ends
would not also catch. Coverage is uploaded to Coveralls once per commit, from
the 3.10 job.

## Lint and format

Ruff handles both, configured in `pyproject.toml`. Install the pre-commit hook
once and it runs on every commit:

```bash
pip install pre-commit && pre-commit install
```

The hook formats in place and then lints. To run it over the whole tree without
committing:

```bash
pre-commit run --all-files
```

CI runs the same checks as `ruff check .` and `ruff format --check .`.

## License

MIT. Copyright (c) 2026 Lattice Data Coordination. See [LICENSE](LICENSE).
