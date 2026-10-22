# db2-flattener

[![CI](https://github.com/Lattice-Data/db2-flattener/actions/workflows/ci.yml/badge.svg)](https://github.com/Lattice-Data/db2-flattener/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/Lattice-Data/db2-flattener/badge.svg?branch=main)](https://coveralls.io/github/Lattice-Data/db2-flattener?branch=main)

Flattener utility for the [Lattice Database](https://data.lattice-data.org/).
It gathers a MatrixFileSet from DB2 and writes MAIN, BIOHUB, GEO, SRA_BIOSAMPLE,
SRA_FILE, SAMPLES, and GUIDE_METADATA CSVs.
Ported from [lattice-tools](https://github.com/Lattice-Data/lattice-tools)
`TOOLS-285-GEO-Flattener` at
[`d4de994`](https://github.com/Lattice-Data/lattice-tools/commit/d4de994638e44e79c04dffed13e7a8b213955fd6).

Runs on Python 3.10 and up, through 3.14 (`requires-python = ">=3.10,<3.15"`).

## Install (dev)

From the repository root:

```bash
pip install -e ".[dev]"
```

Pin a release (requires the `v0.1.0` tag on GitHub):

```bash
pip install "db2-flattener @ git+https://github.com/Lattice-Data/db2-flattener.git@v0.1.0"
```

Check the running code with `db2-flattener --version`.

Set Lattice API credentials as environment variables. Names follow the `--mode`
value, which must start with `db2_` (for example `db2_demo`):

```bash
export DB2_DEMO_KEY=...
export DB2_DEMO_SECRET=...
export DB2_DEMO_SERVER=https://lattice-api-dev.demo.lattice-data.org/
```

## Flatten from this directory

After the editable install, any of these work from the repository root. Output
is written under a directory in the current working directory named
`{CRO_orders}_{timestamp}` (unique `sequence_file_sets_CRO_order` values from
MAIN, sorted and joined by `_`). If no CRO_order values are present, the
directory is `MatrixFileSet_{uuid8}_{timestamp}`.

Exactly one of these BIOHUB version flags is required:

- `--biohub-version` / `-b` — explicit version (`1` or `v1`)
- `--s3-order-dir` / `-s` — S3 URI of the order subdirectory; `v#` is taken from
  `{order}/metadata/` (`samples_template(_set#)?_v#.csv` /
  `guide_template(_set#)?_v#.csv`). Missing/empty metadata → `v1`; otherwise
  existing `v# + 1`. Errors if metadata/ has mixed versions.

```bash
python flatten.py -u <matrix-file-set-uuid> -m db2_demo -b v1
db2-flattener -u <matrix-file-set-uuid> -m db2_demo -s s3://bucket/path/to/order/
python -m db2_flattener -u <matrix-file-set-uuid> -m db2_demo -b 2
```

Layout:

```
{CRO_orders}_{timestamp}/
  MAIN.csv
  GEO/
    geo_metadata.csv                     # includes genetic_modification_alias when present
    demultiplexing_samples.csv           # omitted when empty; same alias column when present
    {gm_alias}_guide_metadata.csv        # copy of each BIOHUB guide template
  SRA/
    sra_biosample.csv
    sra_file.csv
  BIOHUB/
    samples_template_v#.csv              # when no GeneticModifications
    samples_template_setN_v#.csv         # one per GeneticModification
    guide_template_setN_v#.csv           # matching guide TabularFile; omitted when empty/missing
```


## Fetch latest schema from Lattice

After the editable install, the following can be used to run
`src/db2_flattener/schema/generate.py`

```bash
db2-flattener-generate-constants [OPTIONS]
```

Use the `--help` flag for further info on args.

By default this fetches the schema profiles from prod and demo, prints status
and general diffs to the terminal, and updates
`src/db2_flattener/schema/data/constants.yaml` if changes are
found.

## Test

```bash
pytest
```

With coverage:

```bash
pytest --cov=db2_flattener --cov-report=term-missing
```

CI runs the suite on Python 3.10 and 3.14 — the floor declared by
`requires-python` and the newest release. Coverage is uploaded to Coveralls
once per commit, from the 3.10 job.

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
