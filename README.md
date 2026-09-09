# db2-flattener

[![CI](https://github.com/Lattice-Data/db2-flattener/actions/workflows/ci.yml/badge.svg)](https://github.com/Lattice-Data/db2-flattener/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/Lattice-Data/db2-flattener/badge.svg?branch=main)](https://coveralls.io/github/Lattice-Data/db2-flattener?branch=main)

Flattener utility for the [Lattice Database](https://data.lattice-data.org/).
It gathers a MatrixFileSet from DB2 and writes MAIN, BIOHUB, GEO, SRA_BIOSAMPLE,
SAMPLES, and GUIDE_METADATA CSVs.
Ported from [lattice-tools](https://github.com/Lattice-Data/lattice-tools)
`TOOLS-285-GEO-Flattener` at
[`d4de994`](https://github.com/Lattice-Data/lattice-tools/commit/d4de994638e44e79c04dffed13e7a8b213955fd6).

Runs on Python 3.10 and up, through 3.14 (`requires-python = ">=3.10,<3.15"`).

## Install (dev)

From the repository root:

```bash
pip install -e ".[dev]"
```

Set Lattice API credentials as environment variables. Names follow the `--mode`
value, which must start with `db2_` (for example `db2_demo`):

```bash
export DB2_DEMO_KEY=...
export DB2_DEMO_SECRET=...
export DB2_DEMO_SERVER=https://lattice-api-dev.demo.lattice-data.org/
```

## Flatten from this directory

After the editable install, any of these work from the repository root. CSVs
are written to the current working directory.

```bash
python flatten.py -u <matrix-file-set-uuid> -m db2_demo
db2-flattener -u <matrix-file-set-uuid> -m db2_demo
python -m db2_flattener -u <matrix-file-set-uuid> -m db2_demo
```

Optional `-o` sets a shared path prefix. CSVs are written as
`{prefix}_MAIN.csv`, `{prefix}_BIOHUB.csv`, `{prefix}_GEO.csv`,
`{prefix}_SRA_BIOSAMPLE.csv`, `{prefix}_SAMPLES.csv`, and
`{prefix}_GUIDE_METADATA.csv`. SAMPLES and GUIDE_METADATA are omitted when
there is nothing to write. Without `-o`, the prefix is
`MatrixFileSet_{uuid}_{timestamp}`.

SRA_BIOSAMPLE is one row per CRO group, collapsing that group's MAIN rows.
`sample_name` holds the library's `CRO_group_identifier`. That identifier spans
every library in a group, so a group's GEX and CRISPR libraries would land on the
same row; only the GEX ones are taken, using the same test as GEO. A group with no
GEX library leaves the sheet.

`library_id` carries the per-library identity the key gives up: every library in
the group by cleaned alias, GEX first — `cohort1_batch2_9_GEX,
cohort1_batch2_9_CRI`. It is the one column read from MAIN before the GEX filter,
since the libraries it names are the rows that filter drops.

The biological sample aliases appear only inside
`sample_name: sample_probe_barcode`, which is the library's alias-to-barcode map
sorted by alias — so the `sample_name` inside that column name means the
biological sample, not the group naming the row. The map is read from the
library's own embedded `samples` field, so this sheet does not depend on the
per-sample merge that fills the `tissues_*` columns in MAIN.

Every cell is an unordered set of the distinct values a group's samples carry.
One value is written bare; several are prefixed `pooled:` — `pooled: 889023040,
889081306`. The prefix marks the cell as a set and nothing more: it does **not**
record which donor had which age or ethnicity, and a sample pooling several donors
contributes all of them.

`*isolate`, `*sex` and `ethnicity` come from the donors, `*age` and the age bounds
from the sample. `*sex` keeps its own wording, `pooled male and female`, with male
first. `ethnicity` is `HumanDonor`-only, so a non-human run has no such column.

`age_lower_bound` and `age_upper_bound` are optional, appearing only when some
sample carries that bound. Both pair the bound with its `age_units`, pluralised
(`29 years`, `1 year`), and a sample without one contributes `not provided` to the
set.

`*tissue` is the sample's `sample_terms`, one entry per term. A cell line or
primary cell culture has no tissue of origin, so it reads `not available`.

`*biomaterial_provider` is the `title` of the sample's `sources`, falling back per
row to that sample's `lab`.

`*collection_date` and `*geo_loc_name` are the sample's `date_obtained` and
`collection_geographical_location`, with `not provided` for a sample that has
neither.

The remaining columns are optional, appearing only when some sample has a value.
Being optional decides only whether the column exists, not how it is filled: once
one group has a value, every group gets a cell.

- `suspension_type` — the sample's `suspension_type`.
- `preservation_method` — on tissues alone, so any other sample type reads
  `not applicable`.
- `intended_cell_type` — the sample's `intended_cell_types`, which only cell
  lines and organoids have, so the others read `not applicable`.
- `suspension_enriched_cell_types` and `suspension_depleted_cell_types` — the
  sample's `enriched_cell_types` and `depleted_cell_types`. Both are on all four
  sample types, and a sample enriched for one type and depleted of another
  reports in both columns.
- `suspension_selection_kits` and `suspension_enrichment_factors` — the sample's
  `selection_kits` and `selection_markers`, taking their sources from BIOHUB. The
  two are independent: a run with kits but no markers gets only the first.
- `genetic_perturbation_strategy` — the `strategy` of the sample's linked
  `GeneticModification`, rewritten through `GENETIC_PERTURBATION_MAP` so it reads
  `CRISPR interference screen`, matching BIOHUB. A sample with no modification
  reads `not applicable`.
- `experimental_perturbation` — the treatment's duration joined with its
  description, `8 hour stimulation`, or `8-24 hour stimulation` when the bounds
  differ. A sample with no treatment contributes `no treatment`.
- `experimental_perturbation_factors` — the terms named by the sample's
  treatments, bracketed when there are several so a reader can see which went
  together, `[IL2_HUMAN, anti-CD2_HUMAN]`. A sample with no treatment
  contributes `na`.


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
