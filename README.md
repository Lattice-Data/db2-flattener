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

SRA_BIOSAMPLE is one row per library, collapsing that library's MAIN rows. The
submitted "sample" is the sequencing library, so `sample_name` holds the library
alias; the biological sample aliases appear only inside
`sample_name: sample_probe_barcode`, which is the library's alias-to-barcode map
sorted by alias. That map is read from the library's own embedded `samples`
field, so this sheet does not depend on the per-sample merge that fills the
`tissues_*` columns in MAIN.

`*isolate` and `*age` describe the library's donors, enumerated `D1..Dn` in donor
id order — `pooled: D1 - 889023040, D2 - 889081306` and
`pooled: D1 - 32 years, D2 - 29 years`. A library with one donor gets the bare
value and no `pooled:` prefix. Both columns share one enumeration, so `Dn` names
the same donor in each. `*sex` summarises the same donors rather than listing
them, so it carries no labels: one sex gives `female`, a mixed pool gives
`pooled male and female`, with male first and any other value such as `unknown`
pooled the same way after it.

`*biomaterial_provider` is the `title` of the sample's `sources`, falling back
per row to the `title` of that sample's `lab`, so a library whose samples differ
can carry both. Note that `sources` declares a multi-type `linkTo`, which
`extract.get_link_to` only recognises as a single string, so `OBJECT_CONFIG` never
lists it as a reference and the gatherer does not resolve it — when it arrives as
a bare `@id` path there is no title to read and the `lab` fallback applies.

`*collection_date` and `*geo_loc_name` are the sample's `date_obtained` and
`collection_geographical_location`, with `not provided` standing in for any sample
that has none — so a library whose samples disagree reads
`2023-01-05; not provided`.

`*tissue` is the sample's `sample_terms`, one entry per term since the field is an
array. A cell line or primary cell culture has no tissue of origin, so it reports
`not available` instead; tissues and organoids report their term.

`experimental_perturbation` is optional and comes from `treatments_*`, joining the
treatment's duration with its description — `8 hour stimulation`, or
`8-24 hour stimulation` when the bounds differ. Units stay verbatim rather than
pluralised, since the phrase reads adjectivally. One perturbation across a library
gives the bare value; several give `pooled: a, b`, and a sample with no treatment
contributes `not provided`, so a partly perturbed library reads
`pooled: 8 hour stimulation, not provided`.

`experimental_perturbation_factors` is optional and comes from
`treatments_ontological_term`. Each sample contributes the set of terms its
treatments name, bracketed when there is more than one so a reader can see which
factors went together — `[IL2_HUMAN, anti-CD2_HUMAN]`. Distinct sets pool across
the library and a sample with no treatment contributes `na`, giving
`pooled: [IL2_HUMAN, anti-CD2_HUMAN], na`. Terms sort case-sensitively, so
uppercase leads.

`suspension_type` is optional: being optional decides only whether the column
exists, not how it is filled. If no sample anywhere has one the column is dropped;
otherwise it behaves exactly like `*collection_date`, joining distinct values with
`; ` and filling a gap with `not provided`.

`ethnicity` is optional and appears only when some donor has one. It shares the
`D1..Dn` enumeration with `*isolate` and `*age`, so the labels line up across all
three — `pooled: D1 - European American, D2 - African American`. It comes from
`HumanDonor` only, so a non-human run has no such column.

`age_lower_bound` and `age_upper_bound` are optional: each appears only when some
donor carries that bound. Both pair the sample's `lower_bound_age` /
`upper_bound_age` with its `age_units`, pluralised (`29 years`, `1 year`). They
carry no `D1..Dn` labels but hold one entry per donor in the same order, so entry
*n* describes the donor `Dn` names in `*isolate` and `*age`. Nothing is sorted or
deduplicated, and a donor with no bound holds its slot with `not provided` —
`not provided; not provided; 40 years` for three donors where only the third has
one.

Age comes from the sample's `developmental_stages` term,
because `HumanDonor` has no age property of its own: a numeric stage renders as a
number and unit (`29-year-old stage` → `29 years`), and a qualitative one keeps
its term name with a trailing ` stage` or ` human stage` removed (`adult stage` →
`adult`, `10th week post-fertilization human stage` →
`10th week post-fertilization`).

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
