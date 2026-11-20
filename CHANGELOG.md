# Changelog

## Unreleased

## 2.0.0

### Breaking changes

- Requires the database schema with the list-valued `HumanDonor.ethnicities`.
  Versions 1.x run against it silently drop ethnicity from GEO, SRA biosample
  and SAMPLES, and report every human as `unknown` on BIOHUB.
- SRA outputs are tab-separated: `sra_file_set{N}.tsv` and
  `sra_biosample_set{N}.tsv` replace the `.csv` files.
- MAIN.csv: `sample_alias` is removed, and `raw_matrix_file_alias` keeps the lab
  prefix.
- Ethnicity cells use a new format (`pooled: a, b` across donors, `a; b` for one
  donor).
- Flatten now refuses (exit 1) orders it used to accept: multiome or ATAC
  libraries, RawMatrixFiles sharing a `geo_file_name`, and gathered objects with
  `status: deleted`.
- `DB2Flattener.flatten_matrix_file_set` takes `biohub_versions` /
  `biohub_version`, and several `flatten/biohub_version.py` helpers are
  replaced.

### Changes

- `HumanDonor.ethnicity` is replaced by the list-valued `HumanDonor.ethnicities`;
  the flattener reads `human_donors_ethnicities_term_id` / `_term_name` and no
  longer reads `ethnicity`. Before this, ethnicity silently dropped off GEO, SRA
  biosample and SAMPLES, and BIOHUB defaulted every human to `unknown`.
  - BIOHUB `self_reported_ethnicity_ontology_term_id`: several ids are sorted and
    joined with ` || `; an empty list on a human donor is `unknown` (`na` for
    non-human).
  - GEO `donor_ethnicity`, SRA biosample `ethnicity` and SAMPLES
    `donor_ethnicity`: when the library's samples reference more than one
    distinct human donor, several terms read `pooled: a, b` (as before); with
    exactly one donor they read `a; b`. Without `human_donors_@id` the donor
    count is unknown and the pooled form is used.

- The GEO and SAMPLES sheets' `processed data file` columns now hold each
  RawMatrixFile's `geo_file_name` instead of its alias. The main and sample
  dataframes gain a `raw_matrix_files_geo_file_name` column (`raw_matrix_file_alias`
  is still there); a file with no `geo_file_name` leaves its `processed data file`
  cell empty, and the run prints a `====WARNING====` block listing every such
  RawMatrixFile by `@id` (and alias). That block is the one report of a missing
  `geo_file_name`: `md5_checksums.tsv` still leaves the `file name` cell blank
  but no longer prints its own `has no geo_file_name` line (it still warns on a
  missing `md5sum`) (TOOLS-380).
- After writing `md5_checksums.tsv`, the run prints a `====WARNING====` block
  listing every RawMatrixFile that is in `md5_checksums.tsv` but not on the GEO
  sheet, compared by uuid, with its file name (or `no geo_file_name`) and why:
  `no GEX/ATAC library` (e.g. CRISPR-only, so on SAMPLES but not GEO) or
  `not on any library` (TOOLS-380).
- `raw_matrix_file_alias` in MAIN.csv now keeps the lab prefix (`lab:rmf1`, not
  `rmf1`); multiple aliases are joined with `; `. `_get_clean_alias` is removed
  (TOOLS-380).
- The `sample_alias` column is removed from the sample dataframe and MAIN.csv.
  It duplicated the per-type `tissues_aliases`, `cell_lines_aliases`, etc.
  columns, which still carry the same values (TOOLS-380).
- Refuse flatten (exit 1) when two or more RawMatrixFiles share a
  `geo_file_name`, compared case-insensitively (`GSM1.h5` and `gsm1.h5` collide).
  Each shared name prints an `====== ERROR ======` block listing the files'
  `@id` and `s3_uri`, before any output is written. Blank names are not
  counted, and no `--force-*` flag overrides it.

- New tab-separated `GEO/md5_checksums.tsv`: one row per RawMatrixFile of the
  MatrixFileSet,
  with the GEO submission template columns `file name` (the calculated
  `geo_file_name`) and `file checksum` (`md5sum`), sorted by file name. A file
  missing either value keeps its row with that cell blank and a
  `====== WARNING ======` line naming it. Omitted when there are no
  RawMatrixFiles. RawMatrixFile fetches now always request `geo_file_name` and
  `md5sum`, whatever the schema snapshot lists.

- New `db2-flattener-download` command (`src/db2_flattener/download/`, also
  `python -m db2_flattener.download`): given a MatrixFileSet `--uuid`, downloads
  each RawMatrixFile from its `s3_uri` into `--out-directory`, saved as its
  `geo_file_name` when set and under its S3 file name otherwise. A file whose
  name already exists in the directory is not downloaded; a
  `====== WARNING ====== For {s3_uri}, there already exists a file
  {geo_file_name}` line is printed instead (`{geo_file_name}` only when set).
  RawMatrixFiles that would be saved under the same name (compared
  case-insensitively), such as Cell Ranger's `filtered_feature_bc_matrix.h5`
  from different prefixes or a repeated `geo_file_name`, are all refused with
  an error listing their `s3_uri`s; the other files still download. Exits 1 if
  any file is refused or fails to download (TOOLS-304).
- `get_report` now splits a report URL over 8000 characters by spreading the
  values of one filter key (the one with the most terms) across as many
  requests as needed, repeating every other filter term in each request, so
  the AND/OR meaning of the filter is unchanged. Previously it split only on
  `&@id=` and only in two: a long `&uuid=` filter produced one request with no
  filter at all (returning every object of that type), and `&@id=` filters
  gained an empty `&@id=` term. A filter that can't be split to fit (a single
  oversized term, or repeated terms already over the limit) now raises instead
  of being sent. A split request sends each value once and returns each object
  (by `@id`) once, even when the split key is a list field such as `samples`
  (TOOLS-304).
- SRA_FILE `instrument_model` now reads `UG 100` instead of
  `Ultima Genomics UG 100`, via the new `SRA_FILE_INSTRUMENT_MODEL_MAP`
  (separate from the GEO map). Other instrument models are unchanged.
- BIOHUB `samples_template_set{N}_v{M}.csv` gains an `intended_cell_types`
  column: the ControlledTerm `term_name`s from `CellLine.intended_cell_types`
  and `Organoid.intended_cell_types`, joined with `|` when there are several.
  Tissue and primary cell culture rows leave it empty.
- Refuse flatten (exit 1) when any gathered object has `status: deleted`,
  printing each object's type, UUID and raw aliases. Independent of audits and
  not overridable by `--force-*` flags.
- Every object fetch now always requests `@id`, `status` and `aliases`,
  whatever the schema-derived field list contains, ControlledTerms included.
- SRA outputs are now tab-separated: `SRA/sra_file_set{N}.tsv` and
  `SRA/sra_biosample_set{N}.tsv` replace the `.csv` files of the same names.
  These sheets routinely carry commas in free text (`*biomaterial_provider`,
  joined donor stages, probe-barcode lists), which forced quoting on most rows;
  tab-separated they need none.
- Tabs, carriage returns and newlines in SRA cells are collapsed to single
  spaces before writing (`sanitize_tabular_text` /
  `sanitize_tabular_dataframe` in `utils.py`), so the output parses under the
  naive tab-split the NCBI SRA/BioSample templates use. Other output files are
  unchanged.
- With `--s3-order-dir` / `-s`, BIOHUB `samples_template` and `guide_template`
  versions are resolved independently per type and set: each is the highest
  `v#` found in `{order}/metadata/` for that type and set, plus 1 (or `v1`
  when none).
- Only `{samples|guide}_template_set{N}_v{M}.csv` names count; names without a
  set number are ignored.
- `--biohub-version` / `-b` is unchanged: one explicit version for every
  template, mutually exclusive with `-s`.
- `DB2Flattener.flatten_matrix_file_set` takes `biohub_versions` (latest `v#`
  per `(kind, set)`) or `biohub_version` (one version for all templates).
- Removed internal helpers `resolve_biohub_version`,
  `resolve_biohub_version_from_s3`, `versions_from_metadata_names` and
  `next_biohub_version_from_versions` from `flatten/biohub_version.py`
  (replaced by `latest_biohub_versions_from_s3`,
  `latest_versions_from_metadata_names` and `next_biohub_version`).
- Write `SRA/sra_sequence_fileset_setN.tsv` alongside each `sra_file_setN.tsv`:
  one row per `library_ID` and SequenceFileSet uuid, repeating the library_ID
  when it has several sets. Tab-separated and sanitized like the other SRA
  sheets. Warn about libraries with no SequenceFileSet.
- Lower the SRA set size cap from 5 TB to 4.8 TB to leave a safety margin.
- Warn when a SequenceFileSet's file slots and its library's gathered
  `derived_from` SequenceFiles disagree: a slot file not in any `derived_from`,
  a slot file missing its `sequence_file_sets` back-link, a file listing a set
  with no slot for it, or a `derived_from` file that reaches no library.
- GEO `*SRA Experiment or Run` and `*BioSample` now come from the GEX library
  only, and both columns are always written (empty when no library has
  dbxrefs), since GEO requires them.
- Refuse flatten (exit 1), before any output is written, when a library's
  `library_construction_technology` is a multiome assay (`EFO:0030059`,
  `EFO:0920134`, `EFO:0920135`), or when any library's `feature_types`
  include ATAC. Each prints the affected GroupIDs first.

## 1.0.1

- Blank Control guide `intended_features_*` fields; warn when unexpected
  non-empty values are cleared.
- Normalize `control`/`CONTROL` to `Control`; warn on unexpected guide roles.
- Use `====WARNING====` for Control intended-feature and GEO title-length alerts.

## 1.0.0

- Refuse flatten when ERROR/WARNING object audits are present (strict by
  default); override with `--force-only-warning` / `--force-all-audit`.
- Flatten non-empty `audit` onto MAIN; skip empty audits.
- ExperimentalCondition fields and GEO/SRA title fallbacks (library,
  summary, pooled, GM strategy); warn when GEO titles exceed 120 characters.
- BIOHUB term IDs, multiplexed raw matrix file support, SRA file metadata
  updates, and BIOHUB/SRA/GEO output directory layout.

## 0.1.0

First versioned release. Flattens a MatrixFileSet to MAIN, BIOHUB, GEO,
SRA_BIOSAMPLE, SAMPLES, and GUIDE_METADATA CSVs.
