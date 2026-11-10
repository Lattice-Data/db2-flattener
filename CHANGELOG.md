# Changelog

## Unreleased

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
