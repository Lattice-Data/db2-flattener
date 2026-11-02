# Changelog

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
