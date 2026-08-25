import re
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd

from db2_flattener.gather import lattice as DB2lattice
from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.schema.constants import (
    BIOHUB_SORT_ONTOLOGY_IDS,
    GENETIC_PERTURBATION_MAP,
    GUIDE_METADATA_COLUMNS,
    PROP_MAP_BIOHUB,
    PROP_MAP_GEO,
    PROP_MAP_SRA_BIOSAMPLE,
    REFORMAT_LIST,
    TISSUE_TYPE_MAP,
    Configs,
)
from db2_flattener.utils import (
    age_from_developmental_stage,
    collapse_dataframe,
    collapse_duplicate_columns,
    combine_bound_columns,
    extract_references_from_field,
    get_config_obj_type,
    get_url_prefix_from_id,
    is_empty,
    normalize_guide_rna_file_refs,
    sort_ontology_term_id_column,
    split_controlled_term_columns,
    strip_author_metadata_column_prefix,
    to_items,
)

# For SRA Biosample sheet creation, utilizing TISSUE_TYPE_MAP to get various sample types
SAMPLE_URL_PREFIXES = tuple(f"{name.replace(' ', '_')}s" for name in TISSUE_TYPE_MAP.values())


def age_with_units(value, units) -> str | None:
    """
    '29 years' from an age bound and its units. None when the bound is empty.

    Pluralises off the bare unit, so 'year' and 'years' in the data both give
    '29 years' and '1 year'.
    """
    if is_empty(value):
        return None
    # lower/upper_bound_age are typed number, so a null anywhere in the column
    # makes pandas store the lot as float
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    text = str(value).strip()
    if is_empty(units):
        return text
    unit = str(units).strip().rstrip("s")
    return f"{text} {unit}" if text == "1" else f"{text} {unit}s"


@dataclass(frozen=True)
class SRABiosampleDonorColumn:
    """
    One SRA_BIOSAMPLE column derived per donor.

    One source is given. 'columns' names MAIN columns directly and coalesces them,
    first non-null winning, for a field that sits on the donor under both a human
    and a non-human name. 'sample_fields' are fields on the sample, each expanded
    across SAMPLE_URL_PREFIXES and coalesced, then passed to 'transform' in order -
    two of them where a value needs its units. Neither means the column is the
    donor id itself.

    'render' is one of:
      labelled   'pooled: D1 - x, D2 - y', the shared donor enumeration
      positional 'x; y; not provided', one unlabelled entry per donor in that
                 same order, so entry n and Dn are the same donor
      summary    'pooled male and female', a summary over the donors, unordered

    An 'optional' column is dropped when no donor has a value, where a required
    one stays as a blank column.
    """

    name: str
    columns: tuple[str, ...] = ()
    sample_fields: tuple[str, ...] = ()
    render: str = "labelled"
    optional: bool = False
    transform: Callable[..., str | None] | None = None


class DB2Flattener:
    def __init__(self, connection: DB2lattice.Connection, configs: Configs):
        # Setup connection
        self.connection: DB2lattice.Connection = connection
        self.configs = configs

        # Initialize gatherer
        self.gatherer = DB2Gatherer(self.connection, configs)

    def flatten_matrix_file_set(self, matrix_file_set_uuid, output_prefix=None):
        """
        Flatten a MatrixFileSet into library-indexed data and save as CSV file
        """
        print(f"Processing MatrixFileSet {matrix_file_set_uuid}")

        # Gather all data
        complete_data = self.gatherer.gather_complete_library_data(matrix_file_set_uuid)

        if not complete_data:
            print("Error: No data gathered")
            return None

        print("Creating DataFrames...")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        prefix = output_prefix or f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}"
        main_output = f"{prefix}_MAIN.csv"
        biohub_output = f"{prefix}_BIOHUB.csv"
        geo_output = f"{prefix}_GEO.csv"
        sra_biosample_output = f"{prefix}_SRA_BIOSAMPLE.csv"
        sample_output = f"{prefix}_SAMPLES.csv"
        guide_output = f"{prefix}_GUIDE_METADATA.csv"

        # Create main DataFrame and sample DataFrame
        main_df, sample_df = self.create_dataframe(complete_data)
        main_df = main_df.dropna(axis=1, how="all")
        sample_df = sample_df.dropna(axis=1, how="all")

        # Split dict columns into _term_id / _term_name before writing CSV
        main_df = split_controlled_term_columns(main_df)
        if sample_df is not None and not sample_df.empty:
            sample_df = split_controlled_term_columns(sample_df)

        # Save main DataFrame to CSV
        print(f"Saving main DataFrame to {main_output}...")
        main_df.to_csv(main_output, index=False)

        print(f"✅ Main CSV file created: {main_output}")
        print(f"   Rows: {len(main_df)}")
        print(f"   Columns: {len(main_df.columns)}")

        # Create Biohub DataFrame from main and sample df
        biohub_df = self.create_biohub_dataframe(main_df)
        print(f"Saving biohub DataFrame to {biohub_output}...")
        biohub_df.to_csv(biohub_output, index=False)
        print(f"✅ Biohub CSV file created: {biohub_output}")

        # Create GEO DataFrame from main DataFrame
        geo_df = self.create_geo_dataframe(main_df)
        print(f"Saving geo DataFrame to {geo_output}...")
        geo_df.to_csv(geo_output, index=False)
        print(f"✅ GEO CSV file created: {geo_output}")

        # Create SRA/BioSample DataFrame from main DataFrame
        sra_biosample_df = self.create_sra_biosample_dataframe(main_df)
        print(f"Saving SRA/BioSample DataFrame to {sra_biosample_output}...")
        sra_biosample_df.to_csv(sra_biosample_output, index=False)
        print(f"✅ SRA/BioSample CSV file created: {sra_biosample_output}")
        print(f"   Rows: {len(sra_biosample_df)}")
        print(f"   Columns: {len(sra_biosample_df.columns)}")

        guide_file = self._resolve_guide_rna_file(complete_data)
        guide_df = self.create_guide_metadata_dataframe(guide_file)
        if guide_df is not None and not guide_df.empty:
            print(f"Saving guide metadata DataFrame to {guide_output}...")
            guide_df.to_csv(guide_output, index=False)
            print(f"✅ Guide metadata CSV file created: {guide_output}")
            print(f"   Rows: {len(guide_df)}")
            print(f"   Columns: {len(guide_df.columns)}")

        # Save sample DataFrame if it exists
        if sample_df is not None and not sample_df.empty:
            print(f"Saving sample DataFrame to {sample_output}...")
            sample_df.to_csv(sample_output, index=True)

            print(f"✅ Sample CSV file created: {sample_output}")
            print(f"   Rows: {len(sample_df)}")
            print(f"   Columns: {len(sample_df.columns)}")

        return main_output

    def create_dataframe(self, complete_data):
        """Create main library/raw-file DataFrame and sample DataFrame keyed by raw matrix file"""
        libraries_data = complete_data["libraries"]
        resolved_controlled_terms = complete_data["resolved_objects"].get("ControlledTerm", {})

        print("Creating DataFrame by raw matrix file...")

        rows = []
        sample_df = None
        new_sample_df = None

        # Raw matrix file-based rows - samples field is always present
        # First, collect all raw matrix files and their associated libraries
        raw_file_to_libraries = {}
        sample_metadata = {}  # One row per unique sample, keyed by sample_alias

        for lib_data in libraries_data.values():
            library = lib_data["library"]
            samples = lib_data["samples"]
            raw_matrix_files = lib_data["raw_matrix_files"]

            for raw_file in raw_matrix_files:
                raw_file_id = raw_file.get("@id")
                if raw_file_id not in raw_file_to_libraries:
                    raw_file_to_libraries[raw_file_id] = {"library_entries": [], "all_samples": []}

                # Add this library's data to the raw matrix file
                existing_lib_ids = {
                    entry["library"].get("@id")
                    for entry in raw_file_to_libraries[raw_file_id]["library_entries"]
                }
                if library.get("@id") not in existing_lib_ids:
                    raw_file_to_libraries[raw_file_id]["library_entries"].append(
                        {"library": library, "raw_file": raw_file}
                    )
                raw_file_to_libraries[raw_file_id]["all_samples"].extend(samples)

                # Collect per-sample metadata for later merge with raw matrix files
                for sample_ref in raw_file.get("samples", []):
                    sample_obj = next((s for s in samples if s.get("@id") == sample_ref), None)
                    if sample_obj:
                        sample_alias = self._get_clean_alias(sample_obj)
                        if sample_alias in sample_metadata:
                            continue

                        sample_metadata[sample_alias] = {}

                        sample_type = get_config_obj_type(sample_obj, self.configs)
                        sample_config = self.configs.OBJECT_CONFIG[sample_type]
                        sample_refs = sample_config.get("references", {})
                        for field in sample_config.get("fields", []):
                            value = sample_obj.get(field)
                            # Special handling for author_metadata dictionary
                            if field == "author_metadata" and isinstance(value, dict):
                                for key, val in value.items():
                                    field_name = f"{sample_type}_{field}_{'_'.join(key.split(' '))}"
                                    sample_metadata[sample_alias][field_name] = val
                            else:
                                if self._is_controlled_term_field(field, sample_refs):
                                    value = self._resolve_controlled_term_value(
                                        value, resolved_controlled_terms
                                    )
                                field_name = f"{sample_type}_{field}"
                                sample_metadata[sample_alias][field_name] = value

                        self._flatten_resolved_references(
                            sample_obj,
                            lib_data,
                            sample_metadata,
                            sample_alias,
                            resolved_controlled_terms,
                        )

        # Create one row per (raw matrix file, library)
        for file_data in raw_file_to_libraries.values():
            library_entries = file_data["library_entries"]
            samples = file_data["all_samples"]

            # For raw_matrix_file_alias and raw_file_samples, it doesn't matter which raw file copy we use
            # As they will all have the same value
            representative_raw_file = library_entries[0]["raw_file"]
            sample_aliases = []
            for sample_ref in representative_raw_file.get("samples", []):
                sample_obj = next((s for s in samples if s.get("@id") == sample_ref), None)
                if sample_obj:
                    sample_aliases.append(self._get_clean_alias(sample_obj))

            shared = {
                "raw_matrix_file_alias": self._get_clean_alias(representative_raw_file),
                "raw_file_samples": self._join_unique(sample_aliases),
            }

            for entry in library_entries:
                lib = entry["library"]
                # This libraries scoped copy of the raw matrix file
                # Which has the library specific sequence_file and sequence_file_set metadata
                raw_file = entry["raw_file"]
                row = dict(shared)
                lib_type = get_config_obj_type(lib, self.configs)
                for field in self.configs.OBJECT_CONFIG[lib_type].get("fields", []):
                    row[f"{lib_type}_{field}"] = lib.get(field)

                for obj_type in ("sequence_files", "sequence_file_sets", "tabular_files"):
                    obj_config = self.configs.OBJECT_CONFIG.get(obj_type, {})
                    for obj_field in obj_config.get("fields", []):
                        values = [
                            obj.get(obj_field)
                            for obj in raw_file.get(obj_type, [])
                            if obj.get(obj_field) not in (None, "")
                        ]
                        row[f"{obj_type}_{obj_field}"] = (
                            self._join_unique(values) if values else None
                        )
                rows.append(row)

        main_df = pd.DataFrame(rows)

        # Merge sample metadata with main DataFrame columns
        # Samples is one row per unique (raw matrix file + sample)
        sample_df = None
        if sample_metadata and not main_df.empty:
            per_sample_df = pd.DataFrame.from_dict(sample_metadata, orient="index")
            per_sample_df.index.name = "sample_alias"
            per_sample_df = per_sample_df.reset_index()

            rmf_sample_keys = main_df[
                ["raw_matrix_file_alias", "raw_file_samples"]
            ].drop_duplicates()

            sample_df = per_sample_df.merge(
                rmf_sample_keys, left_on="sample_alias", right_on="raw_file_samples", how="right"
            )

            new_sample_df = sample_df.set_index("raw_matrix_file_alias")

            main_df = main_df.merge(
                per_sample_df, left_on="raw_file_samples", right_on="sample_alias", how="left"
            )

            print(
                f"Creating sample DataFrame with {len(sample_df)} rows "
                f"({sample_df['raw_matrix_file_alias'].nunique()} raw matrix files)..."
            )

        return main_df, new_sample_df

    def _row_is_gex(self, row) -> bool:
        """Filter df to only GEX libraries"""
        droplet_ft = row.get("droplet_based_libraries_feature_types")
        plate_ft = row.get("plate_based_libraries_feature_types")

        def has_gex(ft):
            if ft is None or (isinstance(ft, float) and pd.isna(ft)):
                return None  # missing
            if isinstance(ft, str):
                return "Gene Expression" in ft
            if isinstance(ft, list):
                return "Gene Expression" in ft
            return "Gene Expression" in str(ft)

        droplet = has_gex(droplet_ft)
        plate = has_gex(plate_ft)

        if droplet is True or plate is True:
            return True
        if droplet is False or plate is False:
            return False

        # Missing feature_types: plate assumed GEX, droplet assumed non-GEX
        if pd.notna(row.get("plate_based_libraries_@id")) or pd.notna(
            row.get("plate_based_libraries_CRO_group_identifier")
        ):
            return True
        if pd.notna(row.get("droplet_based_libraries_@id")) or pd.notna(
            row.get("droplet_based_libraries_CRO_group_identifier")
        ):
            return False
        return True  # if no GEX found, keep all

    def create_geo_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build GEO submission dataframe from already-split main_df, taking GEX libraries only

        Expects _term_name columns (not raw dict columns).
        """
        gex_mask = main_df.apply(self._row_is_gex, axis=1)
        geo_source = main_df[gex_mask].copy()
        print(f"GEO: filtered to {len(geo_source)} GEX rows out of {len(main_df)} MAIN rows")

        subset_keys = [k for k in PROP_MAP_GEO if k in geo_source.columns]
        geo_df = geo_source[subset_keys].copy()
        geo_df.rename(columns=PROP_MAP_GEO, inplace=True)

        group_col = "*library name"
        return collapse_dataframe(geo_df, group_col=group_col)

    def _sample_probe_barcode_map(self, library_samples):
        """
        Build one library's 'sample_name : barcode|barcode' map.

        library_samples is a library's embedded 'samples' field, each item
        carrying 'aliases' and 'multiplexing_barcodes'. Entries are sorted by
        cleaned alias, so the string does not depend on the API's ordering.
        Returns None when no sample carries barcodes.
        """
        if not isinstance(library_samples, list):
            return None

        entries = []
        for sample in library_samples:
            if not isinstance(sample, dict):
                continue
            alias = self._get_clean_alias(sample)
            if not alias:
                continue
            barcodes = sample.get("multiplexing_barcodes") or []
            if isinstance(barcodes, str):
                barcodes = [barcodes]
            entries.append((alias, "|".join(str(barcode) for barcode in barcodes)))

        if not any(barcodes for _, barcodes in entries):
            return None

        entries.sort(key=lambda entry: entry[0])
        return ", ".join(f"{alias} : {barcodes}" for alias, barcodes in entries)

    def create_sra_biosample_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build the SRA/BioSample dataframe from main_df: one row per library.

        Grouped on the library alias, which PROP_MAP_SRA_BIOSAMPLE renames to
        'sample_name' - the submitted 'sample' is the sequencing library.
        """
        alias_cols = [
            col
            for col in ("droplet_based_libraries_aliases", "plate_based_libraries_aliases")
            if col in main_df.columns
        ]
        if not alias_cols:
            print("Warning: MAIN has no library aliases column; SRA_BIOSAMPLE will be empty")
            return pd.DataFrame()

        columns_to_keep = [k for k in PROP_MAP_SRA_BIOSAMPLE if k in main_df.columns]
        sra_df = main_df[columns_to_keep].copy()

        # Strip the lab prefix before the rename merges droplet and plate into one
        for alias_col in alias_cols:
            sra_df[alias_col] = sra_df[alias_col].map(self._clean_alias_cell)

        sra_df.rename(columns=PROP_MAP_SRA_BIOSAMPLE, inplace=True)
        sra_df = collapse_duplicate_columns(sra_df)

        # Built while sample_name is still row-aligned with main_df, attached after
        # the collapse - one value per library already, so nothing to aggregate
        donor = self._donor_cells_by_library(main_df, sra_df["sample_name"])
        provider = self._biomaterial_provider_by_library(main_df, sra_df["sample_name"])
        collection_date = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "date_obtained"
        )
        geo_loc_name = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "collection_geographical_location"
        )
        tissue = self._tissue_by_library(main_df, sra_df["sample_name"])
        suspension_type = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "suspension_type", optional=True
        )
        perturbation = self._perturbation_by_library(main_df, sra_df["sample_name"])
        perturbation_factors = self._perturbation_factors_by_library(main_df, sra_df["sample_name"])

        # Barcodes are per row, so they do go through the collapse: that is what
        # turns rows disagreeing on the map into a list
        barcode_map = None
        for samples_col in ("droplet_based_libraries_samples", "plate_based_libraries_samples"):
            if samples_col not in main_df.columns:
                continue
            mapped = main_df[samples_col].map(self._sample_probe_barcode_map)
            barcode_map = mapped if barcode_map is None else barcode_map.fillna(mapped)
        if barcode_map is not None:
            sra_df["sample_name: sample_probe_barcode"] = barcode_map

        unnamed = int(sra_df["sample_name"].isna().sum())
        if unnamed:
            print(
                f"Warning: dropping {unnamed} of {len(sra_df)} MAIN row(s) with no library "
                "alias from SRA_BIOSAMPLE"
            )
            sra_df = sra_df[sra_df["sample_name"].notna()]

        if sra_df.empty:
            return sra_df

        if len(sra_df.columns) == 1:
            # groupby().agg({}) raises, so a frame of nothing but the key dedupes
            sra_df = sra_df.drop_duplicates().reset_index(drop=True)
        else:
            sra_df = collapse_dataframe(sra_df, group_col="sample_name")

        # Optional donor columns are absent from 'donor' when nothing had a value
        for spec in self.DONOR_COLUMNS:
            if donor[spec.name]:
                sra_df[spec.name] = sra_df["sample_name"].map(donor[spec.name])
        if tissue:
            sra_df["*tissue"] = sra_df["sample_name"].map(tissue)
        if provider:
            sra_df["*biomaterial_provider"] = sra_df["sample_name"].map(provider)
        sra_df["*collection_date"] = sra_df["sample_name"].map(collection_date)
        sra_df["*geo_loc_name"] = sra_df["sample_name"].map(geo_loc_name)

        # Optional: left out entirely when no library has a value
        if suspension_type:
            sra_df["suspension_type"] = sra_df["sample_name"].map(suspension_type)
        if perturbation:
            sra_df["experimental_perturbation"] = sra_df["sample_name"].map(perturbation)
        if perturbation_factors:
            sra_df["experimental_perturbation_factors"] = sra_df["sample_name"].map(
                perturbation_factors
            )

        return sra_df

    @staticmethod
    def _donor_id_text(value) -> str:
        """Render a donor id as text without pandas' int-to-float artifacts."""
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        return str(value)

    @staticmethod
    def _donor_sort_key(donor_id: str):
        """Numeric ids sort numerically, so D10 does not land between D1 and D2."""
        return (0, int(donor_id), "") if donor_id.isdigit() else (1, 0, donor_id)

    @staticmethod
    def _format_pooled_donor_cell(labelled_values) -> str | None:
        """
        'pooled: D1 - x, D2 - y' across donors, or the bare value for a lone donor.

        labelled_values is [(label, value)] already in donor order. Entries with
        no value drop out, but the surviving labels keep their numbers so they
        still name the same donor as in the other donor column.
        """
        entries = [(label, value) for label, value in labelled_values if not is_empty(value)]
        if not entries:
            return None
        if len(labelled_values) == 1:
            return str(entries[0][1])
        return "pooled: " + ", ".join(f"{label} - {value}" for label, value in entries)

    # Anything not listed sorts alphabetically after these
    SEX_POOL_ORDER = ("male", "female")

    @classmethod
    def _format_pooled_sex(cls, sexes) -> str | None:
        """
        'male' for one sex, 'pooled male and female' across a mixed pool.

        Values outside the pair pool the same way: 'pooled male and unknown'.
        """
        unique = {str(sex).strip() for sex in sexes if not is_empty(sex)}
        if not unique:
            return None

        def sort_key(sex):
            if sex in cls.SEX_POOL_ORDER:
                return (0, cls.SEX_POOL_ORDER.index(sex), "")
            return (1, 0, sex)

        ordered = sorted(unique, key=sort_key)
        if len(ordered) == 1:
            return ordered[0]
        return "pooled " + ", ".join(ordered[:-1]) + " and " + ordered[-1]

    @staticmethod
    def _provider_titles(value) -> list[str]:
        """
        Provider titles from a 'sources' or 'lab' cell: a dict, or a list of them.

        A bare '@id' path yields nothing: 'sources' declares a multi-type linkTo,
        which extract.get_link_to() ignores, so the gatherer never resolves it to
        an object carrying a title.
        """
        titles = []
        for item in to_items(value):
            if isinstance(item, dict):
                title = item.get("title") or item.get("name")
                if title and str(title).strip():
                    titles.append(str(title).strip())
        return titles

    @staticmethod
    def _duration_text(lower, upper, units) -> str:
        """
        '8 hour' from a duration and its units, '8-24 hour' when the bounds differ.

        Units stay verbatim rather than pluralised - '8 hour stimulation' reads
        adjectivally, unlike the '29 years' the age bound columns want.
        """
        bounds = []
        for bound in (lower, upper):
            if is_empty(bound):
                continue
            # durations are typed number, so a null anywhere makes the column float
            if isinstance(bound, float) and bound.is_integer():
                bound = int(bound)
            bounds.append(str(bound).strip())
        if not bounds:
            return ""
        span = bounds[0] if len(set(bounds)) == 1 else "-".join(bounds)
        unit = "" if is_empty(units) else str(units).strip()
        return f"{span} {unit}".strip()

    @staticmethod
    def _format_pooled_values(values) -> str | None:
        """'pooled: a, b' across several values, or the bare value for a lone one."""
        entries = [value for value in values if value]
        if not entries:
            return None
        if len(entries) == 1:
            return entries[0]
        return "pooled: " + ", ".join(entries)

    def _perturbation_by_library(self, main_df, library_key):
        """
        Build the 'experimental_perturbation' cell per library, keyed by library alias.

        Reads the treatment's duration and description straight off 'treatments_*',
        which is already the object prefix and so needs no sample expansion. A
        sample with no treatment contributes 'not provided', which is why a library
        that only partly perturbed reads 'pooled: 8 hour stimulation, not provided'.
        Optional, so nothing anywhere drops the column.
        """
        columns = [
            main_df[name] if name in main_df.columns else None
            for name in (
                "treatments_lower_bound_duration",
                "treatments_upper_bound_duration",
                "treatments_duration_units",
                "treatments_description",
            )
        ]
        if all(column is None for column in columns):
            return {}

        empty = pd.Series(None, index=main_df.index, dtype=object)
        lower, upper, units, descriptions = (empty if c is None else c for c in columns)

        found: dict[str, set[str]] = {}
        missing: set[str] = set()
        for library, low, high, unit, description in zip(
            library_key, lower, upper, units, descriptions, strict=True
        ):
            if not isinstance(library, str):
                continue
            found.setdefault(library, set())
            duration = self._duration_text(low, high, unit)
            text = " ".join(
                part
                for part in (duration, "" if is_empty(description) else str(description).strip())
            ).strip()
            if text:
                found[library].add(text)
            else:
                missing.add(library)

        if not any(found.values()):
            return {}

        return {
            library: self._format_pooled_values(
                sorted(present) + (["not provided"] if library in missing else [])
            )
            for library, present in found.items()
        }

    def _perturbation_factors_by_library(self, main_df, library_key):
        """
        Build the 'experimental_perturbation_factors' cell per library, keyed by alias.

        Each sample contributes the set of ontological terms its treatments name,
        bracketed when there is more than one so a reader can tell which factors
        went together: '[anti-CD2_HUMAN, IL2_HUMAN]'. Distinct sets are pooled
        across the library, and a sample with no treatment contributes 'na'.
        Optional, so nothing anywhere drops the column.
        """
        terms = main_df.get("treatments_ontological_term_term_name")
        if terms is None:
            return {}

        found: dict[str, set[str]] = {}
        missing: set[str] = set()
        for library, value in zip(library_key, terms, strict=True):
            if not isinstance(library, str):
                continue
            found.setdefault(library, set())
            names = sorted({str(name).strip() for name in to_items(value)})
            if not names:
                missing.add(library)
            elif len(names) == 1:
                found[library].add(names[0])
            else:
                found[library].add(f"[{', '.join(names)}]")

        if not any(found.values()):
            return {}

        return {
            library: self._format_pooled_values(
                sorted(present) + (["na"] if library in missing else [])
            )
            for library, present in found.items()
        }

    def _biomaterial_provider_by_library(self, main_df, library_key):
        """
        Build the 'biomaterial_provider' cell per library, keyed by library alias.

        Falls back per row, not per column: a sample with no usable 'sources' uses
        its own 'lab' even where a sibling has one. Distinct titles join with '; '.
        """
        prefixes = SAMPLE_URL_PREFIXES
        sources = self._coalesce_columns(main_df, [f"{prefix}_sources" for prefix in prefixes])
        labs = self._coalesce_columns(main_df, [f"{prefix}_lab" for prefix in prefixes])
        if sources is None and labs is None:
            print(
                "Warning: MAIN has no sample sources or lab column; "
                "SRA_BIOSAMPLE omits biomaterial_provider"
            )
            return {}

        empty = pd.Series(None, index=main_df.index, dtype=object)
        sources = empty if sources is None else sources
        labs = empty if labs is None else labs

        titles_by_library: dict[str, set[str]] = {}
        for library, source, lab in zip(library_key, sources, labs, strict=True):
            if not isinstance(library, str):
                continue
            titles = self._provider_titles(source) or self._provider_titles(lab)
            if titles:
                titles_by_library.setdefault(library, set()).update(titles)

        return {library: "; ".join(sorted(titles)) for library, titles in titles_by_library.items()}

    @staticmethod
    def _format_positional_donor_cell(values) -> str | None:
        """
        'x; y; not provided' - one entry per donor, in donor order.

        values is one per donor, already in order. Nothing is deduplicated or
        sorted, so entry n always describes the same donor as Dn does in the
        labelled columns; a donor with no value holds its slot with the literal.
        None when no donor has a value, leaving the column out.
        """
        entries = list(values)
        if not any(entry for entry in entries):
            return None
        return "; ".join(entry if entry else "not provided" for entry in entries)

    def _sample_field_by_library(self, main_df, library_key, field, optional=False):
        """
        Per-library cell for a plain string sample field, keyed by library alias.

        A sample lacking a value contributes 'not provided' after any real ones,
        so a library whose samples disagree reads '2023-01-05; not provided'.

        'optional' decides only whether the column exists: nothing anywhere means
        nothing is returned and the column drops. Once one library has a value,
        every library gets a cell, filled the same way as a required column.
        """
        values = self._coalesce_columns(
            main_df, [f"{prefix}_{field}" for prefix in SAMPLE_URL_PREFIXES]
        )
        if values is None:
            values = pd.Series(None, index=main_df.index, dtype=object)

        found: dict[str, set[str]] = {}
        missing: set[str] = set()
        for library, value in zip(library_key, values, strict=True):
            if not isinstance(library, str):
                continue
            found.setdefault(library, set())
            if is_empty(value):
                missing.add(library)
            else:
                found[library].add(str(value).strip())

        if optional and not any(found.values()):
            return {}

        return {
            library: "; ".join(
                sorted(present) + (["not provided"] if library in missing or not present else [])
            )
            for library, present in found.items()
        }

    # Sample types with no tissue of origin to report
    TISSUELESS_SAMPLE_TYPES = ("cell_lines", "primary_cell_cultures")

    def _tissue_by_library(self, main_df, library_key):
        """
        Per-library 'tissue' cell, keyed by library alias.

        A cell line or primary cell culture has no tissue, so it reports
        'not available'. A tissue or organoid reports its sample_terms, which is
        an array - each term becomes its own entry. A sample with no term at all
        contributes nothing, which the schema should not allow anyway.
        """
        tissueless = [
            main_df[f"{prefix}_@id"]
            for prefix in self.TISSUELESS_SAMPLE_TYPES
            if f"{prefix}_@id" in main_df.columns
        ]
        terms = self._coalesce_columns(
            main_df,
            [
                f"{prefix}_sample_terms_term_name"
                for prefix in SAMPLE_URL_PREFIXES
                if prefix not in self.TISSUELESS_SAMPLE_TYPES
            ],
        )
        if terms is None and not tissueless:
            return {}

        empty = pd.Series(None, index=main_df.index, dtype=object)
        terms = empty if terms is None else terms
        # True where the row's sample is one of the tissueless types
        is_tissueless = [
            any(not is_empty(column.iloc[row]) for column in tissueless)
            for row in range(len(main_df))
        ]

        by_library: dict[str, set[str]] = {}
        for library, term, tissueless_row in zip(library_key, terms, is_tissueless, strict=True):
            if not isinstance(library, str):
                continue
            if tissueless_row:
                by_library.setdefault(library, set()).add("not available")
                continue
            for item in to_items(term):
                by_library.setdefault(library, set()).add(str(item).strip())

        return {library: "; ".join(sorted(values)) for library, values in by_library.items()}

    @staticmethod
    def _coalesce_columns(main_df, columns):
        """First non-null value across columns, or None when none of them exist."""
        out = None
        for col in columns:
            if col not in main_df.columns:
                continue
            out = main_df[col] if out is None else out.fillna(main_df[col])
        return out

    # Every per-donor column on the SRA Biosample sheet. Add one here and nothing else needs
    # touching: the D1..Dn enumeration, the empty-column rule and the assignment
    # all read off these specs.
    DONOR_COLUMNS = (
        SRABiosampleDonorColumn("*isolate"),
        SRABiosampleDonorColumn(
            "*age",
            sample_fields=("developmental_stages_term_name",),
            transform=age_from_developmental_stage,
        ),
        SRABiosampleDonorColumn(
            "*sex", columns=("human_donors_sex", "non_human_donors_sex"), render="summary"
        ),
        # ethnicity is on HumanDonor only, so a non-human run simply lacks it
        SRABiosampleDonorColumn(
            "ethnicity", columns=("human_donors_ethnicity_term_name",), optional=True
        ),
        SRABiosampleDonorColumn(
            "age_lower_bound",
            sample_fields=("lower_bound_age", "age_units"),
            render="positional",
            optional=True,
            transform=age_with_units,
        ),
        SRABiosampleDonorColumn(
            "age_upper_bound",
            sample_fields=("upper_bound_age", "age_units"),
            render="positional",
            optional=True,
            transform=age_with_units,
        ),
    )

    def _donor_cells_by_library(self, main_df, library_key):
        """
        Build the per-library donor cells, keyed by library alias then column name.

        Donors are enumerated D1..Dn in donor id order and every pooled cell reads
        from that one enumeration, so Dn names the same donor in each.
        """
        donor_ids = self._coalesce_columns(
            main_df, ("human_donors_cxg_donor_id", "non_human_donors_cxg_donor_id")
        )
        if donor_ids is None:
            print("Warning: MAIN has no donor id column; SRA_BIOSAMPLE omits the donor columns")
            return {spec.name: {} for spec in self.DONOR_COLUMNS}

        empty = pd.Series(None, index=main_df.index, dtype=object)

        def coalesced(candidates):
            found = self._coalesce_columns(main_df, candidates)
            return empty if found is None else found

        attributes = {}
        for spec in self.DONOR_COLUMNS:
            if spec.sample_fields:
                sources = [
                    coalesced([f"{p}_{field}" for p in SAMPLE_URL_PREFIXES])
                    for field in spec.sample_fields
                ]
            elif spec.columns:
                sources = [coalesced(spec.columns)]
            else:
                continue  # sourced from the donor id, seeded below
            if spec.transform:
                attributes[spec.name] = [spec.transform(*row) for row in zip(*sources, strict=True)]
            else:
                attributes[spec.name] = list(sources[0])

        # {library: {donor_id: {column: value}}} - first value per donor wins
        names = list(attributes)
        by_library: dict[str, dict[str, dict[str, str | None]]] = {}
        for library, donor_id, *values in zip(
            library_key, donor_ids, *attributes.values(), strict=True
        ):
            if not isinstance(library, str) or is_empty(donor_id):
                continue
            donors = by_library.setdefault(library, {})
            donor = donors.setdefault(self._donor_id_text(donor_id), dict.fromkeys(names, None))
            for name, value in zip(names, values, strict=True):
                if is_empty(donor[name]):
                    donor[name] = value

        cells: dict[str, dict[str, str | None]] = {spec.name: {} for spec in self.DONOR_COLUMNS}
        for library, donors in by_library.items():
            ordered = sorted(donors, key=self._donor_sort_key)
            labelled = [(f"D{number}", donor) for number, donor in enumerate(ordered, start=1)]
            for spec in self.DONOR_COLUMNS:
                if not spec.columns and not spec.sample_fields:
                    value = self._format_pooled_donor_cell(labelled)
                elif spec.render == "positional":
                    value = self._format_positional_donor_cell(
                        donors[donor][spec.name] for donor in ordered
                    )
                elif spec.render == "summary":
                    value = self._format_pooled_sex(donors[donor][spec.name] for donor in ordered)
                else:
                    value = self._format_pooled_donor_cell(
                        [(label, donors[donor][spec.name]) for label, donor in labelled]
                    )
                # An optional cell with no value is left out, so a column with no
                # data anywhere stays completely false and is skipped. A required one keeps
                # its None, so the column still appears blank when there are donors
                # This will allow us to catch missing required data
                if value is None and spec.optional:
                    continue
                cells[spec.name][library] = value

        return cells

    def _clean_alias_cell(self, value):
        """Strip the lab prefix from a raw aliases cell, which may be a list or a str."""
        if isinstance(value, list):
            return self._get_clean_alias({"aliases": value})
        if isinstance(value, str) and value:
            return self._get_clean_alias({"aliases": [value]})
        return value

    def create_guide_metadata_dataframe(self, file_info):
        """
        Build a guide-metadata DataFrame from the single shared guide RNA TabularFile.

        file_info comes from _resolve_guide_rna_file. None (no file, or more than
        one unique TabularFile) returns None.
        """
        if not file_info:
            return None

        guide_df = DB2lattice.read_tabular_file(file_info, self.connection)
        present = [col for col in GUIDE_METADATA_COLUMNS if col in guide_df.columns]
        if not present:
            print("Warning: guide RNA file has none of the expected GUIDE_METADATA columns")
            return None
        return guide_df[present].copy()

    def _resolve_guide_rna_file(self, complete_data):
        """Return the single gathered guide TabularFile, or None."""
        files = self._unique_guide_rna_files(complete_data)
        if not files:
            return None
        if len(files) > 1:
            labels = [
                file_info.get("@id") or file_info.get("s3_uri") or str(file_info)
                for file_info in files
            ]
            print(
                f"Warning: found {len(files)} unique guide RNA TabularFiles; "
                f"expected one. Skipping GUIDE_METADATA. Files: {labels}"
            )
            return None
        return files[0]

    @staticmethod
    def _unique_guide_rna_files(complete_data):
        """Collect unique TabularFile refs, preferring objects already gathered."""
        resolved = complete_data.get("resolved_objects", {})
        genetic_modifications = resolved.get("GeneticModification", {})
        gathered = resolved.get("TabularFile", {})
        unique = {}
        for gm in genetic_modifications.values():
            for file_info in normalize_guide_rna_file_refs(gm.get("guide_rna_files")):
                key = file_info.get("@id") or file_info.get("s3_uri")
                if not key:
                    continue
                unique.setdefault(key, gathered.get(file_info.get("@id"), file_info))
        return list(unique.values())

    def create_biohub_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build Biohub samples dataframe from main_df
        """
        columns_to_keep = [k for k in PROP_MAP_BIOHUB if k in main_df.columns]
        columns_to_keep.extend([k for k in main_df.columns if re.search("_author_metadata_", k)])
        biohub_df = main_df[columns_to_keep].copy()
        biohub_df.rename(columns=PROP_MAP_BIOHUB, inplace=True)
        biohub_df = strip_author_metadata_column_prefix(biohub_df)
        biohub_df = collapse_duplicate_columns(biohub_df)

        # Deduplicate rows, must join lists as they are not hashable
        list_cols = [
            col
            for col in biohub_df.columns
            if biohub_df[col].dropna().map(lambda x: isinstance(x, list)).any()
        ]
        for col in list_cols:
            biohub_df[col] = biohub_df[col].apply(lambda x: tuple(x) if isinstance(x, list) else x)
        biohub_df.drop_duplicates(inplace=True)

        # Add columns default values if not present
        if "disease" not in biohub_df.columns:
            biohub_df["disease"] = "normal"
        else:
            biohub_df["disease"] = biohub_df["disease"].apply(
                lambda v: pd.NA if (v is None or v == "" or v == [] or v == ()) else v
            )
            biohub_df["disease"] = biohub_df["disease"].fillna("normal")

        if "self_reported_ethnicity" not in biohub_df.columns:
            biohub_df["self_reported_ethnicity"] = np.where(
                biohub_df["organism"] == "Homo sapiens", "unknown", "na"
            )
        else:
            biohub_df.loc[
                biohub_df["self_reported_ethnicity"].isna()
                & (biohub_df["organism"] == "Homo sapiens"),
                "self_reported_ethnicity",
            ] = "unknown"
            biohub_df.loc[
                biohub_df["self_reported_ethnicity"].isna()
                & (biohub_df["organism"] != "Homo sapiens"),
                "self_reported_ethnicity",
            ] = "na"

        for col in BIOHUB_SORT_ONTOLOGY_IDS:
            if col in biohub_df.columns:
                biohub_df = sort_ontology_term_id_column(biohub_df, col)
        if "self_reported_ethnicity_ontology_term_id" in biohub_df.columns:
            biohub_df.drop(columns=["self_reported_ethnicity_ontology_term_id"], inplace=True)

        # Combine multiple columns into one
        biohub_df = combine_bound_columns(
            biohub_df,
            lower_col="experimental_conditions_lower_bound_duration",
            upper_col="experimental_conditions_upper_bound_duration",
            units_col="experimental_conditions_duration_units",
            out_col="experimental_perturbation_time_point",
        )
        biohub_df = combine_bound_columns(
            biohub_df,
            lower_col="tissues_lower_bound_age",
            upper_col="tissues_upper_bound_age",
            units_col="tissues_age_units",
            out_col="age",
        )

        # Update values to match schema
        biohub_df["tissue_type"] = biohub_df["tissue_type"].apply(
            lambda x: TISSUE_TYPE_MAP.get(x[0], x)
        )
        if "genetic_perturbation_strategy" in biohub_df.columns:
            biohub_df["genetic_perturbation_strategy"] = biohub_df[
                "genetic_perturbation_strategy"
            ].replace(GENETIC_PERTURBATION_MAP)
        for field in REFORMAT_LIST:
            if field in biohub_df.columns:
                biohub_df[field] = biohub_df[field].apply(
                    lambda x: "|".join(map(str, x)) if isinstance(x, (list, tuple)) else x
                )

        return biohub_df

    def _flatten_resolved_references(
        self, sample_obj, lib_data, sample_metadata, sample_alias, resolved_controlled_terms
    ):
        """Flatten resolved reference objects into sample_metadata columns."""
        sample_type = get_config_obj_type(sample_obj, self.configs)
        config = self.configs.OBJECT_CONFIG.get(sample_type, {})

        refs_by_prefix = {}
        for field_name, ref_types in config.get("references", {}).items():
            ref_type_list = [ref_types] if isinstance(ref_types, str) else ref_types
            if "controlled_terms" in ref_type_list:
                continue

            refs = extract_references_from_field(
                sample_obj.get(field_name), field_name, self.configs
            )
            for ref in refs:
                prefix = get_url_prefix_from_id(ref, self.configs)
                if prefix:
                    refs_by_prefix.setdefault(prefix, []).append(ref)

        for url_prefix, refs in refs_by_prefix.items():
            seen = set()
            unique_refs = [r for r in refs if not (r in seen or seen.add(r))]

            resolved_objs = [
                obj
                for ref in unique_refs
                for obj in lib_data.get(url_prefix, [])
                if obj.get("@id") == ref
            ]

            obj_config = self.configs.OBJECT_CONFIG.get(url_prefix, {})
            obj_refs = obj_config.get("references", {})

            for obj_field in obj_config.get("fields", []):
                col = f"{url_prefix}_{obj_field}"
                field_ref_types = obj_refs.get(obj_field, [])
                if isinstance(field_ref_types, str):
                    field_ref_types = [field_ref_types]
                is_controlled_term = "controlled_terms" in field_ref_types

                values = []
                for obj in resolved_objs:
                    value = obj.get(obj_field)
                    if value is None:
                        continue
                    if is_controlled_term:
                        value = self._resolve_controlled_term_value(
                            value, resolved_controlled_terms
                        )
                    if value not in (None, ""):
                        values.append(value)

                # Specific handling for author metadata fields
                if obj_field == "author_metadata":
                    values_by_key = {}
                    for value in values:
                        if isinstance(value, dict):
                            for key, val in value.items():
                                sanitized_key = "_".join(key.split(" "))
                                values_by_key.setdefault(sanitized_key, []).append(val)
                    for key, vals in values_by_key.items():
                        col = f"{url_prefix}_{obj_field}_{key}"
                        sample_metadata[sample_alias][col] = self._join_unique(vals)

                elif is_controlled_term:
                    # _join_unique() would str() these dicts into the cell
                    sample_metadata[sample_alias][col] = self._dedupe_terms(values)

                else:
                    sample_metadata[sample_alias][col] = (
                        self._join_unique(values) if values else None
                    )

    def _resolve_controlled_term(self, term_ref, resolved_controlled_terms):
        """
        Resolve a controlled term reference to the ControlledTerm object

        Returns the object ({'@id': ..., 'term_name': ...}) rather than the bare
        term id, so split_controlled_term_columns() can emit both a _term_id and
        a _term_name column. Some reference fields already arrive as embedded
        dicts carrying term_name; those are passed through untouched, so
        both reference shapes end up identical downstream.
        """
        if not term_ref:
            return None

        if isinstance(term_ref, dict):
            # Already embedded with a term_name - nothing to look up
            if term_ref.get("term_name"):
                return term_ref
            ref_id = term_ref.get("@id", "")
        else:
            ref_id = term_ref

        return resolved_controlled_terms.get(ref_id)

    @staticmethod
    def _is_controlled_term_field(field_name, references):
        """True if OBJECT_CONFIG marks this field as a controlled_terms reference"""
        ref_types = references.get(field_name, [])
        if isinstance(ref_types, str):
            ref_types = [ref_types]
        return "controlled_terms" in ref_types

    def _resolve_controlled_term_value(self, value, resolved_controlled_terms):
        """
        Resolve a controlled term field value to a dict, or list of dicts

        Deliberately not passed through _join_unique(): keeping dict shape is
        what lets split_controlled_term_columns() split the column into
        _term_id / _term_name. Array-typed fields keep list shape;
        split_term_cell() collapses a single-element list to scalars itself.
        """
        if value is None or value == []:
            return None

        refs = to_items(value)
        resolved = [
            term
            for term in (
                self._resolve_controlled_term(ref, resolved_controlled_terms) for ref in refs
            )
            if term
        ]

        if not resolved:
            return None
        return resolved if isinstance(value, list) else resolved[0]

    @staticmethod
    def _dedupe_terms(values):
        """
        Flatten and de-duplicate resolved controlled terms by @id,
        while keeping them as dicts
        """
        flat = []
        for value in values:
            flat.extend(to_items(value))

        by_id = {term["@id"]: term for term in flat if isinstance(term, dict) and term.get("@id")}
        terms = sorted(by_id.values(), key=lambda term: term["@id"])

        if not terms:
            return None
        return terms[0] if len(terms) == 1 else terms

    def _join_unique(self, items):
        """Join unique non-empty items with semicolon"""
        # Flatten any list-valued items (e.g. array-typed fields) before stringifying
        flattened_items = []
        for item in items:
            if isinstance(item, list):
                flattened_items.extend(item)
            else:
                flattened_items.append(item)

        # Filter out empty/None values. Do not use truthiness: False (and 0)
        # must be kept so boolean fields like is_pilot_order survive flattening.
        filtered_items = [
            str(item).strip()
            for item in flattened_items
            if item is not None and str(item).strip() != ""
        ]

        if not filtered_items:
            return

        # Remove duplicates and sort
        unique_items = sorted(set(filtered_items))
        return "; ".join(unique_items)

    def _get_clean_alias(self, obj_or_list):
        """Extract alias from an object or list of objects, cleaning the result"""
        if not obj_or_list:
            return

        aliases = []

        # Handle both single objects and lists
        if isinstance(obj_or_list, list):
            # List of objects
            for obj in obj_or_list:
                if isinstance(obj, dict):
                    obj_aliases = obj.get("aliases", [])
                    if obj_aliases:
                        aliases.extend(obj_aliases)
        else:
            # Single object
            if isinstance(obj_or_list, dict):
                obj_aliases = obj_or_list.get("aliases", [])
                if obj_aliases:
                    aliases.extend(obj_aliases)

        # If no aliases found, return
        if not aliases:
            return

        # Join unique aliases and clean
        result = self._join_unique(aliases)

        # Clean the result (remove prefix ending with ':')
        if ":" in result:
            return result.split(":", 1)[1]  # Take everything after the first ':'
        return result
