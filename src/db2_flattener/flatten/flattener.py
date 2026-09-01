import re
from datetime import datetime

import numpy as np
import pandas as pd

from db2_flattener.gather import lattice as DB2lattice
from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.schema.constants import (
    BIOHUB_SORT_ONTOLOGY_IDS,
    GENETIC_PERTURBATION_MAP,
    GEO_EXPERIMENTAL_CONDITION_COLS,
    GEO_FLEX_LIBRARY_PROTOCOLS,
    GEO_INSTRUMENT_MODEL_MAP,
    GEO_LIBRARY_CARDINALITY_MAP,
    GEO_LIBRARY_STRATEGY_FEATURE_COL,
    GEO_LIBRARY_STRATEGY_MAP,
    GEO_LIBRARY_STRATEGY_PLATE_FEATURE_COL,
    GEO_LIBRARY_STRATEGY_SOURCE_COLS,
    GEO_SUSPENSION_TYPE_COLS,
    GEO_TITLE_TREATMENT_COLS,
    GEO_TREATMENT_COLS,
    GUIDE_METADATA_COLUMNS,
    PROP_MAP_BIOHUB,
    PROP_MAP_GEO,
    REFORMAT_LIST,
    TISSUE_TYPE_MAP,
    Configs,
)
from db2_flattener.utils import (
    collapse_dataframe,
    collapse_duplicate_columns,
    combine_bound_columns,
    expand_list_column,
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

    def _row_is_geo_library(self, row) -> bool:
        """Filter df to GEX or ATAC libraries."""
        droplet_ft = row.get("droplet_based_libraries_feature_types")
        plate_ft = row.get("plate_based_libraries_feature_types")

        def has_geo_feature(ft):
            if ft is None or (isinstance(ft, float) and pd.isna(ft)):
                return None  # missing
            if isinstance(ft, str):
                return "Gene Expression" in ft or "ATAC" in ft
            if isinstance(ft, list):
                return "Gene Expression" in ft or "ATAC" in ft
            text = str(ft)
            return "Gene Expression" in text or "ATAC" in text

        droplet = has_geo_feature(droplet_ft)
        plate = has_geo_feature(plate_ft)

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
        return True  # if no GEX/ATAC found, keep all

    def create_geo_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build GEO submission dataframe from already-split main_df, taking GEX and ATAC libraries.

        Expects _term_name columns (not raw dict columns).
        """
        geo_mask = main_df.apply(self._row_is_geo_library, axis=1)
        geo_source = main_df[geo_mask].copy()
        print(f"GEO: filtered to {len(geo_source)} GEX/ATAC rows out of {len(main_df)} MAIN rows")

        subset_keys = [k for k in PROP_MAP_GEO if k in geo_source.columns]
        subset_keys.extend(k for k in GEO_EXPERIMENTAL_CONDITION_COLS if k in geo_source.columns)
        subset_keys.extend(k for k in GEO_LIBRARY_STRATEGY_SOURCE_COLS if k in geo_source.columns)
        subset_keys.extend(k for k in GEO_TREATMENT_COLS if k in geo_source.columns)
        subset_keys.extend(
            k for k in GEO_TITLE_TREATMENT_COLS if k in geo_source.columns and k not in subset_keys
        )
        subset_keys.extend(k for k in geo_source.columns if re.search("_author_metadata_", k))
        geo_df = geo_source[subset_keys].copy()
        geo_df.rename(columns=PROP_MAP_GEO, inplace=True)
        geo_df = strip_author_metadata_column_prefix(geo_df)
        geo_df = collapse_duplicate_columns(geo_df)
        geo_df = self._summarize_geo_experimental_condition(geo_df)
        geo_df = self._summarize_geo_treatment(geo_df)
        geo_df = self._add_geo_molecule(geo_df)
        geo_df = self._add_geo_library_strategy(geo_df)

        group_col = "*library name"
        geo_df = collapse_dataframe(geo_df, group_col=group_col)
        col = "genetic_modifications_strategy"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].replace(GENETIC_PERTURBATION_MAP)
        col = "donor_sex"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].map(self._pool_mixed_geo_sex)
        geo_df = self._add_geo_title(geo_df)
        col = "single or paired-end"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].replace(GEO_LIBRARY_CARDINALITY_MAP)
        col = "*instrument model"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].replace(GEO_INSTRUMENT_MODEL_MAP)
        return expand_list_column(geo_df, "processed data file")

    @staticmethod
    def _pool_mixed_geo_sex(val):
        """Rewrite a list of male and female (either order) to a pooled label."""
        if isinstance(val, (list, tuple)):
            sexes = {str(item).strip().lower() for item in val if not is_empty(item)}
            if sexes == {"male", "female"}:
                return "pooled male and female"
        return val

    @staticmethod
    def _summarize_geo_experimental_condition(geo_df: pd.DataFrame) -> pd.DataFrame:
        """Build experimental_condition from source cols, then drop those sources."""
        source_cols = [c for c in GEO_EXPERIMENTAL_CONDITION_COLS if c in geo_df.columns]
        if not source_cols:
            return geo_df

        geo_df = combine_bound_columns(
            geo_df,
            lower_col="experimental_conditions_lower_bound_duration",
            upper_col="experimental_conditions_upper_bound_duration",
            units_col="experimental_conditions_duration_units",
            out_col="_exp_duration",
        )

        def _cell_text(val) -> str:
            if val is pd.NA:
                return ""
            try:
                if is_empty(val):
                    return ""
            except (TypeError, ValueError):
                return ""
            return str(val).strip()

        def _summarize_row(row: pd.Series):
            condition = _cell_text(row.get("experimental_conditions_condition"))
            text_value = _cell_text(row.get("experimental_conditions_text_value"))
            duration = _cell_text(row.get("_exp_duration"))
            right = " ".join(part for part in (text_value, duration) if part)
            if condition and right:
                return f"{condition}; {right}"
            return condition or right or pd.NA

        geo_df["experimental_condition"] = geo_df.apply(_summarize_row, axis=1)
        drop_cols = [
            c
            for c in (
                *GEO_EXPERIMENTAL_CONDITION_COLS,
                "_exp_duration",
            )
            if c in geo_df.columns
        ]
        return geo_df.drop(columns=drop_cols)

    @staticmethod
    def _summarize_geo_treatment(geo_df: pd.DataFrame) -> pd.DataFrame:
        """Build treatment from source cols. Empty rows become 'no treatment' if any row has a value."""
        source_cols = [c for c in GEO_TREATMENT_COLS if c in geo_df.columns]
        if not source_cols:
            return geo_df

        def _cell_text(val) -> str:
            if val is pd.NA:
                return ""
            try:
                if is_empty(val):
                    return ""
            except (TypeError, ValueError):
                return ""
            return str(val).strip()

        def _duration(row: pd.Series, *, include_upper: bool) -> str:
            lower = _cell_text(row.get("treatments_lower_bound_duration"))
            upper = _cell_text(row.get("treatments_upper_bound_duration")) if include_upper else ""
            if include_upper and lower and upper:
                bounds = lower if lower == upper else f"{lower}-{upper}"
            else:
                bounds = lower or (upper if include_upper else "")
            units = _cell_text(row.get("treatments_duration_units"))
            return " ".join(part for part in (bounds, units) if part)

        def _summarize_row(row: pd.Series):
            term = _cell_text(row.get("treatments_ontological_term_term_name"))
            description = _cell_text(row.get("treatments_description"))
            duration = _duration(row, include_upper=False)
            right = " ".join(part for part in (description, duration) if part)
            if term and right:
                return f"{term}; {right}"
            return term or right or None

        def _title_summarize_row(row: pd.Series):
            description = _cell_text(row.get("treatments_description"))
            duration = _duration(row, include_upper=True)
            return " ".join(part for part in (description, duration) if part) or None

        geo_df = geo_df.copy()
        geo_df["treatment"] = geo_df.apply(_summarize_row, axis=1)
        geo_df["_title_treatment"] = geo_df.apply(_title_summarize_row, axis=1)
        has_treatment = geo_df["treatment"].map(lambda v: bool(_cell_text(v)))
        if has_treatment.any():
            empty_treatment = ~has_treatment
            geo_df.loc[empty_treatment, "treatment"] = "no treatment"
            geo_df.loc[empty_treatment, "_title_treatment"] = "no treatment"
        return geo_df

    @staticmethod
    def _add_geo_title(geo_df: pd.DataFrame) -> pd.DataFrame:
        """Build title from library, strategy, treatment, and genetic modification fields."""

        def _cell_text(val) -> str:
            if val is pd.NA:
                return ""
            try:
                if is_empty(val):
                    return ""
            except (TypeError, ValueError):
                return ""
            return str(val).strip()

        def _is_pooled(samples) -> bool:
            if isinstance(samples, (list, tuple)):
                return len([item for item in samples if not is_empty(item)]) > 1
            text = _cell_text(samples)
            if not text:
                return False
            return len([part for part in text.split("; ") if part]) > 1

        def _treatment_segment(row: pd.Series) -> str:
            summarized = row.get("_title_treatment")
            if isinstance(summarized, (list, tuple)):
                items = [item for item in summarized if not is_empty(item)]
                return str(list(items)) if items else ""
            return _cell_text(summarized)

        def _summarize_row(row: pd.Series):
            left = " ".join(
                part
                for part in (
                    _cell_text(row.get("*library name")),
                    _cell_text(row.get("library_strategy")),
                )
                if part
            )
            groups = [left]
            if _is_pooled(row.get("samples")):
                groups.append("pooled")
            groups.extend(
                [
                    _treatment_segment(row),
                    _cell_text(row.get("genetic_modifications_strategy")),
                ]
            )
            return "; ".join(group for group in groups if group) or None

        geo_df = geo_df.copy()
        geo_df["title"] = geo_df.apply(_summarize_row, axis=1)
        drop_cols = [
            c
            for c in (*GEO_TREATMENT_COLS, *GEO_TITLE_TREATMENT_COLS, "_title_treatment")
            if c in geo_df.columns
        ]
        return geo_df.drop(columns=drop_cols) if drop_cols else geo_df

    @staticmethod
    def _add_geo_molecule(geo_df: pd.DataFrame) -> pd.DataFrame:
        """Build *molecule from feature_types and library_protocol."""

        def _has_feature(ft, name: str) -> bool:
            if is_empty(ft):
                return False
            if isinstance(ft, list):
                return name in ft
            return name in str(ft)

        def _map_row(row: pd.Series):
            ft = row.get(GEO_LIBRARY_STRATEGY_FEATURE_COL)
            if _has_feature(ft, "ATAC"):
                return "genomic DNA"
            if _has_feature(ft, "Gene Expression"):
                protocol = row.get("library_protocol")
                if protocol in GEO_FLEX_LIBRARY_PROTOCOLS:
                    return "total RNA"
                return "polyA RNA"
            return None

        geo_df = geo_df.copy()
        geo_df["*molecule"] = geo_df.apply(_map_row, axis=1)
        return geo_df

    @staticmethod
    def _add_geo_library_strategy(geo_df: pd.DataFrame) -> pd.DataFrame:
        """Build library_strategy from feature_types and suspension_type, then drop sources."""
        source_cols = [c for c in GEO_LIBRARY_STRATEGY_SOURCE_COLS if c in geo_df.columns]
        if not source_cols:
            return geo_df

        def _feature_types(row: pd.Series):
            droplet_ft = row.get(GEO_LIBRARY_STRATEGY_FEATURE_COL)
            if not is_empty(droplet_ft):
                return droplet_ft
            return row.get(GEO_LIBRARY_STRATEGY_PLATE_FEATURE_COL)

        def _has_feature(ft, name: str) -> bool:
            if is_empty(ft):
                return False
            if isinstance(ft, list):
                return name in ft
            return name in str(ft)

        def _suspension_types(row: pd.Series) -> list:
            items = []
            for col in GEO_SUSPENSION_TYPE_COLS:
                if col in row.index:
                    items.extend(to_items(row.get(col)))
            return list(dict.fromkeys(items))

        def _map_row(row: pd.Series):
            ft = _feature_types(row)
            if _has_feature(ft, "Gene Expression"):
                strategies = []
                for susp in _suspension_types(row):
                    mapped = GEO_LIBRARY_STRATEGY_MAP.get(("Gene Expression", susp))
                    if mapped:
                        strategies.append(mapped)
                unique = list(dict.fromkeys(strategies))
                if not unique:
                    return None
                return unique[0] if len(unique) == 1 else unique
            if _has_feature(ft, "ATAC"):
                return GEO_LIBRARY_STRATEGY_MAP[("ATAC", None)]
            return None

        geo_df = geo_df.copy()
        geo_df["library_strategy"] = geo_df.apply(_map_row, axis=1)
        return geo_df.drop(columns=source_cols)

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

        refs = value if isinstance(value, list) else [value]
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
            flat.extend(value if isinstance(value, list) else [value])

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
