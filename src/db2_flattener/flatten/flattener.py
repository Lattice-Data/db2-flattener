import re
from datetime import datetime
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd

from db2_flattener.gather import lattice as DB2lattice
from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.schema.constants import (
    BIOHUB_SORT_ONTOLOGY_IDS,
    CRO_ORDER_COL,
    EXPERIMENTAL_CONDITION_SUMMARY_COL,
    GENETIC_PERTURBATION_MAP,
    GEO_FLEX_LIBRARY_PROTOCOLS,
    GEO_INSTRUMENT_MODEL_MAP,
    GEO_LIBRARY_STRATEGY_FEATURE_COL,
    GEO_LIBRARY_STRATEGY_MAP,
    GEO_LIBRARY_STRATEGY_PLATE_FEATURE_COL,
    GEO_LIBRARY_STRATEGY_SOURCE_COLS,
    GEO_RUN_CARDINALITY_MAP,
    GEO_SUSPENSION_TYPE_COLS,
    GUIDE_METADATA_COLUMN_MAP,
    PROP_MAP_BIOHUB,
    PROP_MAP_GEO,
    PROP_MAP_SAMPLES,
    PROP_MAP_SRA_BIOSAMPLE,
    PROP_MAP_SRA_FILE,
    REFORMAT_LIST,
    SEX_ONTOLOGY_MAP,
    SRA_FILE_LIBRARY_LAYOUT_MAP,
    SRA_FILE_LIBRARY_SELECTION_MAP,
    SRA_FILE_LIBRARY_SOURCE_MAP,
    SRA_FILE_LIBRARY_STRATEGY_MAP,
    SRA_FILE_PLATFORM_PREFIXES,
    TISSUE_TYPE_MAP,
    Configs,
)
from db2_flattener.utils import (
    ages_from_developmental_stages,
    collapse_dataframe,
    collapse_duplicate_columns,
    combine_bound_columns,
    expand_list_column,
    extract_references_from_field,
    get_config_obj_type,
    get_url_prefix_from_id,
    is_empty,
    join_sequence_column,
    normalize_guide_rna_file_refs,
    numeric_text,
    sort_ontology_term_id_column,
    split_controlled_term_columns,
    strip_author_metadata_column_prefix,
    to_items,
)

# For SRA Biosample sheet creation, utilizing TISSUE_TYPE_MAP to get various sample types
SAMPLE_URL_PREFIXES = tuple(f"{name.replace(' ', '_')}s" for name in TISSUE_TYPE_MAP.values())

# SequenceFileSet slots that link to SequenceFiles, in filename column order
SEQUENCE_FILE_SET_FILE_SLOTS = (
    "index1",
    "index2",
    "read1",
    "read2",
    "read3",
    "trimmed_cram",
    "untrimmed_cram",
)

# Max total SequenceFile bytes per SRA file/biosample chunk (5 TB)
SRA_CHUNK_MAX_BYTES = 5 * 1000**4  # 5_000_000_000_000

# MAIN / biohub-source columns used to assign rows to a GeneticModification set
BIOHUB_LIBRARY_ID_COLS = (
    "droplet_based_libraries_@id",
    "plate_based_libraries_@id",
)


def age_with_units(value, units) -> str | None:
    """
    '29 years' from an age bound and its units. None when the bound is empty.

    Pluralises off the bare unit, so 'year' and 'years' in the data both give
    '29 years' and '1 year'.
    """
    if is_empty(value):
        return None
    text = numeric_text(value)
    if is_empty(units):
        return text
    unit = str(units).strip().rstrip("s")
    return f"{text} {unit}" if text == "1" else f"{text} {unit}s"


class DB2Flattener:
    def __init__(self, connection: DB2lattice.Connection, configs: Configs):
        # Setup connection
        self.connection: DB2lattice.Connection = connection
        self.configs = configs

        # Initialize gatherer
        self.gatherer = DB2Gatherer(self.connection, configs)

    def flatten_matrix_file_set(self, matrix_file_set_uuid, *, biohub_version: str):
        """
        Flatten a MatrixFileSet into library-indexed data and save CSVs under
        `{CRO_orders}_{timestamp}/` (or `MatrixFileSet_{uuid8}_{timestamp}/`
        when no CRO_order values are present).

        biohub_version is the resolved template suffix (e.g. ``v1``) used for
        BIOHUB samples_template / guide_template filenames.
        """
        print(f"Processing MatrixFileSet {matrix_file_set_uuid}")

        # Gather all data
        complete_data = self.gatherer.gather_complete_library_data(matrix_file_set_uuid)

        if not complete_data:
            print("Error: No data gathered")
            return None

        print("Creating DataFrames...")

        # Fail before creating output dirs when any GM lacks exactly one alias.
        self._validate_genetic_modification_aliases(complete_data)

        # Create main DataFrame and sample DataFrame
        main_df, sample_df = self.create_dataframe(complete_data)
        main_df = main_df.dropna(axis=1, how="all")
        if sample_df is not None:
            sample_df = sample_df.dropna(axis=1, how="all")

        # Split dict columns into _term_id / _term_name before writing CSV.
        # Biohub takes its own copy of the per-sample frame here, before
        # create_samples_dataframe pares that frame down for the SAMPLES sheet.
        main_df = split_controlled_term_columns(main_df)
        biohub_samples = None
        if sample_df is not None and not sample_df.empty:
            sample_df = split_controlled_term_columns(sample_df)
            biohub_samples = sample_df.copy()
            sample_df = self._attach_plate_library_id(main_df, sample_df)
            sample_df = self.create_samples_dataframe(sample_df)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        cro_orders = self._unique_cro_orders(main_df)
        if cro_orders:
            dir_name = f"{'_'.join(cro_orders)}_{timestamp}"
        else:
            dir_name = f"MatrixFileSet_{matrix_file_set_uuid[:8]}_{timestamp}"

        outdir = Path(dir_name)
        geo_dir = outdir / "GEO"
        sra_dir = outdir / "SRA"
        biohub_dir = outdir / "BIOHUB"
        for subdir in (geo_dir, sra_dir, biohub_dir):
            subdir.mkdir(parents=True, exist_ok=True)

        main_output = outdir / "MAIN.csv"
        geo_output = geo_dir / "geo_metadata.csv"
        sample_output = geo_dir / "demultiplexing_samples.csv"

        # Save main DataFrame to CSV
        print(f"Saving main DataFrame to {main_output}...")
        main_df.to_csv(main_output, index=False)

        print(f"✅ Main CSV file created: {main_output}")
        print(f"   Rows: {len(main_df)}")
        print(f"   Columns: {len(main_df.columns)}")

        # Multiplexed files contribute one Biohub row per sample. Other files
        # still contribute their main_df rows. Split by GeneticModification set.
        self._write_biohub_outputs(
            biohub_dir,
            geo_dir,
            main_df,
            biohub_samples,
            complete_data,
            biohub_version=biohub_version,
        )

        # Create GEO DataFrame from main DataFrame
        geo_df = self.create_geo_dataframe(main_df)
        print(f"Saving geo DataFrame to {geo_output}...")
        geo_df.to_csv(geo_output, index=False)
        print(f"✅ GEO CSV file created: {geo_output}")

        # SRA sheets: chunk by sample_name so each set stays under 5 TB of files
        sra_biosample_df = self.create_sra_biosample_dataframe(main_df)
        sequence_file_set_df = self.create_sequence_file_set_dataframe(complete_data)
        sra_files_df = self.create_sra_files_dataframe(
            main_df, sequence_file_set_df=sequence_file_set_df
        )
        self._write_sra_chunked_outputs(sra_dir, sra_files_df, sra_biosample_df)

        # Save sample DataFrame if it exists
        if sample_df is not None and not sample_df.empty:
            print(f"Saving sample DataFrame to {sample_output}...")
            sample_df.to_csv(sample_output, index=False)

            print(f"✅ Sample CSV file created: {sample_output}")
            print(f"   Rows: {len(sample_df)}")
            print(f"   Columns: {len(sample_df.columns)}")

        print(f"✅ Output directory created: {outdir}")
        return str(outdir)

    @staticmethod
    def _unique_cro_orders(main_df):
        """Sorted unique CRO_order values from MAIN ``sequence_file_sets_CRO_order``."""
        col = "sequence_file_sets_CRO_order"
        if col not in main_df.columns:
            return []
        orders = []
        for val in main_df[col]:
            if is_empty(val):
                continue
            for part in str(val).split("; "):
                part = part.strip()
                if part:
                    orders.append(part)
        return sorted(set(orders))

    def create_dataframe(self, complete_data):
        """Create main library/raw-file DataFrame and sample DataFrame keyed by raw matrix file"""
        libraries_data = complete_data["libraries"]
        resolved_controlled_terms = complete_data["resolved_objects"].get("ControlledTerm", {})

        print("Creating DataFrame by raw matrix file...")

        rows = []
        sample_df = None
        # One (raw matrix file, sample) pair per sample on each file
        file_sample_pairs = []

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
                                if field == "sources":
                                    titles = self._source_titles(value)
                                    if titles:
                                        sample_metadata[sample_alias][
                                            f"{sample_type}_sources_title"
                                        ] = self._join_unique(titles)

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

            # File-level fields are the same on every library copy of this raw file.
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
            # Direct copies, not _join_unique(): that stringifies booleans and integers.
            for field in (
                "is_multiplexed",
                "file_format",
                "file_size",
                "software",
                "software_version",
            ):
                shared[f"raw_matrix_files_{field}"] = representative_raw_file.get(field)

            # Recorded here, where the aliases are still separate values, rather
            # than by splitting raw_file_samples back apart on '; ' afterwards
            file_sample_pairs.extend(
                {
                    "raw_matrix_file_alias": shared["raw_matrix_file_alias"],
                    "sample_alias": alias,
                }
                for alias in sample_aliases
            )

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

        # SAMPLES is one row per (raw matrix file, sample), which is what makes a
        # multiplexed file readable. MAIN stays one row per (raw matrix file,
        # library), so it takes a copy of that frame rolled up to one row per
        # file: joining the exploded frame straight on would multiply every
        # library row by the file's sample count.
        sample_df = None
        if sample_metadata and file_sample_pairs and not main_df.empty:
            per_sample_df = pd.DataFrame.from_dict(sample_metadata, orient="index")
            per_sample_df.index.name = "sample_alias"
            per_sample_df = per_sample_df.reset_index()

            pairs_df = pd.DataFrame(file_sample_pairs).drop_duplicates()
            sample_df = pairs_df.merge(per_sample_df, on="sample_alias", how="left")

            rolled_up = sample_df.groupby("raw_matrix_file_alias", as_index=False).agg(
                {
                    column: self._combine_sample_values
                    for column in sample_df.columns
                    if column != "raw_matrix_file_alias"
                }
            )
            main_df = main_df.merge(rolled_up, on="raw_matrix_file_alias", how="left")

            print(
                f"Creating sample DataFrame with {len(sample_df)} rows "
                f"({sample_df['raw_matrix_file_alias'].nunique()} raw matrix files)..."
            )

        return main_df, sample_df

    def create_samples_dataframe(self, sample_df) -> pd.DataFrame:
        """
        Build the SAMPLES sheet from post-merge sample_df: keep mapped columns,
        and rename them.
        """
        sample_df = sample_df.copy()
        plate_experiment = self._frame_has_plate_library(sample_df)

        columns_to_keep = [k for k in PROP_MAP_SAMPLES if k in sample_df.columns]
        columns_to_keep.extend(
            c
            for c in ("human_donors_ethnicity_term_name", "tissues_@id")
            if c in sample_df.columns and c not in columns_to_keep
        )
        columns_to_keep.extend(
            c
            for c in sample_df.columns
            if re.search("_author_metadata_", c) and c not in columns_to_keep
        )
        if not plate_experiment:
            rt_sources = {key for key, name in PROP_MAP_SAMPLES.items() if name == "RT_index"}
            columns_to_keep = [c for c in columns_to_keep if c not in rt_sources]
        sample_df = sample_df[columns_to_keep].copy()

        sample_df.rename(columns=PROP_MAP_SAMPLES, inplace=True)
        sample_df = strip_author_metadata_column_prefix(sample_df)
        sample_df = collapse_duplicate_columns(sample_df)

        col = "genetic_modifications_strategy"
        if col in sample_df.columns:
            sample_df[col] = sample_df[col].replace(GENETIC_PERTURBATION_MAP)

        # An ExperimentalCondition can carry several controlled terms
        sample_df = join_sequence_column(sample_df, "experimental_condition_term", sep="; ")

        sample_df = self._add_samples_donor_ethnicity(sample_df)
        sample_df = join_sequence_column(sample_df, "sample_probe_barcode")
        sample_df = join_sequence_column(sample_df, "RT_index", sep=",")
        if plate_experiment and "sample_probe_barcode" in sample_df.columns:
            sample_df = sample_df.rename(columns={"sample_probe_barcode": "hash_index"})
        return sample_df

    def _attach_plate_library_id(self, main_df, sample_df):
        """Copy the plate library id onto sample rows so the SAMPLES sheet can see it.

        Rolled up to one value per raw matrix file, so two libraries do not
        repeat each sample. create_samples_dataframe drops the column after
        using it to choose the barcode header.
        """
        plate_id = "plate_based_libraries_@id"
        if (
            sample_df is None
            or sample_df.empty
            or plate_id not in main_df.columns
            or plate_id in sample_df.columns
            or "raw_matrix_file_alias" not in main_df.columns
            or "raw_matrix_file_alias" not in sample_df.columns
        ):
            return sample_df

        file_level = main_df.loc[:, ["raw_matrix_file_alias", plate_id]]
        file_level = file_level.groupby("raw_matrix_file_alias", as_index=False).agg(
            {plate_id: self._combine_sample_values}
        )
        return sample_df.merge(file_level, on="raw_matrix_file_alias", how="left")

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
        subset_keys.extend(k for k in GEO_LIBRARY_STRATEGY_SOURCE_COLS if k in geo_source.columns)
        if CRO_ORDER_COL in geo_source.columns:
            subset_keys.append(CRO_ORDER_COL)
        subset_keys.extend(k for k in geo_source.columns if re.search("_author_metadata_", k))
        geo_df = geo_source[subset_keys].copy()
        geo_df.rename(columns=PROP_MAP_GEO, inplace=True)
        geo_df = strip_author_metadata_column_prefix(geo_df)
        geo_df = collapse_duplicate_columns(geo_df)
        if "genetic_modification_alias" in geo_df.columns:
            geo_df["genetic_modification_alias"] = geo_df["genetic_modification_alias"].map(
                self._map_geo_genetic_modification_alias
            )
        geo_df = self._add_geo_molecule(geo_df)
        geo_df = self._add_geo_library_strategy(geo_df)

        ethnicity = {}
        if "*library name" in geo_df.columns:
            ethnicity = self._sample_field_by_library(
                geo_source,
                geo_df["*library name"],
                "ethnicity_term_name",
                prefixes=("human_donors",),
                optional=True,
                split_joined=True,
                formatter=self._format_pooled_values,
                include=self._tissue_row_mask(geo_source),
            )

        group_col = "*library name"
        geo_df = collapse_dataframe(geo_df, group_col=group_col)
        if ethnicity:
            geo_df["donor_ethnicity"] = geo_df["*library name"].map(ethnicity)
        col = "genetic_modifications_strategy"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].map(self._map_geo_genetic_strategy)
        col = "donor_sex"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].map(self._pool_mixed_geo_sex)
        col = "*library strategy"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].map(self._pool_mixed_geo_strategy)
        col = PROP_MAP_GEO[EXPERIMENTAL_CONDITION_SUMMARY_COL]
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].map(self._join_distinct_summaries)
        geo_df = self._add_geo_title(geo_df, source_col=col, drop_source=False)
        col = "*single or paired-end"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].replace(GEO_RUN_CARDINALITY_MAP)
        col = "*instrument model"
        if col in geo_df.columns:
            geo_df[col] = geo_df[col].replace(GEO_INSTRUMENT_MODEL_MAP)
        col = "*SRA Experiment or Run"
        if col in geo_df.columns:
            # Both columns come from the same dbxrefs cell, so read BioSample out
            # before the SRA accessions overwrite it.
            biosample = geo_df[col].map(self._extract_geo_biosample)
            geo_df[col] = geo_df[col].map(self._extract_geo_sra_accession)
            geo_df.insert(geo_df.columns.get_loc(col) + 1, "*BioSample", biosample)
        return expand_list_column(geo_df, "processed data file")

    @staticmethod
    def _extract_dbxref_accessions(val, prefix: str):
        """Take unique accessions after a given prefix from a dbxref cell."""
        accessions = []
        for item in to_items(val):
            text = str(item).strip()
            if text.startswith(prefix):
                accession = text[len(prefix) :]
                if accession:
                    accessions.append(accession)
        unique = list(dict.fromkeys(accessions))
        if not unique:
            return pd.NA
        return unique[0] if len(unique) == 1 else unique

    @classmethod
    def _extract_geo_sra_accession(cls, val):
        """Take unique accessions after an SRA: prefix from a dbxref cell."""
        return cls._extract_dbxref_accessions(val, "SRA:")

    @classmethod
    def _extract_geo_biosample(cls, val):
        """Take unique accessions after a Biomaterial: prefix from a dbxref cell."""
        return cls._extract_dbxref_accessions(val, "Biomaterial:")

    @staticmethod
    def _map_geo_genetic_strategy(val):
        """Translate each '; '-joined genetic strategy; leave non-strings unchanged."""
        if not isinstance(val, str):
            return val
        parts = [part.strip() for part in val.split("; ") if part.strip()]
        if not parts:
            return val
        mapped = [GENETIC_PERTURBATION_MAP.get(part, part) for part in parts]
        return "; ".join(mapped)

    @classmethod
    def _map_geo_genetic_modification_alias(cls, val):
        """Sanitize a GEO genetic_modification_alias cell; empty stays empty."""
        if is_empty(val):
            return val
        return cls._genetic_modification_alias_label(val)

    @staticmethod
    def _pool_mixed_geo_strategy(val):
        """Rewrite a mix of scRNA-seq, snRNA-seq, and RNA-Seq to RNA-Seq."""
        rna_strategies = {"scRNA-seq", "snRNA-seq", "RNA-Seq"}
        if isinstance(val, (list, tuple)):
            strategies = {str(item).strip() for item in val if not is_empty(item)}
            if len(strategies) >= 2 and strategies <= rna_strategies:
                return "RNA-Seq"
        return val

    @staticmethod
    def _pool_mixed_geo_sex(val):
        """Rewrite male and female, as a list or a '; '-joined string, to a pooled label."""
        if isinstance(val, (list, tuple)):
            parts = [part for item in val if not is_empty(item) for part in str(item).split("; ")]
        elif isinstance(val, str):
            parts = val.split("; ")
        else:
            return val
        sexes = {part.strip().lower() for part in parts if part.strip()}
        if sexes == {"male", "female"}:
            return "pooled male and female"
        return val

    @staticmethod
    def _add_geo_title(
        geo_df: pd.DataFrame,
        *,
        title_col: str = "*title",
        source_col: str = EXPERIMENTAL_CONDITION_SUMMARY_COL,
        drop_source: bool = True,
        name_col: str = "*library name",
        strategy_col: str = "*library strategy",
    ) -> pd.DataFrame:
        """
        Build title as '<CRO order>; <library name> <library strategy>; <summary>;
        pooled; <genetic modification strategy>'.

        The library name leads every title, which keeps GEO titles unique
        (GEO requires that) and the format the same whether or not a library
        has an ExperimentalCondition. Empty parts are skipped; see _geo_title.
        drop_source=False keeps the summary column, for a sheet that also
        reports it under its own name. The CRO order column is always dropped,
        since neither sheet reports it on its own.
        """
        geo_df = geo_df.copy()
        geo_df[title_col] = geo_df.apply(
            DB2Flattener._geo_title,
            axis=1,
            summary_col=source_col,
            name_col=name_col,
            strategy_col=strategy_col,
        )
        drop_cols = [CRO_ORDER_COL] if CRO_ORDER_COL in geo_df.columns else []
        if drop_source and source_col in geo_df.columns:
            drop_cols.append(source_col)
        return geo_df.drop(columns=drop_cols) if drop_cols else geo_df

    @staticmethod
    def _geo_title(
        row: pd.Series,
        summary_col: str = EXPERIMENTAL_CONDITION_SUMMARY_COL,
        name_col: str = "*library name",
        strategy_col: str = "*library strategy",
    ):
        """
        One row's title; see _add_geo_title for the order.

        A library with several CRO orders gets them sorted and joined with '; '.
        'pooled' only when the row has more than one distinct sample. Any cell
        may be a list after collapse_dataframe, or a '; '-joined string, so each
        is read through _distinct_parts. A library whose samples carry several
        distinct summaries gets them sorted and joined with '; '. Genetic
        modification strategies are mapped through GENETIC_PERTURBATION_MAP,
        since _map_geo_genetic_strategy leaves list cells unmapped.
        """
        parts = DB2Flattener._distinct_parts
        strategies = [
            GENETIC_PERTURBATION_MAP.get(part, part)
            for part in parts(row.get("genetic_modifications_strategy"))
        ]

        left = " ".join(
            part
            for part in (
                "; ".join(parts(row.get(name_col))),
                "; ".join(parts(row.get(strategy_col))),
            )
            if part
        )
        groups = [
            "; ".join(parts(row.get(CRO_ORDER_COL))),
            left,
            "; ".join(parts(row.get(summary_col))),
            "pooled" if len(parts(row.get("samples"))) > 1 else "",
            "; ".join(dict.fromkeys(strategies)),
        ]
        return "; ".join(group for group in groups if group) or None

    @staticmethod
    def _distinct_parts(value) -> list[str]:
        """
        Distinct non-empty items of a cell, sorted.

        A cell may be a scalar, a list (after collapse_dataframe), or a
        '; '-joined string (how create_dataframe joins several values), so
        every item is split on '; ' first.
        """
        items = {
            part.strip()
            for item in to_items(value)
            if not is_empty(item)
            for part in str(item).split("; ")
        } - {""}
        return sorted(items)

    @staticmethod
    def _join_distinct_summaries(value):
        """
        Distinct summaries in a cell, sorted and joined with '; '.

        Each item is split on '; ' first, because create_dataframe joins a
        sample's several conditions that way: ['a; b', 'a'] becomes 'a; b'.
        """
        return "; ".join(DB2Flattener._distinct_parts(value)) or None

    @staticmethod
    def _add_geo_molecule(geo_df: pd.DataFrame) -> pd.DataFrame:
        """Build *molecule from feature_types and library_protocol."""

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

        def _map_row(row: pd.Series):
            ft = _feature_types(row)
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
        """Build *library strategy from feature_types and suspension_type, then drop sources."""
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
        geo_df["*library strategy"] = geo_df.apply(_map_row, axis=1)
        return geo_df.drop(columns=source_cols)

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
        Build the SRA/BioSample dataframe from main_df: one row per CRO group.

        Grouped on the library's CRO_group_identifier, which
        PROP_MAP_SRA_BIOSAMPLE renames to 'sample_name'. That identifier spans a
        group's libraries, so the GEX libraries are taken and the rest dropped,
        as on GEO, one row per group rather than per library.
        """
        group_cols = [
            col
            for col in (
                "droplet_based_libraries_CRO_group_identifier",
                "plate_based_libraries_CRO_group_identifier",
            )
            if col in main_df.columns
        ]
        if not group_cols:
            print(
                "Warning: MAIN has no library CRO group identifier column; "
                "SRA_BIOSAMPLE will be empty"
            )
            return pd.DataFrame()

        if main_df.empty:
            return pd.DataFrame()

        # Filtered before anything is built: every cell below is keyed by a
        # sample_name taken from this frame, so the two have to stay row-aligned
        gex_mask = main_df.apply(self._row_is_gex, axis=1)
        print(
            f"SRA_BIOSAMPLE: filtered to {int(gex_mask.sum())} GEX rows "
            f"out of {len(main_df)} MAIN rows"
        )
        # The one cell taken from the unfiltered frame, so it has to be built here
        library_id = self._library_ids_by_group(main_df, gex_mask)
        main_df = main_df[gex_mask].reset_index(drop=True)
        if main_df.empty:
            return pd.DataFrame()

        columns_to_keep = [k for k in PROP_MAP_SRA_BIOSAMPLE if k in main_df.columns]
        sra_df = main_df[columns_to_keep].copy()

        sra_df.rename(columns=PROP_MAP_SRA_BIOSAMPLE, inplace=True)
        sra_df = collapse_duplicate_columns(sra_df)

        # Every cell below is built while sample_name is still row-aligned with
        # main_df, then attached after the collapse: one value per library already,
        # so there is nothing left to aggregate. Donor cells split their MAIN value
        # on '; ' because create_dataframe collapses a multi-donor sample that way.
        donor_prefixes = ("human_donors", "non_human_donors")
        isolate = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "cxg_donor_id",
            prefixes=donor_prefixes,
            split_joined=True,
            formatter=self._format_pooled_values,
        )
        sex = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "sex",
            prefixes=donor_prefixes,
            split_joined=True,
            formatter=self._format_pooled_sex,
        )
        ethnicity = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "ethnicity_term_name",
            prefixes=("human_donors",),
            optional=True,
            split_joined=True,
            formatter=self._format_pooled_values,
            include=self._tissue_row_mask(main_df),
        )
        age = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "developmental_stages_term_name",
            transform=ages_from_developmental_stages,
            formatter=self._format_pooled_values,
        )
        age_lower = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "lower_bound_age",
            "age_units",
            optional=True,
            transform=age_with_units,
            formatter=self._format_pooled_values,
        )
        age_upper = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "upper_bound_age",
            "age_units",
            optional=True,
            transform=age_with_units,
            formatter=self._format_pooled_values,
        )
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
        # Distinct ExperimentalCondition summaries across the library's samples.
        # A sample with no condition contributes nothing. Pooled with '; '
        # rather than ', ', since summaries themselves contain commas.
        condition_summary = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "summary",
            optional=True,
            prefixes=("experimental_conditions",),
            gap=None,
            split_joined=True,
            formatter=lambda values: self._format_pooled_values(values, sep="; "),
        )
        # intended_cell_types is only on cell lines and organoids, so a tissue or
        # primary cell culture genuinely cannot have one, hence 'not applicable'
        intended_cell_type = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "intended_cell_types_term_name",
            optional=True,
            prefixes=("cell_lines", "organoids"),
            gap="not applicable",
        )
        # The suspension columns below are on all four sample types, and take
        # their sources from BIOHUB: kits to selection_kits, enrichment factors
        # to selection_markers.
        enriched_cell_types = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "enriched_cell_types_term_name", optional=True
        )
        depleted_cell_types = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "depleted_cell_types_term_name", optional=True
        )
        selection_kits = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "selection_kits", optional=True
        )
        enrichment_factors = self._sample_field_by_library(
            main_df, sra_df["sample_name"], "selection_markers", optional=True
        )
        # preservation_method is on tissues alone, so any other sample type has
        # none to report, hence 'not applicable'
        preservation_method = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "preservation_method",
            optional=True,
            prefixes=("tissues",),
            gap="not applicable",
        )
        # 'strategy' lives on the linked GeneticModification, whose url_prefix is
        # itself the MAIN column prefix. A sample with no modification has no
        # strategy to report, hence 'not applicable'.
        genetic_strategy = self._sample_field_by_library(
            main_df,
            sra_df["sample_name"],
            "strategy",
            optional=True,
            prefixes=("genetic_modifications",),
            gap="not applicable",
            translate=GENETIC_PERTURBATION_MAP,
            split_joined=True,
        )

        # Barcodes are per row, so they do go through the collapse: that is what
        # turns rows disagreeing on the map into a list
        barcode_map = None
        for samples_col in ("droplet_based_libraries_samples", "plate_based_libraries_samples"):
            if samples_col not in main_df.columns:
                continue
            mapped = main_df[samples_col].map(self._sample_probe_barcode_map)
            barcode_map = mapped if barcode_map is None else barcode_map.fillna(mapped)
        if barcode_map is not None and barcode_map.map(lambda v: not is_empty(v)).any():
            barcode_column = (
                "sample_name: hash_index"
                if self._frame_has_plate_library(main_df)
                else "sample_name: sample_probe_barcode"
            )
            sra_df[barcode_column] = barcode_map

        unnamed = int(sra_df["sample_name"].isna().sum())
        if unnamed:
            print(
                f"Warning: dropping {unnamed} of {len(sra_df)} MAIN row(s) with no library "
                "CRO group identifier from SRA_BIOSAMPLE"
            )
            sra_df = sra_df[sra_df["sample_name"].notna()]

        if sra_df.empty:
            return sra_df

        if len(sra_df.columns) == 1:
            # groupby().agg({}) raises, so a frame of nothing but the key dedupes
            sra_df = sra_df.drop_duplicates().reset_index(drop=True)
        else:
            sra_df = collapse_dataframe(sra_df, group_col="sample_name")

        # Beside sample_name rather than at the end: this is what identifies the
        # row's libraries, now the key names their group instead. Positional, so
        # it cannot join the append-ordered lists below.
        if library_id:
            sra_df.insert(1, "library_id", sra_df["sample_name"].map(library_id))

        # Required: every builder above seeds a cell for every group, so these
        # are always assigned - a blank one would be rejected by BioSample
        required = (
            ("*isolate", isolate),
            ("*age", age),
            ("*sex", sex),
            ("*tissue", tissue),
            ("*biomaterial_provider", provider),
            ("*collection_date", collection_date),
            ("*geo_loc_name", geo_loc_name),
        )
        # Optional: the column is left out entirely when no group has a value,
        # but once one does, every group gets a cell
        optional = (
            ("suspension_type", suspension_type),
            ("experimental_condition_summary", condition_summary),
            ("intended_cell_type", intended_cell_type),
            ("suspension_enriched_cell_types", enriched_cell_types),
            ("suspension_depleted_cell_types", depleted_cell_types),
            ("suspension_selection_kits", selection_kits),
            ("suspension_enrichment_factors", enrichment_factors),
            ("genetic_perturbation_strategy", genetic_strategy),
            ("preservation_method", preservation_method),
            ("ethnicity", ethnicity),
            ("age_lower_bound", age_lower),
            ("age_upper_bound", age_upper),
        )
        for column, cells in required:
            sra_df[column] = sra_df["sample_name"].map(cells)
        for column, cells in optional:
            if cells:
                sra_df[column] = sra_df["sample_name"].map(cells)

        return sra_df

    def create_sequence_file_set_dataframe(self, complete_data) -> pd.DataFrame:
        """
        One row per SequenceFileSet: library link and filenames from linked files.

        Filenames come from index1, index2, read1, read2, read3, trimmed_cram,
        untrimmed_cram in that order. Each SequenceFile prefers sra_file_names,
        else the s3_uri basename.
        Intermediate only; not written as a CSV.
        """
        if not complete_data:
            return pd.DataFrame()

        resolved = complete_data.get("resolved_objects") or {}
        seq_files_by_id = resolved.get("SequenceFile") or {}
        rows_by_set_id = {}

        for lib_data in (complete_data.get("libraries") or {}).values():
            for raw_file in lib_data.get("raw_matrix_files") or []:
                local_files = {
                    sf.get("@id"): sf
                    for sf in (raw_file.get("sequence_files") or [])
                    if isinstance(sf, dict) and sf.get("@id")
                }
                for file_set in raw_file.get("sequence_file_sets") or []:
                    if not isinstance(file_set, dict):
                        continue
                    set_id = file_set.get("@id")
                    if not set_id or set_id in rows_by_set_id:
                        continue
                    library_ref = file_set.get("library", "")
                    if isinstance(library_ref, dict):
                        library_atid = library_ref.get("@id") or None
                    else:
                        library_atid = library_ref or None
                    names = []
                    total_size = 0
                    for slot in SEQUENCE_FILE_SET_FILE_SLOTS:
                        ref = file_set.get(slot)
                        if is_empty(ref):
                            continue
                        file_id = ref.get("@id") if isinstance(ref, dict) else ref
                        if not file_id:
                            continue
                        seq_file = local_files.get(file_id) or seq_files_by_id.get(file_id)
                        if not seq_file:
                            continue
                        names.extend(self._sequence_file_filenames(seq_file))
                        total_size += self._coerce_file_size(seq_file.get("file_size"))
                    names = list(dict.fromkeys(names))
                    row = {
                        "sequence_file_set_@id": set_id,
                        "library_@id": library_atid,
                        "file_size": total_size,
                    }
                    for i, name in enumerate(names):
                        col = "filename" if i == 0 else f"filename{i + 1}"
                        row[col] = name
                    rows_by_set_id[set_id] = row

        if not rows_by_set_id:
            return pd.DataFrame()
        return pd.DataFrame(list(rows_by_set_id.values()))

    @classmethod
    def _sequence_file_filenames(cls, seq_file: dict) -> list[str]:
        """Filenames for one SequenceFile: prefer sra_file_names, else s3_uri basename."""
        sra_names = seq_file.get("sra_file_names")
        if not is_empty(sra_names):
            return cls._sra_filename_parts(sra_names, basename=False)
        s3_uri = seq_file.get("s3_uri")
        if not is_empty(s3_uri):
            return cls._sra_filename_parts(s3_uri, basename=True)
        return []

    @staticmethod
    def _coerce_file_size(val) -> int:
        """Parse a single file_size cell to a non-negative int; missing/invalid → 0."""
        if is_empty(val):
            return 0
        try:
            return max(0, int(float(val)))
        except (TypeError, ValueError):
            return 0

    @classmethod
    def _sra_file_size_from_sequence_file_sets(cls, sra_rows, sequence_file_set_df) -> pd.Series:
        """Sum SequenceFileSet file_size totals per library (same source as CRISPR)."""
        if (
            sra_rows.empty
            or "_library_atid" not in sra_rows.columns
            or sequence_file_set_df is None
            or sequence_file_set_df.empty
            or "library_@id" not in sequence_file_set_df.columns
            or "file_size" not in sequence_file_set_df.columns
        ):
            return pd.Series(0, index=sra_rows.index, dtype=int)

        sizes = sequence_file_set_df.copy()
        sizes["file_size"] = sizes["file_size"].map(cls._coerce_file_size)
        by_library = sizes.groupby("library_@id", sort=False)["file_size"].sum()
        return (
            sra_rows["_library_atid"]
            .map(by_library)
            .fillna(0)
            .map(cls._coerce_file_size)
            .astype(int)
        )

    @staticmethod
    def _drop_empty_sra_filename_columns(df: pd.DataFrame) -> pd.DataFrame:
        """Drop filename / filenameN columns that are empty for every row."""
        if df is None or df.empty:
            return df
        drop_cols = []
        for col in df.columns:
            if col != "filename" and not re.fullmatch(r"filename\d+", str(col)):
                continue
            if all(is_empty(v) for v in df[col].tolist()):
                drop_cols.append(col)
        return df.drop(columns=drop_cols) if drop_cols else df

    def _write_sra_chunked_outputs(self, sra_dir, sra_files_df, sra_biosample_df):
        """Write sra_file_setN.csv / sra_biosample_setN.csv under ≤5 TB sample_name packs."""
        chunks = self._chunk_sra_by_sample_name(sra_files_df)
        for set_n, file_chunk in enumerate(chunks, start=1):
            write_files = file_chunk.drop(columns=["file_size"], errors="ignore")
            write_files = self._drop_empty_sra_filename_columns(write_files)
            sample_names = set()
            if "sample_name" in write_files.columns:
                sample_names = {
                    name for name in write_files["sample_name"].tolist() if not is_empty(name)
                }
            if (
                sra_biosample_df is not None
                and not sra_biosample_df.empty
                and "sample_name" in sra_biosample_df.columns
            ):
                bio_chunk = sra_biosample_df[
                    sra_biosample_df["sample_name"].isin(sample_names)
                ].copy()
            else:
                bio_chunk = pd.DataFrame()

            file_output = sra_dir / f"sra_file_set{set_n}.csv"
            bio_output = sra_dir / f"sra_biosample_set{set_n}.csv"
            print(f"Saving SRA file DataFrame to {file_output}...")
            write_files.to_csv(file_output, index=False)
            print(f"✅ SRA file CSV file created: {file_output}")
            print(f"   Rows: {len(write_files)}")
            print(f"   Columns: {len(write_files.columns)}")
            print(f"Saving SRA/BioSample DataFrame to {bio_output}...")
            bio_chunk.to_csv(bio_output, index=False)
            print(f"✅ SRA/BioSample CSV file created: {bio_output}")
            print(f"   Rows: {len(bio_chunk)}")
            print(f"   Columns: {len(bio_chunk.columns)}")

    @classmethod
    def _chunk_sra_by_sample_name(
        cls, sra_files_df, max_bytes: int = SRA_CHUNK_MAX_BYTES
    ) -> list[pd.DataFrame]:
        """Pack SRA_FILE rows by sample_name into chunks totaling ≤ max_bytes."""
        if sra_files_df is None or sra_files_df.empty:
            return [pd.DataFrame() if sra_files_df is None else sra_files_df.copy()]

        df = sra_files_df.copy()
        if "sample_name" not in df.columns:
            print("Warning: SRA_FILE has no sample_name; writing as a single chunk")
            return [df]

        if "file_size" not in df.columns:
            df["file_size"] = 0
        df["file_size"] = df["file_size"].map(cls._coerce_file_size)

        order = list(dict.fromkeys(df["sample_name"].tolist()))
        groups = []
        for name in order:
            part = df.loc[df["sample_name"] == name]
            groups.append((name, part, int(part["file_size"].sum())))

        chunks: list[pd.DataFrame] = []
        current_parts: list[pd.DataFrame] = []
        current_size = 0
        for name, part, size in groups:
            if current_parts and current_size + size > max_bytes:
                chunks.append(pd.concat(current_parts, ignore_index=True))
                current_parts = []
                current_size = 0
            if not current_parts and size > max_bytes:
                print(
                    f"Warning: sample_name {name!r} totals {size} bytes "
                    f"(> {max_bytes}); writing alone in its own SRA set"
                )
                chunks.append(part.reset_index(drop=True))
                continue
            current_parts.append(part)
            current_size += size
        if current_parts:
            chunks.append(pd.concat(current_parts, ignore_index=True))
        return chunks or [df.iloc[0:0].copy()]

    def create_sra_files_dataframe(self, main_df, sequence_file_set_df=None) -> pd.DataFrame:
        """
        Build the SRA file dataframe from main_df: one row per CRO group × library.

        sample_name is the same CRO_group_identifier as SRA_BIOSAMPLE. Unlike that
        sheet, every library in the group gets its own row, including non-GEX.
        Duplicate MAIN rows (same library on more than one raw matrix file)
        collapse to one. CRISPR libraries expand to one row per SequenceFileSet
        when sequence_file_set_df is provided; non-CRISPR keep MAIN filenames.
        """
        group_cols = [
            col
            for col in (
                "droplet_based_libraries_CRO_group_identifier",
                "plate_based_libraries_CRO_group_identifier",
            )
            if col in main_df.columns
        ]
        if not group_cols:
            print(
                "Warning: MAIN has no library CRO group identifier column; SRA_FILE will be empty"
            )
            return pd.DataFrame()

        if main_df.empty:
            return pd.DataFrame()

        title_source_cols = []
        seen = set()
        for col in (
            "raw_file_samples",
            "genetic_modifications_strategy",
            EXPERIMENTAL_CONDITION_SUMMARY_COL,
            CRO_ORDER_COL,
        ):
            if col in main_df.columns and col not in seen:
                seen.add(col)
                title_source_cols.append(col)

        if "sequence_files_sra_file_names" in main_df.columns:
            filename_source = "sequence_files_sra_file_names"
            filename_basename = False
        elif "sequence_files_s3_uri" in main_df.columns:
            filename_source = "sequence_files_s3_uri"
            filename_basename = True
        else:
            filename_source = None
            filename_basename = False

        columns_to_keep = [k for k in PROP_MAP_SRA_FILE if k in main_df.columns]
        for col in title_source_cols:
            if col not in columns_to_keep:
                columns_to_keep.append(col)
        if filename_source and filename_source not in columns_to_keep:
            columns_to_keep.append(filename_source)
        sra_df = main_df[columns_to_keep].copy()
        sra_df.rename(columns=PROP_MAP_SRA_FILE, inplace=True)
        if "raw_file_samples" in sra_df.columns:
            sra_df.rename(columns={"raw_file_samples": "samples"}, inplace=True)
        sra_df = collapse_duplicate_columns(sra_df)

        aliases = self._coalesce_columns(
            main_df, [f"{kind}_based_libraries_aliases" for kind in ("droplet", "plate")]
        )
        library_atid = self._coalesce_columns(
            main_df, [f"{kind}_based_libraries_@id" for kind in ("droplet", "plate")]
        )
        if aliases is None:
            print("Warning: MAIN has no library aliases column; SRA_FILE omits library_ID")
        else:
            sra_df["library_ID"] = aliases.map(
                lambda alias: self._get_clean_alias({"aliases": to_items(alias)})
            )

        feature_types = self._coalesce_columns(
            main_df,
            [f"{kind}_based_libraries_feature_types" for kind in ("droplet", "plate")],
        )
        if feature_types is None:
            sra_df["library_strategy"] = None
            sra_df["library_source"] = None
            sra_df["library_selection"] = None
            sra_df["_feature_types"] = None
        else:
            sra_df["library_strategy"] = feature_types.map(self._map_sra_library_strategy)
            sra_df["library_source"] = feature_types.map(self._map_sra_library_source)
            sra_df["library_selection"] = feature_types.map(self._map_sra_library_selection)
            sra_df["_feature_types"] = feature_types

        run_cardinality = main_df.get("sequence_file_sets_run_cardinality")
        if run_cardinality is None:
            sra_df["library_layout"] = None
        else:
            sra_df["library_layout"] = run_cardinality.map(self._map_sra_library_layout)

        sequencing_platform = main_df.get("sequence_file_sets_sequencing_platform")
        if sequencing_platform is None:
            sra_df["platform"] = None
        else:
            sra_df["platform"] = sequencing_platform.map(self._map_sra_platform)
        if "instrument_model" not in sra_df.columns:
            sra_df["instrument_model"] = None
        if "filetype" not in sra_df.columns:
            sra_df["filetype"] = None

        barcode_map = None
        for samples_col in ("droplet_based_libraries_samples", "plate_based_libraries_samples"):
            if samples_col not in main_df.columns:
                continue
            mapped = main_df[samples_col].map(self._sample_probe_barcode_map)
            barcode_map = mapped if barcode_map is None else barcode_map.fillna(mapped)
        if barcode_map is not None and barcode_map.map(lambda v: not is_empty(v)).any():
            barcode_column = (
                "sample_name: hash_index"
                if self._frame_has_plate_library(main_df)
                else "sample_name: sample_probe_barcode"
            )
            sra_df[barcode_column] = barcode_map

        if library_atid is not None:
            dedup_key = library_atid
            if "library_ID" in sra_df.columns:
                dedup_key = dedup_key.fillna(sra_df["library_ID"])
            sra_df["_dedup"] = dedup_key
        elif "library_ID" in sra_df.columns:
            sra_df["_dedup"] = sra_df["library_ID"]

        unnamed = int(sra_df["sample_name"].isna().sum())
        if unnamed:
            print(
                f"Warning: dropping {unnamed} of {len(sra_df)} MAIN row(s) with no library "
                "CRO group identifier from SRA_FILE"
            )
            sra_df = sra_df[sra_df["sample_name"].notna()]

        if sra_df.empty:
            return sra_df.drop(columns=["_dedup"], errors="ignore").reset_index(drop=True)

        if "_dedup" in sra_df.columns:
            sra_df["_dedup"] = sra_df["_dedup"].fillna(sra_df["sample_name"])
            sra_df = collapse_dataframe(sra_df, group_col="_dedup", sort=False)
            sra_df = sra_df.rename(columns={"_dedup": "_library_atid"})
        else:
            sra_df = collapse_dataframe(sra_df, group_col="sample_name", sort=False)

        col = "genetic_modifications_strategy"
        if col in sra_df.columns:
            sra_df[col] = sra_df[col].map(self._map_geo_genetic_strategy)

        sra_df = self._add_geo_title(
            sra_df,
            title_col="title",
            name_col="sample_name",
            strategy_col="library_strategy",
        )

        sra_df["design_description"] = ""
        if "_feature_types" in sra_df.columns:
            sra_df["_is_crispr"] = sra_df["_feature_types"].map(self._sra_feature_has_crispr)
            if bool(sra_df["_is_crispr"].any()):
                sra_df["options"] = sra_df["_is_crispr"].map(
                    lambda is_crispr: "--ignoreNames" if is_crispr else ""
                )
        else:
            sra_df["_is_crispr"] = False

        drop_cols = [
            c
            for c in (
                "samples",
                "_feature_types",
            )
            if c in sra_df.columns
        ]
        if drop_cols:
            sra_df = sra_df.drop(columns=drop_cols)

        ordered = [
            c
            for c in (
                "sample_name",
                "library_ID",
                "library_strategy",
                "library_source",
                "library_selection",
                "library_layout",
                "platform",
                "instrument_model",
                "filetype",
                "title",
                "design_description",
                "options",
                "sample_name: sample_probe_barcode",
                "sample_name: hash_index",
                "donor_ids",
                "donor_dev_stage",
                "tissue",
                "genetic_modifications_strategy",
                "_library_atid",
                "_is_crispr",
            )
            if c in sra_df.columns
        ]
        if filename_source and filename_source in sra_df.columns:
            ordered.append(filename_source)
        sra_df = sra_df[ordered]
        sra_df = self._apply_sra_filenames(
            sra_df,
            sequence_file_set_df=sequence_file_set_df,
            filename_source=filename_source,
            filename_basename=filename_basename,
        )
        sort_cols = [c for c in ("sample_name", "sequence_file_set_@id") if c in sra_df.columns]
        if sort_cols:
            sra_df = sra_df.sort_values(sort_cols, kind="mergesort")
        helper_cols = [
            c
            for c in (
                "_library_atid",
                "_is_crispr",
                "library_@id",
                "sequence_file_set_@id",
            )
            if c in sra_df.columns
        ]
        if helper_cols:
            sra_df = sra_df.drop(columns=helper_cols)
        if "file_size" in sra_df.columns:
            sra_df["file_size"] = sra_df["file_size"].map(self._coerce_file_size)
        else:
            sra_df["file_size"] = 0
        return sra_df.reset_index(drop=True)

    @staticmethod
    def _map_sra_feature_type(ft, mapping):
        """Map the first listed feature type through mapping, or None if unmapped."""
        for item in to_items(ft):
            mapped = mapping.get(item)
            if mapped:
                return mapped
        return None

    @classmethod
    def _map_sra_library_strategy(cls, ft):
        """Map library feature_types to an SRA library_strategy, or None if unmapped."""
        return cls._map_sra_feature_type(ft, SRA_FILE_LIBRARY_STRATEGY_MAP)

    @classmethod
    def _map_sra_library_source(cls, ft):
        """Map library feature_types to an SRA library_source, or None if unmapped."""
        return cls._map_sra_feature_type(ft, SRA_FILE_LIBRARY_SOURCE_MAP)

    @classmethod
    def _map_sra_library_selection(cls, ft):
        """Map library feature_types to an SRA library_selection, or None if unmapped."""
        return cls._map_sra_feature_type(ft, SRA_FILE_LIBRARY_SELECTION_MAP)

    @staticmethod
    def _sra_feature_has_crispr(ft) -> bool:
        """True when feature_types include CRISPR Guide Capture."""
        return any("CRISPR Guide Capture" in str(item) for item in to_items(ft))

    @classmethod
    def _map_sra_library_layout(cls, val):
        """Map sequence file set run_cardinality to an SRA library_layout, or None."""
        if is_empty(val):
            return None
        mapped = []
        for item in to_items(val):
            for part in str(item).split(";"):
                part = part.strip()
                if not part:
                    continue
                mapped.append(SRA_FILE_LIBRARY_LAYOUT_MAP.get(part))
        if not mapped:
            return None
        unique = {item for item in mapped if item is not None}
        if None not in mapped and unique in ({"paired"}, {"single"}):
            return next(iter(unique))
        print(
            f"Warning: SRA_FILE cannot map run_cardinality {val!r} to a single "
            "library_layout (paired or single)"
        )
        return None

    @staticmethod
    def _map_sra_platform(val):
        """Map sequencing_platform to an SRA platform by vendor prefix, or None."""
        for item in to_items(val):
            if is_empty(item):
                continue
            text = str(item)
            for prefix, platform in SRA_FILE_PLATFORM_PREFIXES:
                if text.startswith(prefix):
                    return platform
        return None

    @classmethod
    def _sra_filename_parts(cls, val, *, basename: bool) -> list[str]:
        """Split a joined or collapsed filename cell into unique names."""
        names = []
        for item in to_items(val):
            for part in str(item).split("; "):
                part = part.strip()
                if not part:
                    continue
                names.append(PurePosixPath(part).name if basename else part)
        return list(dict.fromkeys(names))

    @classmethod
    def _add_sra_filenames(cls, sra_df, source_col, *, basename: bool) -> pd.DataFrame:
        """Replace source_col with filename, filename2, filename3, … on the right."""
        if source_col not in sra_df.columns:
            return sra_df
        names = sra_df[source_col].map(lambda val: cls._sra_filename_parts(val, basename=basename))
        sra_df = sra_df.drop(columns=[source_col])
        max_n = int(names.map(len).max()) if len(names) else 0
        if max_n == 0:
            return sra_df
        for i in range(max_n):
            col = "filename" if i == 0 else f"filename{i + 1}"
            sra_df[col] = names.map(lambda xs, i=i: xs[i] if i < len(xs) else pd.NA)
        return sra_df

    @classmethod
    def _apply_sra_filenames(
        cls,
        sra_df,
        *,
        sequence_file_set_df,
        filename_source,
        filename_basename: bool,
    ) -> pd.DataFrame:
        """
        Attach filenames and file_size for SRA rows.

        Non-CRISPR: MAIN filenames; file_size is the sum of matching SequenceFileSet
        totals (same SequenceFile object sums as CRISPR), one row per library.
        CRISPR: expand to one row per SequenceFileSet when matches exist; keep one
        row with no filenames when none do.
        """
        if sra_df.empty:
            return sra_df

        is_crispr = (
            sra_df["_is_crispr"].fillna(False).astype(bool)
            if "_is_crispr" in sra_df.columns
            else pd.Series(False, index=sra_df.index)
        )
        non_crispr = sra_df.loc[~is_crispr].copy()
        crispr = sra_df.loc[is_crispr].copy()

        if not non_crispr.empty:
            if filename_source and filename_source in non_crispr.columns:
                non_crispr = cls._add_sra_filenames(
                    non_crispr, filename_source, basename=filename_basename
                )
            non_crispr["file_size"] = cls._sra_file_size_from_sequence_file_sets(
                non_crispr, sequence_file_set_df
            )

        if not crispr.empty:
            if filename_source and filename_source in crispr.columns:
                crispr = crispr.drop(columns=[filename_source])
            if (
                sequence_file_set_df is not None
                and not sequence_file_set_df.empty
                and "_library_atid" in crispr.columns
                and "library_@id" in sequence_file_set_df.columns
            ):
                crispr = crispr.merge(
                    sequence_file_set_df,
                    left_on="_library_atid",
                    right_on="library_@id",
                    how="left",
                )
                crispr = cls._suffix_crispr_library_ids(crispr)
            if "file_size" not in crispr.columns:
                crispr["file_size"] = 0
            else:
                crispr["file_size"] = crispr["file_size"].map(cls._coerce_file_size)

        parts = [part for part in (non_crispr, crispr) if not part.empty]
        if not parts:
            return sra_df
        return pd.concat(parts, ignore_index=True, sort=False)

    @staticmethod
    def _suffix_crispr_library_ids(crispr: pd.DataFrame) -> pd.DataFrame:
        """Append .1, .2, … to library_ID when a CRISPR library has multiple SequenceFileSets.

        A library with only one matched SequenceFileSet keeps the unsuffixed library_ID.
        """
        if crispr.empty or "library_ID" not in crispr.columns:
            return crispr
        if "sequence_file_set_@id" not in crispr.columns:
            return crispr

        crispr = crispr.copy()
        group_col = "_library_atid" if "_library_atid" in crispr.columns else None
        sort_cols = [c for c in (group_col, "sequence_file_set_@id") if c]
        crispr = crispr.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

        matched = crispr["sequence_file_set_@id"].map(lambda v: not is_empty(v))
        if not matched.any():
            return crispr
        matched_rows = crispr.loc[matched]
        if group_col:
            group_sizes = matched_rows.groupby(group_col, sort=False)[group_col].transform("size")
            suffixes = matched_rows.groupby(group_col, sort=False).cumcount().add(1)
            multi = group_sizes.gt(1)
        else:
            suffixes = pd.Series(range(1, len(matched_rows) + 1), index=matched_rows.index)
            multi = pd.Series(len(matched_rows) > 1, index=matched_rows.index)
        suffix_idx = multi[multi].index
        if suffix_idx.empty:
            return crispr
        base = matched_rows.loc[suffix_idx, "library_ID"].map(
            lambda v: "" if is_empty(v) else str(v)
        )
        crispr.loc[suffix_idx, "library_ID"] = base + "." + suffixes.loc[suffix_idx].astype(str)
        return crispr

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
    def _format_pooled_values(values, sep: str = ", ") -> str | None:
        """'pooled: a, b' across several values, or the bare value for a lone one."""
        entries = [value for value in values if value]
        if not entries:
            return None
        if len(entries) == 1:
            return entries[0]
        return "pooled: " + sep.join(entries)

    def _library_ids_by_group(self, main_df, gex_mask):
        """
        Per-group 'library_id' cell: every library in the group, GEX ones first.

        The one builder taking the unfiltered frame, because the non-GEX
        libraries it reports are exactly the rows the GEX filter removes. Both
        halves of a group carry the same CRO identifier, so grouping still works
        on the wider frame. gex_mask is passed in rather than recomputed to keep
        one answer to what counts as GEX.

        Aliases arrive with the lab prefix in MAIN, so they are cleaned. Sorted
        within each half, but GEX before the rest: plain sorting would lead with
        the CRISPR library, whose alias suffix happens to come first.
        """
        groups = self._coalesce_columns(
            main_df,
            [f"{kind}_based_libraries_CRO_group_identifier" for kind in ("droplet", "plate")],
        )
        aliases = self._coalesce_columns(
            main_df, [f"{kind}_based_libraries_aliases" for kind in ("droplet", "plate")]
        )
        if aliases is None:
            print("Warning: MAIN has no library aliases column; SRA_BIOSAMPLE omits library_id")
            return {}

        # {group: (gex aliases, everything else)}
        by_group: dict[str, tuple[set[str], set[str]]] = {}
        for group, alias, is_gex in zip(groups, aliases, gex_mask, strict=True):
            if not isinstance(group, str):
                continue
            cleaned = self._get_clean_alias({"aliases": to_items(alias)})
            if not cleaned:
                continue
            gex, other = by_group.setdefault(group, (set(), set()))
            (gex if is_gex else other).add(cleaned)

        return {
            group: ", ".join(sorted(gex) + sorted(other))
            for group, (gex, other) in by_group.items()
        }

    def _biomaterial_provider_by_library(self, main_df, library_key):
        """
        Build the 'biomaterial_provider' cell per library, keyed by CRO group.

        Falls back per row, not per column: a sample with no usable 'sources' uses
        its own 'lab' even where a sibling has one. Distinct titles join with '; '.
        Missing columns are still warned about, but do not drop the column: this
        one is required, so it fills with 'not provided' instead.
        """
        prefixes = SAMPLE_URL_PREFIXES
        sources = self._coalesce_columns(main_df, [f"{prefix}_sources" for prefix in prefixes])
        labs = self._coalesce_columns(main_df, [f"{prefix}_lab" for prefix in prefixes])
        if sources is None and labs is None:
            print(
                "Warning: MAIN has no sample sources or lab column; "
                "SRA_BIOSAMPLE reports biomaterial_provider as 'not provided'"
            )

        empty = pd.Series(None, index=main_df.index, dtype=object)
        sources = empty if sources is None else sources
        labs = empty if labs is None else labs

        titles_by_library: dict[str, set[str]] = {}
        for library, source, lab in zip(library_key, sources, labs, strict=True):
            if not isinstance(library, str):
                continue
            # Seeded for every group, so a group with nothing still gets a cell:
            # this column is required, and a blank one would be rejected
            found = titles_by_library.setdefault(library, set())
            titles = self._provider_titles(source) or self._provider_titles(lab)
            found.update(titles or ["not provided"])

        return {library: "; ".join(sorted(titles)) for library, titles in titles_by_library.items()}

    @staticmethod
    def _tissue_row_mask(df) -> pd.Series:
        """True where the row is a Tissue (ethnicity is tissue-only)."""
        col = "tissues_@id"
        if col not in df.columns:
            return pd.Series(False, index=df.index)
        return df[col].map(lambda v: not is_empty(v))

    def _add_samples_donor_ethnicity(self, sample_df) -> pd.DataFrame:
        """SRA-style pooled ethnicity for tissue rows only; drop source columns."""
        source = "human_donors_ethnicity_term_name"
        tissue = self._tissue_row_mask(sample_df)
        drop = [c for c in (source, "tissues_@id") if c in sample_df.columns]
        if source not in sample_df.columns:
            return sample_df.drop(columns=drop) if drop else sample_df

        def _cell(val):
            texts = self._field_texts((val,), None, True, None)
            return self._format_pooled_values(sorted(texts)) if texts else None

        values = sample_df[source].map(_cell).where(tissue)
        sample_df = sample_df.drop(columns=drop)
        if not values.notna().any():
            return sample_df
        values = values.mask(tissue & values.isna(), "not provided")
        sample_df["donor_ethnicity"] = values
        return sample_df

    def _sample_field_by_library(
        self,
        main_df,
        library_key,
        *fields,
        optional=False,
        prefixes=SAMPLE_URL_PREFIXES,
        gap="not provided",
        translate=None,
        split_joined=False,
        transform=None,
        formatter=None,
        include=None,
    ):
        """
        Per-library cell for a field on the sample or on something it links to.

        Every cell on this sheet is an unordered set of the distinct values its
        library's samples carry. A sample with no value contributes 'gap', unless
        gap is None. An array field contributes each of its items separately.

        'prefixes' narrows which MAIN column prefixes are read - the sample types
        for a field on the sample itself, or a linked object's prefix.
        'split_joined' splits items on '; ', which is how create_dataframe
        collapses a field when a sample references several of an object; only set
        it where that is the source, so a value legitimately containing '; ' is
        never torn apart.
        'transform' combines one row's fields into a single value, for a value
        that needs its units. 'translate' rewrites each value through a lookup.
        'formatter' renders the finished list, defaulting to a '; ' join.
        'optional' decides only whether the column exists: nothing anywhere means
        nothing is returned and the column drops.
        'include' is an optional row mask: False skips that row so it neither
        contributes a value nor the gap.
        """
        empty = pd.Series(None, index=main_df.index, dtype=object)
        sources = []
        for field in fields:
            found_column = self._coalesce_columns(
                main_df, [f"{prefix}_{field}" for prefix in prefixes]
            )
            sources.append(empty if found_column is None else found_column)

        include_flags = (
            [True] * len(main_df) if include is None else [bool(flag) for flag in include]
        )
        found: dict[str, set[str]] = {}
        missing: set[str] = set()
        for library, keep, *row in zip(library_key, include_flags, *sources, strict=True):
            if not keep or not isinstance(library, str):
                continue
            found.setdefault(library, set())
            texts = self._field_texts(row, transform, split_joined, translate)
            if texts:
                found[library].update(texts)
            else:
                missing.add(library)

        if optional and not any(found.values()):
            return {}

        render = formatter or "; ".join
        cells = {}
        for library, present in found.items():
            values = sorted(present)
            if gap is not None and (library in missing or not values):
                values.append(gap)
            cells[library] = render(values)
        return cells

    @staticmethod
    def _field_texts(row, transform, split_joined, translate) -> list[str]:
        """One row's fields as zero or more cell values."""
        if transform:
            # to_items so a transform may return one value or several: a cell
            # describing two subjects yields an age each, not just the first
            return [str(text).strip() for text in to_items(transform(*row)) if not is_empty(text)]

        # numeric_text before the split, or str() stamps a float's '.0' on first
        texts = [numeric_text(item) for item in to_items(row[0])]
        if split_joined:
            texts = [part for text in texts for part in text.split("; ")]
        if translate:
            texts = [translate.get(text, text) for text in texts]
        return [text for text in texts if text]

    # Sample types with no tissue of origin to report
    TISSUELESS_SAMPLE_TYPES = ("cell_lines", "primary_cell_cultures")

    def _tissue_by_library(self, main_df, library_key):
        """
        Per-library 'tissue' cell, keyed by CRO group.

        A cell line or primary cell culture has no tissue, so it reports
        'not applicable'. A tissue or organoid reports its sample_terms, which is
        an array - each term becomes its own entry. A sample with no term at all
        contributes 'not provided', a different thing from 'not applicable': the
        schema should not allow it, but this column is required and cannot be
        left blank. A MAIN with no sample_terms column at all fills the same way,
        rather than dropping a required column.

        Both markers are INSDC missing-value terms, which is what BioSample
        validates a required attribute against.
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
            found = by_library.setdefault(library, set())
            if tissueless_row:
                found.add("not applicable")
                continue
            found.update([str(item).strip() for item in to_items(term)] or ["not provided"])

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

    def create_guide_metadata_dataframe(self, file_info):
        """
        Build a guide-metadata DataFrame from one guide RNA TabularFile.

        file_info is the TabularFile object (or embed) for a single BIOHUB set.
        None returns None.
        """
        if not file_info:
            return None

        guide_df = DB2lattice.read_tabular_file(file_info, self.connection)
        present = [col for col in GUIDE_METADATA_COLUMN_MAP if col in guide_df.columns]
        if not present:
            print("Warning: guide RNA file has none of the expected GUIDE_METADATA columns")
            return None
        return guide_df[present].rename(columns=GUIDE_METADATA_COLUMN_MAP).copy()

    def _write_biohub_outputs(
        self,
        biohub_dir,
        geo_dir,
        main_df,
        biohub_samples,
        complete_data,
        *,
        biohub_version: str,
    ):
        """Write BIOHUB samples/guide templates, split by GeneticModification set.

        Guide frames are also copied to GEO as ``{alias}_guide_metadata.csv``.
        """
        source = self.biohub_source_dataframe(main_df, biohub_samples)
        sets = self._biohub_guide_sets(complete_data)

        if not sets:
            biohub_df = self.create_biohub_dataframe(source)
            biohub_output = biohub_dir / f"samples_template_set1_{biohub_version}.csv"
            print(f"Saving biohub DataFrame to {biohub_output}...")
            biohub_df.to_csv(biohub_output, index=False)
            print(f"✅ Biohub CSV file created: {biohub_output}")
            return

        assigned_index = set()
        set_subsets = []
        for biohub_set in sets:
            library_ids = biohub_set["library_ids"]
            if source.empty:
                subset = source
            else:
                mask = source.apply(
                    lambda row, libs=library_ids: bool(self._row_library_ids(row) & libs),
                    axis=1,
                )
                subset = source.loc[mask]
                assigned_index.update(subset.index.tolist())
            set_subsets.append((biohub_set, subset))

        if not source.empty:
            unassigned = len(source.index.difference(assigned_index))
            if unassigned:
                message = (
                    f"{unassigned} BIOHUB source row(s) have no matching GeneticModification set"
                )
                print(f"Warning: {message}")
                raise ValueError(message)

        for biohub_set, subset in set_subsets:
            set_n = biohub_set["set_number"]
            biohub_df = self.create_biohub_dataframe(subset)
            samples_output = biohub_dir / f"samples_template_set{set_n}_{biohub_version}.csv"
            print(f"Saving biohub DataFrame to {samples_output}...")
            biohub_df.to_csv(samples_output, index=False)
            print(f"✅ Biohub CSV file created: {samples_output}")

            file_info = biohub_set["file_info"]
            if not file_info:
                continue
            guide_df = self.create_guide_metadata_dataframe(file_info)
            if guide_df is None or guide_df.empty:
                continue
            guide_output = biohub_dir / f"guide_template_set{set_n}_{biohub_version}.csv"
            print(f"Saving guide metadata DataFrame to {guide_output}...")
            guide_df.to_csv(guide_output, index=False)
            print(f"✅ Guide metadata CSV file created: {guide_output}")
            print(f"   Rows: {len(guide_df)}")
            print(f"   Columns: {len(guide_df.columns)}")

            alias = biohub_set["alias"]
            geo_guide_output = geo_dir / f"{alias}_guide_metadata.csv"
            print(f"Saving GEO guide metadata DataFrame to {geo_guide_output}...")
            guide_df.to_csv(geo_guide_output, index=False)
            print(f"✅ GEO guide metadata CSV file created: {geo_guide_output}")

    def _biohub_guide_sets(self, complete_data):
        """
        One BIOHUB set per GeneticModification (sorted by @id).

        Each set expects exactly one guide TabularFile. Divergence (0 or >1 files
        on a GM, or a TabularFile claimed by two GMs) is warned and that set's
        guide file is skipped; samples are still written for the set.
        """
        resolved = complete_data.get("resolved_objects") or {}
        genetic_modifications = resolved.get("GeneticModification") or {}
        if not genetic_modifications:
            return []

        gathered = resolved.get("TabularFile") or {}
        gm_to_libraries = {gm_id: set() for gm_id in genetic_modifications}
        for lib_data in (complete_data.get("libraries") or {}).values():
            library = lib_data.get("library") or {}
            lib_atid = library.get("@id")
            if not lib_atid:
                continue
            for gm in lib_data.get("genetic_modifications") or []:
                gm_id = gm.get("@id") if isinstance(gm, dict) else None
                if gm_id in gm_to_libraries:
                    gm_to_libraries[gm_id].add(lib_atid)

        library_to_gms = {}
        for gm_id, libs in gm_to_libraries.items():
            for lib_atid in libs:
                library_to_gms.setdefault(lib_atid, set()).add(gm_id)
        for lib_atid, gm_ids in sorted(library_to_gms.items()):
            if len(gm_ids) > 1:
                raise ValueError(
                    f"Library {lib_atid} is linked to multiple GeneticModifications: "
                    f"{sorted(gm_ids)}; expected one GM per library"
                )

        claimed_files = {}
        sets = []
        for set_number, gm_id in enumerate(sorted(genetic_modifications), start=1):
            gm = genetic_modifications[gm_id]
            file_entries = []
            for file_info in normalize_guide_rna_file_refs(gm.get("guide_rna_files")):
                key = file_info.get("@id") or file_info.get("s3_uri")
                if not key:
                    continue
                resolved_file = gathered.get(file_info.get("@id"), file_info)
                file_entries.append((key, resolved_file))

            file_info = None
            if not file_entries:
                print(
                    f"Warning: GeneticModification {gm_id} has no guide_rna_files; "
                    f"skipping guide_template for set{set_number}"
                )
            elif len(file_entries) > 1:
                labels = [key for key, _ in file_entries]
                print(
                    f"Warning: GeneticModification {gm_id} has {len(file_entries)} "
                    f"guide_rna_files; expected one. Skipping guide_template for "
                    f"set{set_number}. Files: {labels}"
                )
            else:
                key, resolved_file = file_entries[0]
                prior_gm = claimed_files.get(key)
                if prior_gm is not None:
                    print(
                        f"Warning: guide TabularFile {key} is shared by "
                        f"{prior_gm} and {gm_id}; skipping guide_template for "
                        f"set{set_number}"
                    )
                else:
                    claimed_files[key] = gm_id
                    file_info = resolved_file

            sets.append(
                {
                    "set_number": set_number,
                    "gm_id": gm_id,
                    "file_info": file_info,
                    "library_ids": gm_to_libraries[gm_id],
                    "alias": self._genetic_modification_alias_label(gm.get("aliases"), gm_id=gm_id),
                }
            )
        return sets

    @classmethod
    def _genetic_modification_alias_label(cls, aliases, *, gm_id=None) -> str:
        """Require exactly one GM alias; strip lab prefix and sanitize for GEO/filenames."""
        parts = []
        for item in to_items(aliases):
            if is_empty(item):
                continue
            for part in str(item).split(";"):
                part = part.strip()
                if part:
                    parts.append(part)
        # Preserve order while dropping duplicates
        unique = list(dict.fromkeys(parts))
        if len(unique) != 1:
            gm_label = gm_id if gm_id is not None else "GeneticModification"
            raise ValueError(f"{gm_label} must have exactly one alias; found {unique!r}")
        label = unique[0]
        if ":" in label:
            label = label.split(":", 1)[1]
        return re.sub(r"[^A-Za-z0-9._-]", "_", label)

    @classmethod
    def _validate_genetic_modification_aliases(cls, complete_data) -> None:
        """Raise if any GeneticModification lacks exactly one alias."""
        resolved = (complete_data or {}).get("resolved_objects") or {}
        genetic_modifications = resolved.get("GeneticModification") or {}
        for gm_id, gm in genetic_modifications.items():
            aliases = gm.get("aliases") if isinstance(gm, dict) else None
            cls._genetic_modification_alias_label(aliases, gm_id=gm_id)

    @classmethod
    def _row_library_ids(cls, row):
        """Library @ids on a biohub-source row (droplet and/or plate)."""
        ids = set()
        for col in BIOHUB_LIBRARY_ID_COLS:
            if col not in row.index:
                continue
            val = row[col]
            if is_empty(val):
                continue
            for item in to_items(val):
                for part in str(item).split("; "):
                    part = part.strip()
                    if part:
                        ids.add(part)
        return ids

    @staticmethod
    def _frame_has_plate_library(df) -> bool:
        """True when any plate library column on this frame has a value."""
        for column in df.columns:
            if not str(column).startswith("plate_based_libraries_"):
                continue
            if df[column].map(lambda value: not is_empty(value)).any():
                return True
        return False

    def biohub_source_dataframe(self, main_df, sample_df):
        """
        Rows create_biohub_dataframe should read.

        Files with raw_matrix_files_is_multiplexed true contribute the per-sample
        frame (one row per sample, sample_name from the cleaned sample_alias).
        Every other file contributes its main_df rows, where sample_name still
        comes from raw_file_samples. A flagged file with no sample rows stays
        on main_df.
        """
        if (
            sample_df is None
            or sample_df.empty
            or "raw_matrix_file_alias" not in main_df.columns
            or "raw_matrix_file_alias" not in sample_df.columns
            or "raw_matrix_files_is_multiplexed" not in main_df.columns
        ):
            return main_df

        multiplexed = set(
            main_df.loc[
                main_df["raw_matrix_files_is_multiplexed"].eq(True),
                "raw_matrix_file_alias",
            ].dropna()
        )
        if not multiplexed:
            return main_df

        exploded = sample_df.loc[sample_df["raw_matrix_file_alias"].isin(multiplexed)].copy()
        if exploded.empty:
            return main_df

        if "sample_alias" in exploded.columns:
            exploded["sample_name"] = exploded["sample_alias"]
        if "raw_file_samples" in exploded.columns:
            exploded = exploded.drop(columns=["raw_file_samples"])

        # Library and file-set fields are not on the per-sample frame. Collapse
        # to one value per file first so two libraries do not repeat each sample.
        # Author-metadata columns use the same selection as create_biohub_dataframe.
        # A column already on the sample rows is left there, so per-sample values
        # are not replaced by the rolled-up MAIN copy.
        extras = [
            column
            for column in PROP_MAP_BIOHUB
            if column != "raw_file_samples"
            and column in main_df.columns
            and column not in exploded.columns
        ]
        extras.extend(
            column
            for column in main_df.columns
            if re.search("_author_metadata_", column)
            and column not in exploded.columns
            and column not in extras
        )
        # Library @ids: plate recognition after per-sample rows replace main, and
        # BIOHUB set splitting. create_biohub_dataframe does not write these.
        for lib_id_col in BIOHUB_LIBRARY_ID_COLS:
            if (
                lib_id_col in main_df.columns
                and lib_id_col not in exploded.columns
                and lib_id_col not in extras
            ):
                extras.append(lib_id_col)
        if extras:
            file_level = main_df.loc[
                main_df["raw_matrix_file_alias"].isin(exploded["raw_matrix_file_alias"]),
                ["raw_matrix_file_alias", *extras],
            ]
            file_level = file_level.groupby("raw_matrix_file_alias", as_index=False).agg(
                {column: self._combine_sample_values for column in extras}
            )
            exploded = exploded.merge(file_level, on="raw_matrix_file_alias", how="left")

        plain = main_df.loc[
            ~main_df["raw_matrix_file_alias"].isin(set(exploded["raw_matrix_file_alias"]))
        ]
        return pd.concat([plain, exploded], ignore_index=True)

    def create_biohub_dataframe(self, main_df) -> pd.DataFrame:
        """
        Build Biohub samples dataframe from main_df
        """
        plate_experiment = self._frame_has_plate_library(main_df)
        columns_to_keep = [k for k in PROP_MAP_BIOHUB if k in main_df.columns]
        columns_to_keep.extend([k for k in main_df.columns if re.search("_author_metadata_", k)])
        # Multiplexed rows already carry sample_name from sample_alias. Keeping
        # it here lets collapse_duplicate_columns coalesce it with the
        # raw_file_samples rename used by every other row.
        if "sample_name" in main_df.columns and "sample_name" not in columns_to_keep:
            columns_to_keep.append("sample_name")
        if not plate_experiment:
            rt_sources = {key for key, name in PROP_MAP_BIOHUB.items() if name == "RT_index"}
            columns_to_keep = [c for c in columns_to_keep if c not in rt_sources]
        biohub_df = main_df[columns_to_keep].copy()
        biohub_df.rename(columns=PROP_MAP_BIOHUB, inplace=True)
        biohub_df = strip_author_metadata_column_prefix(biohub_df)
        biohub_df = collapse_duplicate_columns(biohub_df)
        # A file pooling samples joins their summaries with _join_unique, which
        # leaves 'a; a; b' when one sample already carried 'a; b'
        col = "experimental_condition_summary"
        if col in biohub_df.columns:
            biohub_df[col] = biohub_df[col].map(self._join_distinct_summaries)

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
        if "disease_ontology_term_id" not in biohub_df.columns:
            biohub_df["disease_ontology_term_id"] = "PATO:0000461"
        else:
            biohub_df["disease_ontology_term_id"] = biohub_df["disease_ontology_term_id"].apply(
                lambda v: pd.NA if (v is None or v == "" or v == [] or v == ()) else v
            )
            biohub_df["disease_ontology_term_id"] = biohub_df["disease_ontology_term_id"].fillna(
                "PATO:0000461"
            )

        if "self_reported_ethnicity_ontology_term_id" not in biohub_df.columns:
            biohub_df["self_reported_ethnicity_ontology_term_id"] = np.where(
                biohub_df["organism_ontology_term_id"] == "NCBITaxon:9606", "unknown", "na"
            )
        else:
            biohub_df.loc[
                biohub_df["self_reported_ethnicity_ontology_term_id"].isna()
                & (biohub_df["organism_ontology_term_id"] == "NCBITaxon:9606"),
                "self_reported_ethnicity_ontology_term_id",
            ] = "unknown"
            biohub_df.loc[
                biohub_df["self_reported_ethnicity_ontology_term_id"].isna()
                & (biohub_df["organism_ontology_term_id"] != "NCBITaxon:9606"),
                "self_reported_ethnicity_ontology_term_id",
            ] = "na"

        if "preservation_method" not in biohub_df.columns:
            biohub_df["preservation_method"] = "unknown"
        else:
            biohub_df["preservation_method"] = biohub_df["preservation_method"].apply(
                lambda v: pd.NA if (v is None or v == "" or v == [] or v == ()) else v
            )
            biohub_df["preservation_method"] = biohub_df["preservation_method"].fillna("unknown")

        if "experimental_condition_ontology_term_id" in biohub_df.columns:
            biohub_df["experimental_condition_ontology_term_id"] = biohub_df[
                "experimental_condition_ontology_term_id"
            ].apply(lambda v: pd.NA if (v is None or v == "" or v == [] or v == ()) else v)
            biohub_df["experimental_condition_ontology_term_id"] = biohub_df[
                "experimental_condition_ontology_term_id"
            ].fillna("na")

        for col in BIOHUB_SORT_ONTOLOGY_IDS:
            if col in biohub_df.columns:
                biohub_df = sort_ontology_term_id_column(biohub_df, col)

        # Combine multiple columns into one
        biohub_df = combine_bound_columns(
            biohub_df,
            lower_col="tissues_lower_bound_age",
            upper_col="tissues_upper_bound_age",
            units_col="tissues_age_units",
            out_col="age",
        )

        # Update values to match schema
        if "sex_ontology_term_id" in biohub_df.columns:
            biohub_df["sex_ontology_term_id"] = biohub_df["sex_ontology_term_id"].apply(
                self._map_sex_ontology_term_id
            )
        if "tissue_type" in biohub_df.columns:
            biohub_df["tissue_type"] = biohub_df["tissue_type"].apply(
                lambda x: TISSUE_TYPE_MAP.get(x[0], x) if isinstance(x, (list, tuple)) else x
            )
            cell_line = biohub_df["tissue_type"] == "cell line"
            for col in (
                "development_stage_ontology_term_id",
                "donor_id",
                "sex_ontology_term_id",
                "self_reported_ethnicity_ontology_term_id",
            ):
                if col not in biohub_df.columns:
                    biohub_df[col] = pd.NA
                biohub_df.loc[cell_line, col] = "na"
        if "development_stage_ontology_term_id" in biohub_df.columns:
            biohub_df["development_stage_ontology_term_id"] = (
                biohub_df["development_stage_ontology_term_id"].replace("", pd.NA).fillna("unknown")
            )
        if "genetic_perturbation_strategy" in biohub_df.columns:
            biohub_df["genetic_perturbation_strategy"] = biohub_df[
                "genetic_perturbation_strategy"
            ].replace(GENETIC_PERTURBATION_MAP)
        for field in REFORMAT_LIST:
            biohub_df = join_sequence_column(biohub_df, field)
        biohub_df = join_sequence_column(biohub_df, "RT_index", sep=",")
        if plate_experiment and "sample_probe_barcode" in biohub_df.columns:
            biohub_df = biohub_df.rename(columns={"sample_probe_barcode": "hash_index"})

        return biohub_df

    @staticmethod
    def _map_sex_ontology_term_id(value):
        """Map donor sex to a BIOHUB ontology term id; multi-value cells become unknown."""
        if is_empty(value):
            return value
        if isinstance(value, (list, tuple, set)):
            parts = [str(item).strip() for item in value if not is_empty(item)]
        else:
            parts = [part.strip() for part in str(value).split("; ") if part.strip()]
        if len(parts) > 1:
            return "unknown"
        if len(parts) == 1:
            return SEX_ONTOLOGY_MAP.get(parts[0], parts[0])
        return value

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
    def _source_titles(value):
        """Collect title values from an embedded Source object or list of objects."""
        if value is None or value == []:
            return []
        items = value if isinstance(value, list) else [value]
        titles = []
        for item in items:
            if isinstance(item, dict):
                title = item.get("title")
                if title not in (None, ""):
                    titles.append(title)
        return titles

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

    def _combine_sample_values(self, values):
        """
        One column's values across the samples sharing a raw matrix file.

        Mirrors how _flatten_resolved_references combines a field across several
        referenced objects: controlled terms stay dicts, because _join_unique()
        would str() them into the cell and split_controlled_term_columns() could
        no longer split them into _term_id and _term_name.

        A single distinct value passes through untouched - first, before any
        combining - so a file carrying one sample, which is every file in a run
        that is not multiplexed, keeps exactly the cell it had before: the same
        type, and terms in their original order rather than sorted by id.
        """
        present = [value for value in values if not is_empty(value)]
        if not present:
            return None

        unique = []
        for value in present:
            if value not in unique:
                unique.append(value)
        if len(unique) == 1:
            return unique[0]

        items = [item for value in unique for item in to_items(value)]
        if any(isinstance(item, dict) for item in items):
            # Only controlled terms can be deduped by id; anything else dict-
            # shaped, such as an embedded 'sources', is kept as it is rather
            # than dropped for having no '@id'
            if all(isinstance(item, dict) and item.get("@id") for item in items):
                return self._dedupe_terms(items)
            return items
        return self._join_unique(unique)

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
