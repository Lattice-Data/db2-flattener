"""Command-line entry for flattening a MatrixFileSet."""

import argparse
import sys

from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.gather.lattice import Connection
from db2_flattener.schema.constants import Configs
from db2_flattener.schema.generate import load_and_return_constant_dicts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Flatten MatrixFileSet data for DB2 processing and export to CSV"
    )
    parser.add_argument(
        "--uuid",
        "-u",
        required=True,
        help="UUID of the MatrixFileSet to process",
    )
    parser.add_argument(
        "--mode",
        "-m",
        help="mode/instance to run on",
        default="db2_demo",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output path prefix (optional). Writes {prefix}_MAIN.csv, "
        "{prefix}_BIOHUB.csv, {prefix}_GEO.csv, {prefix}_SAMPLES.csv, and "
        "{prefix}_GUIDE_METADATA.csv. Defaults to MatrixFileSet_{uuid}_{timestamp}",
    )
    args = parser.parse_args()

    try:
        connection = Connection(args.mode)
        field_types, object_config = load_and_return_constant_dicts(args.mode)
        configs = Configs(
            FIELD_TYPES=field_types,
            OBJECT_CONFIG=object_config,
        )
        flattener = DB2Flattener(connection, configs)
        output_file = flattener.flatten_matrix_file_set(args.uuid, args.output)

        if output_file:
            print(f"\nSuccess! CSV file created: {output_file}")
        else:
            print(f"Failed to process MatrixFileSet {args.uuid}")
            sys.exit(1)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
