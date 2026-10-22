"""Command-line entry for flattening a MatrixFileSet."""

import argparse
import sys

from db2_flattener import __version__
from db2_flattener.flatten.flattener import DB2Flattener
from db2_flattener.gather.lattice import Connection
from db2_flattener.schema.constants import Configs
from db2_flattener.schema.generate import load_and_return_constant_dicts


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="db2-flattener",
        description="Flatten MatrixFileSet data for DB2 processing and export to CSV",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
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
    args = parser.parse_args()
    print(f"{parser.prog} {__version__}")

    try:
        connection = Connection(args.mode)
        field_types, object_config = load_and_return_constant_dicts(args.mode)
        configs = Configs(
            FIELD_TYPES=field_types,
            OBJECT_CONFIG=object_config,
        )
        flattener = DB2Flattener(connection, configs)
        output_dir = flattener.flatten_matrix_file_set(args.uuid)

        if output_dir:
            print(f"\nSuccess! Output directory created: {output_dir}")
        else:
            print(f"Failed to process MatrixFileSet {args.uuid}")
            sys.exit(1)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
