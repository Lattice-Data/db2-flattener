"""Command-line entry for flattening a MatrixFileSet."""

import argparse
import sys

from db2_flattener import __version__
from db2_flattener.flatten.biohub_version import resolve_biohub_version
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
    version_group = parser.add_mutually_exclusive_group(required=True)
    version_group.add_argument(
        "--s3-order-dir",
        "-s",
        help="S3 URI of the order subdirectory. BIOHUB v# is taken from "
        "metadata/ under this path (missing/empty → v1; else highest found v# + 1).",
    )
    version_group.add_argument(
        "--biohub-version",
        "-b",
        help="Explicit BIOHUB template version for this run (e.g. 1 or v1). "
        "Mutually exclusive with --s3-order-dir.",
    )
    audit_group = parser.add_mutually_exclusive_group()
    audit_group.add_argument(
        "--force_ALL_AUDIT",
        action="store_true",
        help="Continue flattening even when ERROR or WARNING audits are present.",
    )
    audit_group.add_argument(
        "--force_ONLY_WARNING",
        action="store_true",
        help="Continue flattening when only WARNING audits are present; "
        "still stop on ERROR audits.",
    )
    args = parser.parse_args()
    print(f"{parser.prog} {__version__}")

    try:
        biohub_version = resolve_biohub_version(
            s3_order_dir=args.s3_order_dir,
            biohub_version=args.biohub_version,
        )
        print(f"BIOHUB template version: {biohub_version}")

        if args.force_ALL_AUDIT:
            audit_mode = "all_ok"
        elif args.force_ONLY_WARNING:
            audit_mode = "warnings_ok"
        else:
            audit_mode = "strict"

        connection = Connection(args.mode)
        field_types, object_config = load_and_return_constant_dicts(args.mode)
        configs = Configs(
            FIELD_TYPES=field_types,
            OBJECT_CONFIG=object_config,
        )
        flattener = DB2Flattener(connection, configs)
        output_dir = flattener.flatten_matrix_file_set(
            args.uuid,
            biohub_version=biohub_version,
            audit_mode=audit_mode,
        )

        if output_dir:
            print(f"\nSuccess! Output directory created: {output_dir}")
        else:
            print(f"Failed to process MatrixFileSet {args.uuid}")
            sys.exit(1)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
