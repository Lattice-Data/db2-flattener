"""Command-line entry for flattening a MatrixFileSet."""

import argparse
import sys

from db2_flattener import __version__
from db2_flattener.flatten.biohub_version import latest_biohub_versions_from_s3
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
    parser.add_argument(
        "--s3-order-dir",
        "-s",
        required=True,
        help="S3 URI of the order subdirectory. Each BIOHUB template's v# is "
        "taken independently per type and set from "
        "{samples|guide}_template_set{N}_v{M}.csv names in metadata/ under this "
        "path (none found for that type/set → v1; else highest v# + 1).",
    )
    audit_group = parser.add_mutually_exclusive_group()
    audit_group.add_argument(
        "--force-all-audit",
        action="store_true",
        help="Continue flattening even when ERROR or WARNING audits are present.",
    )
    audit_group.add_argument(
        "--force-only-warning",
        action="store_true",
        help="Continue flattening when only WARNING audits are present; "
        "still stop on ERROR audits.",
    )
    args = parser.parse_args()
    print(f"{parser.prog} {__version__}")

    try:
        biohub_versions = latest_biohub_versions_from_s3(args.s3_order_dir)
        if biohub_versions:
            found = ", ".join(
                f"{kind}_template_set{set_n}_v{version}"
                for (kind, set_n), version in sorted(biohub_versions.items())
            )
            print(f"Latest BIOHUB templates in metadata/: {found}")
        else:
            print("No BIOHUB templates in metadata/; all templates start at v1")

        if args.force_all_audit:
            audit_mode = "all_ok"
        elif args.force_only_warning:
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
            biohub_versions=biohub_versions,
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
