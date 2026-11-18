"""Command-line entry for downloading a MatrixFileSet's RawMatrixFiles."""

import argparse
import sys

from db2_flattener import __version__
from db2_flattener.download.downloader import download_matrix_file_set
from db2_flattener.gather.lattice import Connection


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="db2-flattener-download",
        description="Download all RawMatrixFiles of a MatrixFileSet from S3. Files with a "
        "geo_file_name are saved under that name; others keep their S3 file name. "
        "Files that already exist in the output directory are not downloaded.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--uuid", "-u", required=True, help="UUID of the MatrixFileSet to download")
    parser.add_argument(
        "--out-directory", "-o", required=True, help="Directory to download files into"
    )
    parser.add_argument("--mode", "-m", help="mode/instance to run on", default="db2_demo")
    args = parser.parse_args()
    print(f"{parser.prog} {__version__}")

    try:
        connection = Connection(args.mode)
        downloaded, skipped, failed = download_matrix_file_set(
            args.uuid, args.out_directory, connection
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(
        f"\nDownloaded {len(downloaded)}, skipped {len(skipped)} already present, "
        f"failed {len(failed)} → {args.out_directory}"
    )
    if failed:
        sys.exit(1)
