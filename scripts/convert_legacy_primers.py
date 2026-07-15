#!/usr/bin/env python3
"""Convert a legacy primer JSON file to the normalized primer database format."""

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import primer_checker  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description="Convert legacy primer JSON to normalized primer database JSON.")
    parser.add_argument("--input", required=True, help="Legacy primer JSON path.")
    parser.add_argument("--output", required=True, help="Normalized JSON output path.")
    parser.add_argument("--database-version", required=True, help="Database version to write, e.g. 2026-05-06.")
    parser.add_argument("--source", default="legacy", help="Source label to store in scheme metadata.")
    args = parser.parse_args()

    with open(args.input, encoding="utf-8") as infile:
        legacy = json.load(infile)

    validation = primer_checker.validate_legacy_primer_library(legacy)
    if validation.errors:
        for error in validation.errors:
            print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
    for warning in validation.warnings:
        print(f"Warning: {warning}", file=sys.stderr)

    normalized = primer_checker.legacy_library_to_normalized_database(
        legacy,
        database_version=args.database_version,
        source=args.source,
    )
    normalized_validation = primer_checker.validate_normalized_primer_library(normalized)
    if normalized_validation.errors:
        for error in normalized_validation.errors:
            print(f"Converted database error: {error}", file=sys.stderr)
        sys.exit(1)
    for warning in normalized_validation.warnings:
        print(f"Converted database warning: {warning}", file=sys.stderr)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(normalized, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"Wrote normalized primer database to {output_path}")


if __name__ == "__main__":
    main()
