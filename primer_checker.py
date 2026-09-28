#!/usr/bin/env python3
"""
Command-line entrypoint and compatibility imports for primer checker.

The implementation is split between:
  - primer_analysis.py: primer loading, metadata matching, FASTA/BLAST analysis,
    and CSV writing.
  - primer_report.py: previous-report loading and self-contained HTML report
    generation.
"""

import argparse
import os
import sys

from primer_analysis import *  # noqa: F401,F403
from primer_report import *  # noqa: F401,F403


def main():
    parser = argparse.ArgumentParser(
        description="Use BLASTn to check primer alignments (allowing partial alignments) for a given virus type."
    )
    parser.add_argument(
        "--primers",
        type=str,
        default="primers.json",
        help="Path to the primer library file (JSON format)."
    )
    parser.add_argument(
        "--virus",
        type=str,
        help="Virus type (e.g., 'SARS-CoV-2', 'Influenza', 'RSV-A', 'RSV-B')."
    )
    parser.add_argument(
        "--flu-type",
        type=str,
        help="Influenza type or subtype from the database, e.g. A (all A primers), H5N1, or B/VICTORIA."
    )
    parser.add_argument(
        "--assay-type",
        choices=["all", "pcr", "ngs"],
        default="all",
        help="Select PCR/qPCR schemes, NGS panels, or all assays (default: all).",
    )
    parser.add_argument(
        "--assay-id",
        help="Optional exact PCR scheme_id or NGS panel_id to analyze.",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        nargs="+",
        help="One or more FASTA file paths (subject sequences) to process."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="primer_report.csv",
        help="Output CSV filename (default: primer_report.csv)."
    )
    parser.add_argument(
        "--html-report",
        type=str,
        help="Optional self-contained HTML report filename."
    )
    parser.add_argument(
        "--previous-report-csv",
        action="append",
        default=[],
        help="Optional previous CSV report to embed in the HTML report. Repeat to attach multiple reports."
    )
    parser.add_argument(
        "--metadata-csv",
        type=str,
        help="Optional sample metadata CSV with SampleID plus optional Sample_Date and Ct/Ct_Value columns."
    )
    parser.add_argument(
        "--validate-primers",
        action="store_true",
        help="Validate the primer library and exit without running BLAST."
    )
    args = parser.parse_args()

    primer_records, validation = load_primer_records(args.primers)
    print_validation_messages(validation)
    if load_primer_library(args.primers).get("purpose") == "synthetic-test-only":
        print("DUMMY DATABASE: synthetic software test data only.", file=sys.stderr)
    if args.validate_primers:
        print(f"Primer library validation succeeded for {args.primers}")
        return

    if not args.virus:
        sys.exit("Error: please supply --virus, or use --validate-primers to only validate the primer library.")

    if not args.fasta:
        sys.exit("Error: please supply one or more FASTA files with --fasta, or use --validate-primers.")

    virus_type, selected_primers = select_primer_records(
        primer_records,
        args.virus,
        args.flu_type,
        assay_type=args.assay_type,
        assay_id=args.assay_id,
    )
    metadata_records, metadata_validation = load_metadata_csv(args.metadata_csv)
    if metadata_validation.errors:
        error_text = "\n".join(f"- {error}" for error in metadata_validation.errors)
        sys.exit(f"Metadata CSV validation failed:\n{error_text}")
    print_metadata_validation_messages(metadata_validation)
    ensure_blastn_available()

    all_results = []
    for fasta_file in args.fasta:
        if not os.path.exists(fasta_file):
            print(f"FASTA file '{fasta_file}' not found. Skipping.", file=sys.stderr)
            continue
        print(f"Processing FASTA file: {fasta_file}")
        file_results = process_fasta_file(fasta_file, virus_type, selected_primers, metadata_records=metadata_records)
        all_results.extend(file_results)

    write_csv_report(all_results, args.output)
    if args.html_report:
        previous_reports = load_previous_report_csvs(args.previous_report_csv)
        write_html_report(all_results, args.html_report, previous_reports=previous_reports)


if __name__ == "__main__":
    main()
