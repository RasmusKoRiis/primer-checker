#!/usr/bin/env python3
"""
Run primer analysis over a folder of FASTA files with targets inferred from filenames.

Examples:
    python3 scripts/run_primer_checker_batch.py --input-folder data --primers primer_db/fhi_primers.normalized.json
    python3 scripts/run_primer_checker_batch.py --input-folder data --analysis-only
    python3 scripts/run_primer_checker_batch.py --input-folder data --dry-run
"""

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import primer_analysis  # noqa: E402
import primer_report  # noqa: E402


FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}


def default_primer_db() -> str:
    unified = REPO_ROOT / "primer_db" / "fhi_primers.unified.json"
    normalized = REPO_ROOT / "primer_db" / "fhi_primers.normalized.json"
    if unified.exists():
        return str(unified)
    return str(normalized if normalized.exists() else REPO_ROOT / "primers.json")


def find_fasta_files(input_folder: Path, recursive: bool = False) -> list[Path]:
    pattern = "**/*" if recursive else "*"
    return sorted(
        path for path in input_folder.glob(pattern)
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def describe_target(target: primer_analysis.AnalysisTarget | None) -> str:
    if target is None:
        return "unclassified"
    if target.virus_type.lower() == "influenza":
        return f"influenza --flu-type {target.flu_type}"
    return target.virus_type


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run primer analysis for every FASTA in a folder, choosing virus/subtype from each filename."
    )
    parser.add_argument("--input-folder", required=True, help="Folder containing FASTA files.")
    parser.add_argument(
        "--primers",
        default=default_primer_db(),
        help="Primer JSON file. Defaults to primer_db/fhi_primers.normalized.json when present.",
    )
    parser.add_argument("--output", default="batch_primer_report.csv", help="Combined output CSV path.")
    parser.add_argument(
        "--html-report",
        default="batch_primer_report.html",
        help="Combined self-contained HTML report path. Defaults to batch_primer_report.html.",
    )
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Only run primer analysis and write the CSV; skip HTML report generation.",
    )
    parser.add_argument(
        "--previous-report-csv",
        action="append",
        default=[],
        help="Optional previous CSV report to embed in the HTML report. Repeat for multiple reports.",
    )
    parser.add_argument(
        "--metadata-csv",
        help="Optional sample metadata CSV with SampleID plus optional Sample_Date and Ct/Ct_Value columns.",
    )
    parser.add_argument("--recursive", action="store_true", help="Search the input folder recursively.")
    parser.add_argument("--dry-run", action="store_true", help="Print inferred targets without running BLAST.")
    parser.add_argument(
        "--fail-on-unclassified",
        action="store_true",
        help="Exit with an error if any FASTA filename cannot be classified.",
    )
    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    if not input_folder.is_dir():
        sys.exit(f"Input folder does not exist or is not a directory: {input_folder}")

    fasta_files = find_fasta_files(input_folder, recursive=args.recursive)
    if not fasta_files:
        sys.exit(f"No FASTA files found in {input_folder}")

    primer_records, validation = primer_analysis.load_primer_records(args.primers)
    primer_analysis.print_validation_messages(validation)
    metadata_records, metadata_validation = primer_analysis.load_metadata_csv(args.metadata_csv)
    if metadata_validation.errors:
        error_text = "\n".join(f"- {error}" for error in metadata_validation.errors)
        sys.exit(f"Metadata CSV validation failed:\n{error_text}")
    primer_analysis.print_metadata_validation_messages(metadata_validation)

    planned: list[tuple[Path, primer_analysis.AnalysisTarget]] = []
    unclassified: list[Path] = []
    available_organisms = list(primer_records)
    for fasta_file in fasta_files:
        target = primer_analysis.infer_analysis_target_from_filename(
            str(fasta_file), available_organisms=available_organisms, primer_records=primer_records
        )
        print(f"{fasta_file}: {describe_target(target)}")
        if target is None:
            unclassified.append(fasta_file)
        else:
            planned.append((fasta_file, target))

    if unclassified and args.fail_on_unclassified:
        names = "\n".join(f"  - {path}" for path in unclassified)
        sys.exit(f"Could not classify these FASTA files from their filenames:\n{names}")

    if args.dry_run:
        return

    if not planned:
        sys.exit("No classified FASTA files to process.")

    primer_analysis.ensure_blastn_available()

    all_results = []
    for fasta_file, target in planned:
        virus_type, selected_primers = primer_analysis.select_primer_records(
            primer_records,
            target.virus_type,
            target.flu_type,
        )
        print(f"Processing FASTA file: {fasta_file} as {describe_target(target)}")
        all_results.extend(
            primer_analysis.process_fasta_file(
                str(fasta_file),
                virus_type,
                selected_primers,
                metadata_records=metadata_records,
            )
        )

    primer_analysis.write_csv_report(all_results, args.output)
    if not args.analysis_only:
        previous_reports = primer_report.load_previous_report_csvs(args.previous_report_csv)
        primer_report.write_html_report(all_results, args.html_report, previous_reports=previous_reports)


if __name__ == "__main__":
    main()
