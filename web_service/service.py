"""Upload validation and orchestration; all analytical decisions live in primer_analysis."""

import csv
import hashlib
import html
import io
import json
import os
import re
import subprocess
import tempfile
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import primer_analysis as engine
import primer_report

ROOT = Path(__file__).resolve().parents[1]
APP_VERSION = "0.1.0"
MAX_UPLOAD_BYTES = 3_000_000
MAX_REQUEST_BYTES = 4_000_000
MAX_RESPONSE_BYTES = 4_000_000
MAX_FILES = 10
MAX_RECORDS = 200
MAX_COMPARISONS = 2000
MAX_BLAST_CALLS = 1200
ANALYSIS_SECONDS = 240


class WebError(Exception):
    def __init__(self, message: str, code: str = "invalid_input", status: int = 422):
        super().__init__(message)
        self.message, self.code, self.status = message, code, status


def database_path() -> Path:
    return Path(os.environ.get("PRIMER_DATABASE_PATH", ROOT / "primer_db/fhi_primers.unified.json"))


def load_database():
    try:
        path = database_path()
        records, validation = engine.load_primer_records(str(path))
        # Include loaded primer sequences/metadata, including BED/FASTA assets, in the fingerprint.
        # source_file is installation-dependent, so exclude it from reproducibility fingerprints.
        canonical = json.dumps(
            [{k: v for k, v in asdict(p).items() if k != "source_file"} for group in records.values() for p in group],
            sort_keys=True,
        )
        versions = sorted({p.database_version for group in records.values() for p in group if p.database_version})
        return (
            records,
            {
                "version": ", ".join(versions) or "unversioned",
                "sha256": hashlib.sha256(canonical.encode()).hexdigest(),
            },
            validation,
        )
    except (SystemExit, Exception):  # noqa: BLE001 - never expose local database paths
        raise WebError(
            "The primer database could not be loaded. Contact the deployment maintainer.", "database_error", 503
        ) from None


def catalog() -> dict:
    records, database, _ = load_database()
    viruses = []
    if "Influenza-A" in records or "Influenza-B" in records:
        viruses.append(
            {
                "id": "influenza",
                "name": "Influenza",
                "subtypes": (["A", "H1", "H3"] if "Influenza-A" in records else [])
                + (["B"] if "Influenza-B" in records else []),
            }
        )
    viruses.extend(
        {"id": name, "name": name, "subtypes": []} for name in records if name not in {"Influenza-A", "Influenza-B"}
    )
    assays = {}
    for group in records.values():
        for p in group:
            key = (p.organism, p.scheme_id, p.assay_type)
            if key not in assays:
                assays[key] = {
                    "id": p.scheme_id,
                    "name": p.assay_name or p.scheme_id or "Default primers",
                    "organism": p.organism,
                    "type": p.assay_type,
                    "primers": 0,
                }
            assays[key]["primers"] += 1
    return {
        "application_version": APP_VERSION,
        "database": database,
        "viruses": viruses,
        "assays": list(assays.values()),
        "limits": {
            "upload_bytes": MAX_UPLOAD_BYTES,
            "files": MAX_FILES,
            "records": MAX_RECORDS,
            "comparisons": MAX_COMPARISONS,
        },
    }


def safe_filename(name: str, suffixes: set[str]) -> str:
    name = name.replace("\\", "/").split("/")[-1]
    name = re.sub(r"[^A-Za-z0-9._-]", "_", name).strip(".")[:100]
    if not name or Path(name).suffix.lower() not in suffixes:
        raise WebError("Choose FASTA files (.fasta, .fa, .fna, .fas) and a .csv metadata file.")
    return name


def decode_text(data: bytes, label: str) -> str:
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise WebError(f"{label} must be a UTF-8 text file.") from None
    if any(ord(c) < 32 and c not in "\n\r\t" for c in text):
        raise WebError(f"{label} contains invalid control characters.")
    return text


def validate_fasta(data: bytes) -> str:
    text = decode_text(data, "FASTA")
    seen, current, length = set(), None, 0
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None and not length:
                raise WebError("FASTA contains a header without a sequence.", "invalid_fasta")
            header = line[1:].split()
            if not header or len(header[0]) > 200 or len(line) > 1000:
                raise WebError("FASTA headers need an identifier of at most 200 characters.", "invalid_fasta")
            current, length = header[0], 0
            # BLAST treats pipe-delimited IDs specially; reject reserved parser prefixes.
            if current.startswith(("gi|", "lcl|", "ref|", "gb|")):
                raise WebError(
                    "Use plain FASTA identifiers rather than reserved BLAST ID prefixes (gi|, lcl|, ref|, gb|).",
                    "invalid_fasta",
                )
            if current in seen:
                raise WebError("FASTA identifiers must be unique within each file.", "invalid_fasta")
            seen.add(current)
            if len(seen) > MAX_RECORDS:
                raise WebError(f"Use at most {MAX_RECORDS} sequence records per analysis.", "work_limit", 413)
        else:
            if current is None or not set(line.upper()) <= engine.ALLOWED_SEQUENCE_CODES:
                raise WebError(
                    "Invalid FASTA: use a >header followed by nucleotide sequences with IUPAC bases. FASTQ and gapped alignments are not accepted.",
                    "invalid_fasta",
                )
            length += len(line)
    if not seen or not length:
        raise WebError("FASTA is empty or contains a header without a sequence.", "invalid_fasta")
    # The engine's two FASTA readers differ in whitespace handling. Normalize
    # surrounding whitespace once so validation and both readers see the same records.
    return "\n".join(line.strip() for line in text.splitlines()) + "\n"


def validate_metadata(data: bytes) -> str:
    text = decode_text(data, "Metadata CSV")
    try:
        rows = list(csv.reader(io.StringIO(text), strict=True))
    except csv.Error:
        raise WebError("Malformed metadata CSV. Check quoting and column delimiters.", "invalid_metadata") from None
    if not rows or len(rows) < 2 or len(rows) > 2001:
        raise WebError("Metadata CSV needs a header and 1–2000 data rows.", "invalid_metadata")
    header = rows[0]
    if (
        len(header) > 100
        or len({engine.normalize_column_name(c) for c in header}) != len(header)
        or any(not c.strip() for c in header)
    ):
        raise WebError("Metadata CSV needs unique, nonempty column names (at most 100).", "invalid_metadata")
    if not engine.find_metadata_column(header, {"sampleid", "sample", "id"}):
        raise WebError("Metadata CSV needs a SampleID, Sample_ID, Sample, or ID column.", "invalid_metadata")
    if any(len(row) != len(header) for row in rows[1:] if row):
        raise WebError("Every metadata CSV row must have the same number of columns as the header.", "invalid_metadata")
    if any(len(cell) > 500 for row in rows for cell in row):
        raise WebError("Metadata CSV cells must contain at most 500 characters.", "invalid_metadata")
    return text


def analyze(
    files: list[tuple[str, bytes]],
    metadata: tuple[str, bytes] | None,
    virus: str,
    flu_type: str | None,
    assay_type: str,
    assay_id: str | None,
) -> dict:
    started = time.monotonic()
    if not files or len(files) > MAX_FILES:
        raise WebError(f"Upload between 1 and {MAX_FILES} FASTA files.")
    if sum(len(data) for _, data in files) + (len(metadata[1]) if metadata else 0) > MAX_UPLOAD_BYTES:
        raise WebError(
            "Combined uploads exceed the 3 MB limit. Split the analysis into smaller batches.", "upload_too_large", 413
        )
    records, database, validation = load_database()
    if assay_type not in {"pcr", "ngs", "all"}:
        raise WebError("Select PCR, NGS, or All assays.", "invalid_selection")
    allowed_viruses = set(records) | {"influenza"}
    if virus not in allowed_viruses or (virus == "influenza" and flu_type not in {"A", "H1", "H3", "B"}):
        raise WebError("Select an available virus and an influenza subtype when applicable.", "invalid_selection")
    if virus != "influenza" and flu_type:
        raise WebError("Influenza subtype is only valid for Influenza.", "invalid_selection")
    try:
        selected_virus, primers = engine.select_primer_records(records, virus, flu_type, assay_type, assay_id)
    except SystemExit:
        raise WebError("No primers match this virus, subtype, and assay selection.", "invalid_selection") from None
    if len(primers) * len(files) > MAX_BLAST_CALLS:
        raise WebError(
            "This selection includes too many BLAST searches. Select one assay or fewer files.", "work_limit", 413
        )
    warnings = []
    if validation.warnings:
        warnings.append(
            "The primer database has validation warnings; the maintainer should review it with --validate-primers."
        )
    rows, inputs, total_records = [], [], 0
    with tempfile.TemporaryDirectory(prefix="primer-web-") as folder:
        root = Path(folder)
        metadata_records = []
        if metadata:
            safe_filename(metadata[0], {".csv"})
            metadata_path = root / "metadata.csv"
            metadata_path.write_text(validate_metadata(metadata[1]), encoding="utf-8")
            metadata_records, check = engine.load_metadata_csv(str(metadata_path))
            if check.errors:
                raise WebError(
                    "Metadata could not be loaded. Check its sample ID, date, and Ct columns.", "invalid_metadata"
                )
            if check.warnings:
                warnings.append(
                    "Some metadata fields are missing or sample IDs are repeated. Ambiguous metadata is not attached."
                )
        prepared, used_names = [], set()
        comparisons = 0
        for index, (name, data) in enumerate(files):
            name = safe_filename(name, {".fa", ".fas", ".fna", ".fasta"})
            if name in used_names:
                raise WebError("FASTA filenames must be distinct after removing unsafe characters.")
            used_names.add(name)
            path = root / str(index) / name
            path.parent.mkdir()
            path.write_text(validate_fasta(data), encoding="utf-8")
            subjects = engine.get_subject_ids(str(path))
            total_records += len(subjects)
            comparisons += sum(
                len(engine.filter_subject_ids_for_primer(subjects, p, selected_virus, str(path), quiet=True))
                for p in primers
            )
            if total_records > MAX_RECORDS or comparisons > MAX_COMPARISONS:
                raise WebError(
                    f"Use at most {MAX_RECORDS} sequence records and {MAX_COMPARISONS} primer/record comparisons. Select fewer files or one assay.",
                    "work_limit",
                    413,
                )
            if selected_virus.startswith("Influenza") and any(not engine.get_segment(s) for s in subjects):
                warnings.append(
                    f"{name}: some headers lack influenza segment tokens (for example 01-HA|sample). Segment-specific primers exclude those records."
                )
            prepared.append(path)
            inputs.append({"filename": name, "sha256": hashlib.sha256(data).hexdigest(), "records": len(subjects)})
        if not comparisons:
            raise WebError(
                "No sequence records match the selected primer segments. Influenza headers need segment tokens such as 01-HA|sample or 03-M|sample.",
                "no_comparisons",
            )
        executable = engine.resolve_blastn()
        try:
            blast_version = subprocess.run(
                [executable, "-version"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
                env=engine.blast_environment(executable),
            ).stdout.splitlines()[0]
        except (OSError, subprocess.SubprocessError, IndexError):
            raise engine.BlastUnavailableError("BLAST executable could not start.") from None
        execution = engine.BlastExecution(executable, started + ANALYSIS_SECONDS, strict_errors=True, quiet=True)
        for path in prepared:
            rows.extend(
                engine.process_fasta_file(str(path), selected_virus, primers, metadata_records, execution=execution)
            )

    # The temporary directory and uploads are gone before serializing the response.
    manifest = {
        "analysis_utc": datetime.now(timezone.utc).isoformat(),
        "application_version": APP_VERSION,
        "git_commit": os.environ.get("VERCEL_GIT_COMMIT_SHA", os.environ.get("APP_GIT_COMMIT", "unknown")),
        "database": database,
        "selection": {"virus": virus, "flu_type": flu_type, "assay_type": assay_type, "assay_id": assay_id},
        "files": inputs,
        "metadata_sha256": hashlib.sha256(metadata[1]).hexdigest() if metadata else None,
        "blast": {
            "version": blast_version,
            "reward": 2,
            "penalty": -3,
            "word_size": 4,
            "dust": "yes",
            "max_target_seqs": engine.BLAST_MAX_TARGET_SEQS,
        },
    }
    sample_key = lambda r: (
        ("metadata", r["Metadata_Sample_ID"])
        if r["Metadata_Sample_ID"]
        else (r["Fasta_File"], r["Subject_Sequence_ID"])
    )
    affected = [r for r in rows if r["Hit_Status"] == "hit" and r["Mismatches"] > 0]
    summary = {
        "files": len(files),
        "sequence_records": total_records,
        "samples": len({sample_key(r) for r in rows}),
        "primers": len(primers),
        "comparisons": len(rows),
        "hits": sum(r["Hit_Status"] == "hit" for r in rows),
        "mismatch_comparisons": len(affected),
        "samples_affected": len({sample_key(r) for r in affected}),
    }
    csv_stream = io.StringIO(newline="")
    # Protect web CSV downloads against spreadsheet formula injection without changing raw results/CLI.
    csv_rows = [
        {
            k: ("'" + v if isinstance(v, str) and v.lstrip().startswith(("=", "+", "-", "@")) else v)
            for k, v in row.items()
        }
        for row in rows
    ]
    engine.write_csv_rows(csv_rows, csv_stream)
    report_html = primer_report.build_html_report(rows)
    manifest_json = json.dumps(manifest, indent=2)
    provenance = (
        '<details style="margin:24px"><summary>Analysis provenance</summary><pre>'
        + html.escape(manifest_json)
        + "</pre></details>"
    )
    report_html = report_html.replace("</body>", provenance + "</body>")
    return {
        "rows": rows,
        "summary": summary,
        "manifest": manifest,
        "warnings": list(dict.fromkeys(warnings)),
        "downloads": {"csv": csv_stream.getvalue(), "html": report_html, "manifest": manifest_json},
    }
