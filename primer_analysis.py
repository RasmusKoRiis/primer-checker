#!/usr/bin/env python3
"""Primer loading, metadata matching, FASTA parsing, BLAST analysis, and CSV output."""

import csv
from dataclasses import dataclass, field
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# --- Ambiguous Nucleotide Handling ---
AMBIGUITY_CODES = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'U': {'T'},  # Treat U as T if needed
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}
ALLOWED_SEQUENCE_CODES = set(AMBIGUITY_CODES)
INFLUENZA_SEGMENTS = {"HA", "M", "NS"}
CSV_FIELDNAMES = [
    "Fasta_File",
    "Virus_Type",
    "Assay_Type",
    "Assay_ID",
    "Assay_Name",
    "Primer_Name",
    "Primer_Sequence",
    "Primer_Segment",
    "Primer_Pool",
    "Primer_Start",
    "Primer_End",
    "Subject_Sequence_ID",
    "Subject_Segment",
    "Hit_Status",
    "Percent_Identity",
    "Alignment_Length",
    "Mismatches",
    "Gap_Openings",
    "Query_Start",
    "Query_End",
    "Subject_Start",
    "Subject_End",
    "E_value",
    "Bitscore",
    "Mismatch_Positions",
    "Mismatch_Details",
    "Query_Alignment",
    "Subject_Alignment",
    "Metadata_Sample_ID",
    "Sample_Date",
    "Ct_Value",
    "Ct_Source",
]
BLAST_MAX_TARGET_SEQS = "5000"


class BlastError(RuntimeError):
    """BLAST could not complete; distinct from a successful search with no hits."""


class BlastUnavailableError(BlastError):
    pass


class AnalysisTimeoutError(BlastError):
    pass


@dataclass(frozen=True)
class BlastExecution:
    """Optional web execution policy; does not change alignment parameters."""
    executable: str | None = None
    deadline: float | None = None
    strict_errors: bool = False
    quiet: bool = False


def resolve_blastn() -> str:
    configured = os.environ.get("BLASTN_PATH")
    bundled = Path(__file__).resolve().parent / "bin" / "blastn"
    candidate = configured or (str(bundled) if bundled.is_file() else "blastn")
    resolved = shutil.which(candidate)
    if not resolved:
        raise BlastUnavailableError("BLAST executable unavailable. Install BLAST+ or configure BLASTN_PATH.")
    return resolved


def blast_environment(executable: str) -> dict[str, str] | None:
    """Resolve adjacent packaged libraries without changing process-wide state."""
    libraries = Path(executable).resolve().parent / "lib"
    if sys.platform != "linux" or not libraries.is_dir():
        return None
    existing = os.environ.get("LD_LIBRARY_PATH", "")
    return {**os.environ, "LD_LIBRARY_PATH": str(libraries) + (os.pathsep + existing if existing else "")}


@dataclass(frozen=True)
class PrimerRecord:
    """Internal primer representation used by legacy and future database formats."""
    organism: str
    name: str
    sequence: str
    segment: str = ""
    role: str = "primer"
    subtype_tags: tuple[str, ...] = field(default_factory=tuple)
    scheme_id: str = ""
    scheme_version: str = ""
    database_version: str = ""
    pool: str = ""
    strand: str = ""
    reference_name: str = ""
    start: str = ""
    end: str = ""
    source_file: str = ""
    assay_type: str = "pcr"
    assay_name: str = ""


@dataclass(frozen=True)
class AnalysisTarget:
    """Virus/subtype selection inferred for one input FASTA file."""
    virus_type: str
    flu_type: str | None = None
    reason: str = ""


@dataclass(frozen=True)
class MetadataRecord:
    """Optional per-sample metadata attached to result rows."""
    sample_id: str
    sample_date: str = ""
    ct_value: str = ""
    ct_values: dict[str, str] = field(default_factory=dict)


@dataclass
class ValidationResult:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

def bases_match(base1: str, base2: str) -> bool:
    """Returns True if the two nucleotide bases match, taking ambiguous IUPAC codes into account."""
    base1 = base1.upper()
    base2 = base2.upper()
    set1 = AMBIGUITY_CODES.get(base1, {base1})
    set2 = AMBIGUITY_CODES.get(base2, {base2})
    return bool(set1 & set2)

def count_mismatches(query_aln: str, subject_aln: str) -> int:
    """
    Counts mismatches between two aligned sequences (query and subject),
    using ambiguous nucleotide matching. Gaps ('-') are treated as mismatches.
    Both input strings should be of equal length.
    """
    mismatches = 0
    for a, b in zip(query_aln, subject_aln):
        if a == '-' or b == '-' or not bases_match(a, b):
            mismatches += 1
    return mismatches


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
    return sequence.translate(complement)[::-1].upper()


def read_fasta_sequences(fasta_file: str) -> dict[str, str]:
    """Read FASTA records into an ID -> uppercase sequence mapping."""
    sequences = {}
    current_id = ""
    chunks = []
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(chunks).upper()
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id:
        sequences[current_id] = "".join(chunks).upper()
    return sequences


def _pad_left(sequence: str, expected_length: int) -> str:
    return ("-" * max(0, expected_length - len(sequence))) + sequence


def _pad_right(sequence: str, expected_length: int) -> str:
    return sequence + ("-" * max(0, expected_length - len(sequence)))


def reconstruct_full_primer_alignment(
    primer_seq: str,
    qseq: str,
    sseq: str,
    qstart: int,
    qend: int,
    sstart: int,
    send: int,
    subject_sequence: str,
) -> tuple[str, str]:
    """
    Extend a local BLAST primer hit to primer length using adjacent subject bases.

    BLAST optimizes a local alignment, so it may omit low-scoring primer-end
    bases. For reporting, those omitted bases are important because terminal
    mismatches can affect primer performance.
    """
    missing_start = qstart - 1
    missing_end = len(primer_seq) - qend
    if not subject_sequence:
        return primer_seq, (
            ("-" * missing_start)
            + sseq.upper()
            + ("-" * missing_end)
        )

    subject_sequence = subject_sequence.upper()
    plus_strand = sstart <= send
    if plus_strand:
        left_start = max(0, sstart - missing_start - 1)
        left = subject_sequence[left_start:max(0, sstart - 1)]
        left = _pad_left(left, missing_start)
        right = subject_sequence[send:send + missing_end]
        right = _pad_right(right, missing_end)
    else:
        left = reverse_complement(subject_sequence[sstart:sstart + missing_start])
        left = _pad_left(left, missing_start)
        right_start = max(0, send - missing_end - 1)
        right = reverse_complement(subject_sequence[right_start:max(0, send - 1)])
        right = _pad_right(right, missing_end)

    return primer_seq.upper(), left + sseq.upper() + right


def get_mismatch_positions(qseq: str, sseq: str, qstart: int, qend: int, full_primer: str) -> str:
    """
    Determines the positions (with respect to the full primer, 1-indexed)
    at which mismatches occur. In the aligned region (qseq vs sseq),
    positions where bases do not match (using ambiguous matching) are noted.
    Additionally, any bases in the primer not covered by the alignment (i.e.
    positions before qstart or after qend) are included.
    Returns a comma-separated string of mismatch positions.
    """
    mismatches = []

    # Check aligned region positions.
    for i, (q_base, s_base) in enumerate(zip(qseq, sseq)):
        pos = qstart + i  # position in full primer (1-indexed)
        if q_base == '-' or s_base == '-' or not bases_match(q_base, s_base):
            mismatches.append(pos)

    # Include positions for missing bases at the beginning.
    for pos in range(1, qstart):
        mismatches.append(pos)
    # Include positions for missing bases at the end.
    for pos in range(qend+1, len(full_primer)+1):
        mismatches.append(pos)

    mismatches = sorted(mismatches)
    return ",".join(map(str, mismatches)) if mismatches else ""


def get_mismatch_details(qseq: str, sseq: str, qstart: int, qend: int, full_primer: str) -> str:
    """
    Return comma-separated mismatch details as position:primer_base>subject_base.
    Bases outside partial BLAST alignment are represented as subject gaps.
    """
    details = []
    for i, (q_base, s_base) in enumerate(zip(qseq, sseq)):
        pos = qstart + i
        if q_base == '-' or s_base == '-' or not bases_match(q_base, s_base):
            primer_base = full_primer[pos - 1] if 1 <= pos <= len(full_primer) else q_base
            details.append((pos, f"{pos}:{primer_base.upper()}>{s_base.upper()}"))

    for pos in range(1, qstart):
        details.append((pos, f"{pos}:{full_primer[pos - 1].upper()}>-"))
    for pos in range(qend + 1, len(full_primer) + 1):
        details.append((pos, f"{pos}:{full_primer[pos - 1].upper()}>-"))

    return ",".join(detail for _pos, detail in sorted(details))

# --- End Ambiguity Functions ---


def tokenize_name(name: str) -> list[str]:
    """Split primer names into comparison tokens without losing established naming support."""
    return [token for token in re.split(r"[^A-Za-z0-9]+", name.upper()) if token]


def infer_primer_role(primer_name: str) -> str:
    tokens = tokenize_name(primer_name)
    return "probe" if any("PROBE" in token or token in {"TM", "VIC2", "YAM2"} for token in tokens) else "primer"


def infer_primer_segment(organism: str, primer_name: str) -> str:
    """Infer a legacy primer segment from name tokens for backward-compatible databases."""
    if not organism.upper().startswith("INFLUENZA"):
        return ""
    tokens = tokenize_name(primer_name)
    for segment in ("HA", "NS", "M"):
        if segment in tokens:
            return segment
    return ""


def infer_subtype_tags(primer_name: str) -> tuple[str, ...]:
    tokens = tokenize_name(primer_name)
    return tuple(tag for tag in ("H1", "H3") if tag in tokens)


def infer_analysis_target_from_filename(fasta_file: str, available_organisms: list[str] | None = None) -> AnalysisTarget | None:
    """
    Infer the primer target from a FASTA filename.

    This is intentionally conservative: H1/H3 influenza subtype rules are applied
    before broader Influenza-A rules, and ambiguous RSV names without A/B are left
    unclassified because the primer database stores RSV-A and RSV-B separately.
    """
    stem = os.path.splitext(os.path.basename(fasta_file))[0]
    upper_stem = stem.upper()
    tokens = set(tokenize_name(stem))
    compact = re.sub(r"[^A-Z0-9]+", "", upper_stem)

    has_h3 = "H3" in tokens or "H3N2" in tokens or bool(re.search(r"(^|[^A-Z0-9])H3(N2)?($|[^A-Z0-9])", upper_stem))
    has_h1 = "H1" in tokens or "H1N1" in tokens or bool(re.search(r"(^|[^A-Z0-9])H1(N1)?($|[^A-Z0-9])", upper_stem))
    if has_h3:
        return AnalysisTarget("influenza", "H3", "filename contains H3/H3N2")
    if has_h1:
        return AnalysisTarget("influenza", "H1", "filename contains H1/H1N1")

    if tokens & {"INFB", "FLUB"} or compact in {"INFB", "FLUB", "INFLUENZAB"} or compact.startswith(("INFB", "FLUB")):
        return AnalysisTarget("influenza", "B", "filename contains Influenza-B marker")
    if tokens & {"INFA", "FLUA"} or compact in {"INFA", "FLUA", "INFLUENZAA"} or compact.startswith(("INFA", "FLUA")):
        return AnalysisTarget("influenza", "A", "filename contains Influenza-A marker")

    if tokens & {"RSVA", "RSV-A"} or compact.startswith("RSVA"):
        return AnalysisTarget("RSV-A", None, "filename contains RSV-A marker")
    if tokens & {"RSVB", "RSV-B"} or compact.startswith("RSVB"):
        return AnalysisTarget("RSV-B", None, "filename contains RSV-B marker")

    if tokens & {"SC2", "SARS", "SARS2", "COV2", "COVID", "COVID19"} or compact.startswith(("SC2", "SARS", "COVID")):
        return AnalysisTarget("SARS-CoV-2", None, "filename contains SARS-CoV-2 marker")

    for organism in sorted(available_organisms or [], key=len, reverse=True):
        if organism.upper().startswith("INFLUENZA"):
            continue
        organism_tokens = set(tokenize_name(organism))
        organism_compact = re.sub(r"[^A-Z0-9]+", "", organism.upper())
        if not organism_tokens or not organism_compact:
            continue
        if organism_tokens.issubset(tokens) or organism_compact in compact:
            return AnalysisTarget(organism, None, f"filename contains organism marker '{organism}'")

    return None


def normalize_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def find_metadata_column(fieldnames: list[str], candidates: set[str]) -> str | None:
    for fieldname in fieldnames:
        if normalize_column_name(fieldname) in candidates:
            return fieldname
    return None


def is_ct_metadata_column(fieldname: str, sample_col: str | None, date_col: str | None) -> bool:
    if fieldname in {sample_col, date_col}:
        return False
    normalized = normalize_column_name(fieldname)
    if normalized in {"ct", "ctvalue", "cq", "cp", "cyclethreshold"}:
        return True
    if normalized.startswith(("ct", "cq", "cp")) and len(normalized) > 2:
        return True
    if normalized.endswith(("ct", "cq", "cp")) and len(normalized) > 2:
        return True
    return False


def ct_column_aliases(fieldname: str) -> set[str]:
    normalized = normalize_column_name(fieldname)
    tokens = [token.lower() for token in re.split(r"[^A-Za-z0-9]+", fieldname) if token]
    aliases = {normalized}
    assay_tokens = [token for token in tokens if token not in {"ct", "cq", "cp", "value", "cycle", "threshold"}]
    if assay_tokens:
        aliases.add("".join(assay_tokens))
        aliases.update(assay_tokens)

    replacements = {
        "sc2": {"sarscov2", "sars2", "cov2", "covid19"},
        "sarscov2": {"sc2", "sars2", "cov2", "covid19"},
        "h3": {"h3n2"},
        "h1": {"h1n1"},
        "rsva": {"rsv-a"},
        "rsvb": {"rsv-b"},
        "infa": {"influenzaa", "flua"},
        "infb": {"influenzab", "flub"},
        "triplexinfa": {"infa", "influenzaa", "flua"},
        "triplexinfb": {"infb", "influenzab", "flub"},
    }
    for alias in list(aliases):
        aliases.update(replacements.get(alias, set()))
        if alias.startswith("ct") and len(alias) > 2:
            aliases.add(alias[2:])
        if alias.endswith("ct") and len(alias) > 2:
            aliases.add(alias[:-2])
    return {alias for alias in aliases if alias}


def primer_ct_aliases(virus_type: str, primer: PrimerRecord) -> list[str]:
    aliases = []
    virus = virus_type.upper()
    primer_name = primer.name.upper()
    primer_tokens = set(tokenize_name(primer_name))
    subtype_tags = {tag.upper() for tag in primer.subtype_tags}
    is_triplex_primer = "TRIPLEX" in primer_tokens
    is_influenza_a = primer.organism.upper() == "INFLUENZA-A" or virus in {"INFLUENZA-A", "INFLUENZA-H1", "INFLUENZA-H3"}
    is_influenza_b = primer.organism.upper() == "INFLUENZA-B" or virus == "INFLUENZA-B"

    if is_triplex_primer and is_influenza_a:
        aliases.extend(["triplexinfa", "infa", "influenzaa", "flua"])
    if is_triplex_primer and is_influenza_b:
        aliases.extend(["triplexinfb", "infb", "influenzab", "flub"])
    if "H3" in virus or "H3" in subtype_tags or "H3" in primer_tokens:
        aliases.extend(["h3", "h3n2"])
    if "H1" in virus or "H1" in subtype_tags or "H1" in primer_tokens:
        aliases.extend(["h1", "h1n1"])
    if is_influenza_a:
        aliases.extend(["infa", "influenzaa", "flua"])
    if is_influenza_b:
        aliases.extend(["infb", "influenzab", "flub"])
    if "SARS" in virus or "SC2" in primer_name or "COV" in virus:
        aliases.extend(["triplexsc2", "sc2", "sarscov2", "sars2", "cov2", "covid19"])
    if virus == "RSV-A":
        aliases.extend(["rsva", "rsv-a"])
    if virus == "RSV-B":
        aliases.extend(["rsvb", "rsv-b"])

    normalized_aliases = []
    for alias in aliases:
        normalized = normalize_column_name(alias)
        if normalized and normalized not in normalized_aliases:
            normalized_aliases.append(normalized)
    return normalized_aliases


def select_ct_value(record: MetadataRecord, virus_type: str, primer: PrimerRecord) -> tuple[str, str]:
    if not record:
        return "", ""
    alias_order = primer_ct_aliases(virus_type, primer)
    for alias in alias_order:
        for column_name, value in record.ct_values.items():
            if not value:
                continue
            if alias in ct_column_aliases(column_name):
                return value, column_name
    if record.ct_value:
        return record.ct_value, "Ct_Value"
    for column_name, value in record.ct_values.items():
        if value:
            return value, column_name
    return "", ""


def load_metadata_csv(metadata_file: str | None) -> tuple[list[MetadataRecord], ValidationResult]:
    """
    Load optional sample metadata.

    Accepted column names are intentionally flexible:
      sample id: SampleID, Sample_ID, Sample, ID
      date: Sample_Date, SampleDate, Date, Collection_Date
      Ct: Ct, CT, Ct_Value, Cq, Cp, CT_H3, CT_H1, Triplex-SC2_CT, etc.
    """
    validation = ValidationResult()
    if not metadata_file:
        return [], validation

    try:
        with open(metadata_file, "r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            if not fieldnames:
                validation.errors.append(f"Metadata CSV '{metadata_file}' has no header row.")
                return [], validation

            sample_col = find_metadata_column(fieldnames, {"sampleid", "sample", "id"})
            date_col = find_metadata_column(fieldnames, {"sampledate", "date", "collectiondate", "samplingdate"})
            ct_columns = [fieldname for fieldname in fieldnames if is_ct_metadata_column(fieldname, sample_col, date_col)]
            ct_col = find_metadata_column(ct_columns, {"ct", "ctvalue", "cq", "cp", "cyclethreshold"})
            if not sample_col:
                validation.errors.append(
                    f"Metadata CSV '{metadata_file}' needs a sample ID column such as SampleID or Sample_ID."
                )
                return [], validation
            if not date_col:
                validation.warnings.append(
                    f"Metadata CSV '{metadata_file}' has no recognized Sample_Date column; timeline dates will be unavailable."
                )
            if not ct_columns:
                validation.warnings.append(
                    f"Metadata CSV '{metadata_file}' has no recognized Ct assay columns; Ct timeline values will be unavailable."
                )

            records = []
            seen_sample_ids: dict[str, int] = {}
            for line_number, row in enumerate(reader, start=2):
                sample_id = (row.get(sample_col) or "").strip()
                if not sample_id:
                    validation.warnings.append(f"Metadata CSV '{metadata_file}' row {line_number} has an empty sample ID and was skipped.")
                    continue
                normalized_sample_id = sample_id.lower()
                seen_sample_ids[normalized_sample_id] = seen_sample_ids.get(normalized_sample_id, 0) + 1
                records.append(
                    MetadataRecord(
                        sample_id=sample_id,
                        sample_date=(row.get(date_col) or "").strip() if date_col else "",
                        ct_value=(row.get(ct_col) or "").strip() if ct_col else "",
                        ct_values={
                            column: (row.get(column) or "").strip()
                            for column in ct_columns
                        },
                    )
                )

            for sample_id, count in seen_sample_ids.items():
                if count > 1:
                    validation.warnings.append(
                        f"Metadata CSV '{metadata_file}' contains {count} rows for sample ID '{sample_id}'. "
                        "Ambiguous FASTA matches will not be attached."
                    )
            return records, validation
    except Exception as exc:
        validation.errors.append(f"Error loading metadata CSV '{metadata_file}': {exc}")
        return [], validation


def clean_fasta_subject_id(subject_id: str) -> str:
    return subject_id.strip().lstrip(">").strip()


def split_match_tokens(value: str) -> list[str]:
    return [token for token in re.split(r"[|\\/\s;:,]+", value.strip()) if token]


def metadata_sample_matches_subject(sample_id: str, subject_id: str) -> bool:
    """
    Match metadata IDs to FASTA headers without unsafe substring matching.

    Examples:
      sample ID 454511 matches FASTA header Genome|454511
      sample ID 4545 does not match FASTA header Genome|454511
      sample ID Genome|4545 matches FASTA header Genome|4545, but not Genome|4546
    """
    sample = sample_id.strip()
    subject = clean_fasta_subject_id(subject_id)
    if not sample or not subject:
        return False
    if sample.casefold() == subject.casefold():
        return True
    sample_tokens = split_match_tokens(sample)
    subject_tokens = split_match_tokens(subject)
    if len(sample_tokens) == 1:
        sample_token = sample_tokens[0].casefold()
        return any(sample_token == token.casefold() for token in subject_tokens)
    escaped = re.escape(sample)
    return bool(re.search(rf"(?<![A-Za-z0-9]){escaped}(?![A-Za-z0-9])", subject, flags=re.IGNORECASE))


def choose_metadata_match(subject_id: str, metadata_records: list[MetadataRecord]) -> tuple[MetadataRecord | None, str]:
    matches = [record for record in metadata_records if metadata_sample_matches_subject(record.sample_id, subject_id)]
    if not matches:
        return None, "no_match"

    subject = clean_fasta_subject_id(subject_id).casefold()
    exact_matches = [record for record in matches if record.sample_id.casefold() == subject]
    if len(exact_matches) == 1:
        return exact_matches[0], "matched"
    if len(matches) == 1:
        return matches[0], "matched"
    return None, "ambiguous"


def build_metadata_matches(
    subject_ids: list[str],
    metadata_records: list[MetadataRecord],
    fasta_file: str,
    quiet: bool = False,
) -> dict[str, MetadataRecord]:
    if not metadata_records:
        return {}

    matches = {}
    ambiguous = []
    for subject_id in subject_ids:
        match, status = choose_metadata_match(subject_id, metadata_records)
        if match:
            matches[subject_id] = match
        elif status == "ambiguous":
            ambiguous.append(subject_id)

    if ambiguous and not quiet:
        preview = ", ".join(ambiguous[:5])
        suffix = "..." if len(ambiguous) > 5 else ""
        print(
            f"Warning: {len(ambiguous)} subject header(s) in {fasta_file} matched multiple metadata rows "
            f"and were left without metadata: {preview}{suffix}",
            file=sys.stderr,
        )
    return matches


def metadata_columns_for_subject(
    subject_id: str,
    metadata_matches: dict[str, MetadataRecord],
    virus_type: str,
    primer: PrimerRecord,
) -> dict[str, str]:
    record = metadata_matches.get(subject_id)
    if not record:
        return {
            "Metadata_Sample_ID": "",
            "Sample_Date": "",
            "Ct_Value": "",
            "Ct_Source": "",
        }
    ct_value, ct_source = select_ct_value(record, virus_type, primer)
    return {
        "Metadata_Sample_ID": record.sample_id,
        "Sample_Date": record.sample_date,
        "Ct_Value": ct_value,
        "Ct_Source": ct_source,
    }


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "unknown"


def infer_pool(primer_name: str) -> str:
    tokens = tokenize_name(primer_name)
    return tokens[0].lower() if tokens and tokens[0] in {"TRIPLEX"} else ""


def infer_normalized_role(primer_name: str) -> str:
    role = infer_primer_role(primer_name)
    if role == "probe":
        return "probe"
    tokens = tokenize_name(primer_name)
    if "LEFT" in tokens:
        return "forward_primer"
    if "RIGHT" in tokens:
        return "reverse_primer"
    if any(token.startswith("F") for token in tokens):
        return "forward_primer"
    if any(token.startswith("R") for token in tokens):
        return "reverse_primer"
    return "primer"


def primer_records_to_legacy_dict(records: list[PrimerRecord]) -> dict[str, str]:
    return {record.name: record.sequence for record in records}


def merge_validation_results(*results: ValidationResult) -> ValidationResult:
    merged = ValidationResult()
    for result in results:
        merged.errors.extend(result.errors)
        merged.warnings.extend(result.warnings)
    return merged


def merge_primer_record_maps(*record_maps: dict[str, list[PrimerRecord]]) -> dict[str, list[PrimerRecord]]:
    merged: dict[str, list[PrimerRecord]] = {}
    for record_map in record_maps:
        for organism, records in record_map.items():
            merged.setdefault(organism, []).extend(records)
    return merged


def validate_legacy_primer_library(primer_library: object) -> ValidationResult:
    """
    Validate the legacy {organism: {primer_name: sequence}} primer database format.
    This intentionally does not require metadata that the legacy format cannot store.
    """
    result = ValidationResult()
    if not isinstance(primer_library, dict):
        result.errors.append("Primer library must be a JSON object mapping organism names to primer objects.")
        return result
    if not primer_library:
        result.errors.append("Primer library is empty.")
        return result

    for organism, primers in primer_library.items():
        if not isinstance(organism, str) or not organism.strip():
            result.errors.append("Primer library contains an empty or non-string organism name.")
            continue
        if not isinstance(primers, dict):
            result.errors.append(f"Primer library section '{organism}' must be an object mapping primer names to sequences.")
            continue
        if not primers:
            result.errors.append(f"Primer library section '{organism}' contains no primers.")
            continue

        sequences_seen: dict[str, str] = {}
        for primer_name, sequence in primers.items():
            context = f"{organism}/{primer_name}"
            if not isinstance(primer_name, str) or not primer_name.strip():
                result.errors.append(f"{organism} contains an empty or non-string primer name.")
                continue
            if not isinstance(sequence, str):
                result.errors.append(f"{context} sequence must be a string.")
                continue

            normalized_sequence = sequence.strip().upper()
            if not normalized_sequence:
                result.errors.append(f"{context} sequence is empty.")
                continue

            invalid_codes = sorted(set(normalized_sequence) - ALLOWED_SEQUENCE_CODES)
            if invalid_codes:
                result.errors.append(f"{context} sequence contains unsupported IUPAC code(s): {', '.join(invalid_codes)}.")

            previous_name = sequences_seen.get(normalized_sequence)
            if previous_name:
                result.warnings.append(f"{context} has the same sequence as {organism}/{previous_name}.")
            else:
                sequences_seen[normalized_sequence] = primer_name

            if organism.upper().startswith("INFLUENZA") and not infer_primer_segment(organism, primer_name):
                result.warnings.append(f"{context} has no inferable HA, M, or NS segment in the legacy primer name.")

            if infer_primer_role(primer_name) == "probe":
                result.warnings.append(f"{context} appears to be a probe; the legacy format cannot store role metadata explicitly.")

    return result


def legacy_library_to_records(primer_library: dict) -> dict[str, list[PrimerRecord]]:
    records_by_organism: dict[str, list[PrimerRecord]] = {}
    for organism, primers in primer_library.items():
        records = []
        for primer_name, sequence in primers.items():
            records.append(
                PrimerRecord(
                    organism=organism,
                    name=primer_name.strip(),
                    sequence=sequence.strip().upper(),
                    segment=infer_primer_segment(organism, primer_name),
                    role=infer_primer_role(primer_name),
                    subtype_tags=infer_subtype_tags(primer_name),
                )
            )
        records_by_organism[organism] = records
    return records_by_organism


def is_normalized_primer_library(primer_library: object) -> bool:
    return isinstance(primer_library, dict) and "schema_version" in primer_library and "schemes" in primer_library


def is_panel_primer_library(primer_library: object) -> bool:
    return isinstance(primer_library, dict) and "schema_version" in primer_library and "panels" in primer_library


def is_virus_organized_primer_library(primer_library: object) -> bool:
    return isinstance(primer_library, dict) and "schema_version" in primer_library and "viruses" in primer_library


def _section_list(section: object, key: str) -> list:
    if section is None:
        return []
    if isinstance(section, list):
        return section
    if isinstance(section, dict):
        value = section.get(key, [])
        return value if isinstance(value, list) else []
    return []


def flatten_virus_organized_library(primer_library: dict) -> dict:
    """Convert viruses[].pcr.schemes[] / viruses[].ngs.panels[] into flat schemes[] / panels[]."""
    flat = {
        "schema_version": primer_library.get("schema_version", ""),
        "database_version": primer_library.get("database_version", ""),
        "schemes": [],
        "panels": [],
    }
    scheme_ids_seen = set()
    panel_ids_seen = set()

    for virus in primer_library.get("viruses", []):
        if not isinstance(virus, dict):
            continue
        organism = virus.get("organism", "")
        for scheme in _section_list(virus.get("pcr"), "schemes"):
            if not isinstance(scheme, dict):
                continue
            scheme = dict(scheme)
            scheme.setdefault("organism", organism)
            scheme_id = scheme.get("scheme_id")
            if scheme_id and scheme_id not in scheme_ids_seen:
                flat["schemes"].append(scheme)
                scheme_ids_seen.add(scheme_id)
        for panel in _section_list(virus.get("ngs"), "panels"):
            if not isinstance(panel, dict):
                continue
            panel = dict(panel)
            panel.setdefault("organism", organism)
            panel_id = panel.get("panel_id")
            if panel_id and panel_id not in panel_ids_seen:
                flat["panels"].append(panel)
                panel_ids_seen.add(panel_id)
    return flat


def validate_virus_organized_primer_library(primer_library: object, library_file: str | None = None) -> ValidationResult:
    result = ValidationResult()
    if not isinstance(primer_library, dict):
        result.errors.append("Virus-organized primer library must be a JSON object.")
        return result
    for field_name in ("schema_version", "database_version", "viruses"):
        if field_name not in primer_library:
            result.errors.append(f"Virus-organized primer library is missing required field '{field_name}'.")
    if result.errors:
        return result
    if not isinstance(primer_library["schema_version"], str) or not primer_library["schema_version"].strip():
        result.errors.append("Virus-organized primer library field 'schema_version' must be a non-empty string.")
    if not isinstance(primer_library["database_version"], str) or not primer_library["database_version"].strip():
        result.errors.append("Virus-organized primer library field 'database_version' must be a non-empty string.")
    viruses = primer_library.get("viruses")
    if not isinstance(viruses, list) or not viruses:
        result.errors.append("Virus-organized primer library field 'viruses' must be a non-empty list.")
        return result

    organism_names = set()
    scheme_ids_seen = set()
    panel_ids_seen = set()
    for virus_index, virus in enumerate(viruses, start=1):
        context = f"virus #{virus_index}"
        if not isinstance(virus, dict):
            result.errors.append(f"{context} must be an object.")
            continue
        organism = virus.get("organism", "")
        if not isinstance(organism, str) or not organism.strip():
            result.errors.append(f"{context} field 'organism' must be a non-empty string.")
            continue
        if organism in organism_names:
            result.errors.append(f"{context} has duplicate organism '{organism}'.")
        organism_names.add(organism)
        if "pcr" not in virus and "ngs" not in virus:
            result.errors.append(f"{context} must contain at least one of 'pcr' or 'ngs'.")

        for scheme in _section_list(virus.get("pcr"), "schemes"):
            if not isinstance(scheme, dict):
                continue
            scheme_id = scheme.get("scheme_id", "")
            if scheme_id in scheme_ids_seen:
                result.errors.append(f"{context} has duplicate PCR scheme_id '{scheme_id}'.")
            elif scheme_id:
                scheme_ids_seen.add(scheme_id)
            if scheme.get("organism") and scheme.get("organism") != organism:
                result.errors.append(f"{context} PCR scheme '{scheme_id}' organism does not match '{organism}'.")

        for panel in _section_list(virus.get("ngs"), "panels"):
            if not isinstance(panel, dict):
                continue
            panel_id = panel.get("panel_id", "")
            if panel_id in panel_ids_seen:
                result.errors.append(f"{context} has duplicate NGS panel_id '{panel_id}'.")
            elif panel_id:
                panel_ids_seen.add(panel_id)
            if panel.get("organism") and panel.get("organism") != organism:
                result.errors.append(f"{context} NGS panel '{panel_id}' organism does not match '{organism}'.")

    if result.errors:
        return result

    flat = flatten_virus_organized_library(primer_library)
    flat_validations = []
    if flat.get("schemes"):
        flat_validations.append(validate_normalized_primer_library(flat))
    if flat.get("panels"):
        flat_validations.append(validate_panel_primer_library(flat, library_file=library_file))
    if not flat_validations:
        result.errors.append("Virus-organized primer library does not contain any PCR schemes or NGS panels.")
        return result
    return merge_validation_results(result, *flat_validations)


def resolve_panel_asset_path(library_file: str, asset_path: str) -> str:
    """Resolve panel asset paths relative to the primer database JSON file."""
    if os.path.isabs(asset_path):
        return asset_path
    return os.path.join(os.path.dirname(os.path.abspath(library_file)), asset_path)


def bed_field_index(field_spec: object, default: int) -> int:
    """Return a zero-based BED field index from a schema mapping value."""
    aliases = {
        "chrom": 0,
        "reference": 0,
        "start": 1,
        "end": 2,
        "name": 3,
        "amplicon_id": 3,
        "primer_id": 3,
        "score": 4,
        "pool": 4,
        "strand": 5,
        "sequence": 6,
        "primer_sequence": 6,
    }
    if field_spec is None:
        return default
    if isinstance(field_spec, int):
        return field_spec
    if isinstance(field_spec, str):
        stripped = field_spec.strip()
        if stripped.isdigit():
            return int(stripped)
        return aliases.get(stripped.lower(), default)
    return default


def read_panel_fasta_records(fasta_file: str) -> dict[str, str]:
    return read_fasta_sequences(fasta_file)


def read_panel_bed_records(bed_file: str, mapping: dict | None = None) -> list[dict[str, str]]:
    """
    Read BED-like panel rows.

    ARTIC-style primer BED files commonly use columns:
      chrom, start, end, name, pool, strand, sequence
    The mapping object can override name/pool/sequence indexes with zero-based
    field numbers or aliases such as "amplicon_id" and "sequence".
    """
    mapping = mapping or {}
    name_index = bed_field_index(mapping.get("bed_name_field"), 3)
    pool_index = bed_field_index(mapping.get("pool_from_bed_field"), 4)
    sequence_index = bed_field_index(mapping.get("bed_sequence_field"), 6)
    entries: list[dict[str, str]] = []
    with open(bed_file, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.lower().startswith(("track", "browser")):
                continue
            fields = stripped.split()
            if len(fields) < 3:
                continue
            name = fields[name_index].strip() if name_index < len(fields) else f"{fields[0]}:{fields[1]}-{fields[2]}"
            entries.append(
                {
                    "chrom": fields[0],
                    "start": fields[1],
                    "end": fields[2],
                    "name": name,
                    "pool": fields[pool_index].strip() if pool_index < len(fields) else "",
                    "strand": fields[5].strip() if len(fields) > 5 else "",
                    "sequence": fields[sequence_index].strip().upper() if sequence_index < len(fields) else "",
                    "line_number": str(line_number),
                }
            )
    return entries


def primer_name_candidates_from_pattern(pattern: str, amplicon_id: str) -> dict[str, str]:
    """Build possible primer names from a {amplicon_id}_{role} pattern."""
    if not pattern or "{amplicon_id}" not in pattern or "{role}" not in pattern:
        return {}
    roles = {
        "LEFT": "forward_primer",
        "RIGHT": "reverse_primer",
        "F": "forward_primer",
        "R": "reverse_primer",
        "FORWARD": "forward_primer",
        "REVERSE": "reverse_primer",
    }
    return {pattern.format(amplicon_id=amplicon_id, role=role): normalized_role for role, normalized_role in roles.items()}


def bed_metadata_by_primer_name(entries: list[dict[str, str]], mapping: dict | None = None) -> dict[str, dict[str, str]]:
    mapping = mapping or {}
    metadata: dict[str, dict[str, str]] = {}
    pattern = mapping.get("primer_name_pattern", "")
    for entry in entries:
        metadata[entry["name"]] = entry
        for candidate_name in primer_name_candidates_from_pattern(pattern, entry["name"]):
            metadata.setdefault(candidate_name, entry)
    return metadata


def panel_primer_role(primer_name: str, bed_entry: dict[str, str] | None = None, pattern_roles: dict[str, str] | None = None) -> str:
    if pattern_roles and primer_name in pattern_roles:
        return pattern_roles[primer_name]
    role = infer_normalized_role(primer_name)
    if role != "primer":
        return role
    if bed_entry and bed_entry.get("strand") == "-":
        return "reverse_primer"
    if bed_entry and bed_entry.get("strand") == "+":
        return "forward_primer"
    return role


def validate_panel_primer_library(primer_library: object, library_file: str | None = None) -> ValidationResult:
    """Validate schema v2 panel databases that point to BED/FASTA assets."""
    result = ValidationResult()
    if not isinstance(primer_library, dict):
        result.errors.append("Panel primer library must be a JSON object.")
        return result

    for field_name in ("schema_version", "database_version", "panels"):
        if field_name not in primer_library:
            result.errors.append(f"Panel primer library is missing required field '{field_name}'.")
    if result.errors:
        return result

    if not isinstance(primer_library["schema_version"], str) or not primer_library["schema_version"].strip():
        result.errors.append("Panel primer library field 'schema_version' must be a non-empty string.")
    if not isinstance(primer_library["database_version"], str) or not primer_library["database_version"].strip():
        result.errors.append("Panel primer library field 'database_version' must be a non-empty string.")

    panels = primer_library["panels"]
    if not isinstance(panels, list) or not panels:
        result.errors.append("Panel primer library field 'panels' must be a non-empty list.")
        return result

    panel_ids_seen = set()
    for panel_index, panel in enumerate(panels, start=1):
        context = f"panel #{panel_index}"
        if not isinstance(panel, dict):
            result.errors.append(f"{context} must be an object.")
            continue
        for field_name in ("panel_id", "display_name", "organism", "panel_version", "files"):
            if field_name not in panel:
                result.errors.append(f"{context} is missing required field '{field_name}'.")
        panel_id = panel.get("panel_id", "")
        if isinstance(panel_id, str) and panel_id.strip():
            if panel_id in panel_ids_seen:
                result.errors.append(f"{context} has duplicate panel_id '{panel_id}'.")
            panel_ids_seen.add(panel_id)
        else:
            result.errors.append(f"{context} field 'panel_id' must be a non-empty string.")
        for field_name in ("display_name", "organism", "panel_version"):
            if field_name in panel and (not isinstance(panel[field_name], str) or not panel[field_name].strip()):
                result.errors.append(f"{context} field '{field_name}' must be a non-empty string.")
        files = panel.get("files")
        if not isinstance(files, dict):
            result.errors.append(f"{context} field 'files' must be an object.")
            continue
        if not files.get("bed") and not files.get("primers_fasta"):
            result.errors.append(f"{context} files must include at least 'bed' or 'primers_fasta'.")
        if library_file:
            for file_key in ("bed", "primers_fasta"):
                asset_path = files.get(file_key)
                if asset_path:
                    resolved = resolve_panel_asset_path(library_file, asset_path)
                    if not os.path.exists(resolved):
                        result.errors.append(f"{context} file '{file_key}' does not exist: {resolved}")
    return result


def panel_library_to_records(primer_library: dict, library_file: str, validation: ValidationResult | None = None) -> dict[str, list[PrimerRecord]]:
    records_by_organism: dict[str, list[PrimerRecord]] = {}
    database_version = primer_library.get("database_version", "")
    validation = validation or ValidationResult()

    for panel in primer_library.get("panels", []):
        organism = panel["organism"]
        panel_id = panel["panel_id"]
        panel_version = panel["panel_version"]
        panel_name = panel["display_name"]
        files = panel.get("files", {})
        mapping = panel.get("mapping") or {}
        reference = panel.get("reference") or {}
        records = records_by_organism.setdefault(organism, [])
        by_name: dict[str, PrimerRecord] = {}

        bed_entries: list[dict[str, str]] = []
        bed_path = ""
        if files.get("bed"):
            bed_path = resolve_panel_asset_path(library_file, files["bed"])
            bed_entries = read_panel_bed_records(bed_path, mapping=mapping)
            pattern_roles = {
                candidate_name: role
                for entry in bed_entries
                for candidate_name, role in primer_name_candidates_from_pattern(mapping.get("primer_name_pattern", ""), entry["name"]).items()
            }
            for entry in bed_entries:
                sequence = entry.get("sequence", "")
                if not sequence:
                    continue
                name = entry["name"]
                by_name[name] = PrimerRecord(
                    organism=organism,
                    name=name,
                    sequence=sequence,
                    segment=(panel.get("segment") or "").strip().upper(),
                    role=panel_primer_role(name, entry, pattern_roles),
                    subtype_tags=tuple(panel.get("subtype_tags") or infer_subtype_tags(name)),
                    scheme_id=panel_id,
                    scheme_version=panel_version,
                    database_version=database_version,
                    pool=entry.get("pool", ""),
                    strand=entry.get("strand", ""),
                    reference_name=reference.get("name", entry.get("chrom", "")),
                    start=entry.get("start", ""),
                    end=entry.get("end", ""),
                    source_file=files.get("bed", ""),
                    assay_type="ngs",
                    assay_name=panel_name,
                )

        bed_metadata = bed_metadata_by_primer_name(bed_entries, mapping=mapping)
        pattern_roles = {
            candidate_name: role
            for entry in bed_entries
            for candidate_name, role in primer_name_candidates_from_pattern(mapping.get("primer_name_pattern", ""), entry["name"]).items()
        }

        if files.get("primers_fasta"):
            fasta_path = resolve_panel_asset_path(library_file, files["primers_fasta"])
            for name, sequence in read_panel_fasta_records(fasta_path).items():
                entry = bed_metadata.get(name, {})
                by_name[name] = PrimerRecord(
                    organism=organism,
                    name=name,
                    sequence=sequence,
                    segment=(panel.get("segment") or "").strip().upper(),
                    role=panel_primer_role(name, entry, pattern_roles),
                    subtype_tags=tuple(panel.get("subtype_tags") or infer_subtype_tags(name)),
                    scheme_id=panel_id,
                    scheme_version=panel_version,
                    database_version=database_version,
                    pool=entry.get("pool", ""),
                    strand=entry.get("strand", ""),
                    reference_name=reference.get("name", entry.get("chrom", "")),
                    start=entry.get("start", ""),
                    end=entry.get("end", ""),
                    source_file=files.get("primers_fasta", ""),
                    assay_type="ngs",
                    assay_name=panel_name,
                )

        for record in by_name.values():
            invalid_codes = sorted(set(record.sequence) - ALLOWED_SEQUENCE_CODES)
            if invalid_codes:
                validation.errors.append(
                    f"Panel '{panel_id}' primer '{record.name}' sequence contains unsupported IUPAC code(s): "
                    f"{', '.join(invalid_codes)}."
                )

        if not by_name:
            source = bed_path or files.get("primers_fasta", "")
            validation.errors.append(
                f"Panel '{panel_id}' did not yield any primer sequences from {source}. "
                "Use an ARTIC-style primer BED with sequence in column 7 or provide files.primers_fasta."
            )
        records.extend(by_name[name] for name in sorted(by_name))

    return records_by_organism


def validate_normalized_primer_library(primer_library: object) -> ValidationResult:
    """Validate the normalized schema_version/database_version/schemes primer database format."""
    result = ValidationResult()
    if not isinstance(primer_library, dict):
        result.errors.append("Normalized primer library must be a JSON object.")
        return result

    for field_name in ("schema_version", "database_version", "schemes"):
        if field_name not in primer_library:
            result.errors.append(f"Normalized primer library is missing required field '{field_name}'.")

    if result.errors:
        return result

    if not isinstance(primer_library["schema_version"], str) or not primer_library["schema_version"].strip():
        result.errors.append("Normalized primer library field 'schema_version' must be a non-empty string.")
    if not isinstance(primer_library["database_version"], str) or not primer_library["database_version"].strip():
        result.errors.append("Normalized primer library field 'database_version' must be a non-empty string.")

    schemes = primer_library["schemes"]
    if not isinstance(schemes, list) or not schemes:
        result.errors.append("Normalized primer library field 'schemes' must be a non-empty list.")
        return result

    scheme_ids_seen = set()
    for scheme_index, scheme in enumerate(schemes, start=1):
        context = f"scheme #{scheme_index}"
        if not isinstance(scheme, dict):
            result.errors.append(f"{context} must be an object.")
            continue

        for field_name in ("scheme_id", "display_name", "organism", "version", "primers"):
            if field_name not in scheme:
                result.errors.append(f"{context} is missing required field '{field_name}'.")

        scheme_id = scheme.get("scheme_id", "")
        organism = scheme.get("organism", "")
        if isinstance(scheme_id, str) and scheme_id.strip():
            if scheme_id in scheme_ids_seen:
                result.errors.append(f"{context} has duplicate scheme_id '{scheme_id}'.")
            scheme_ids_seen.add(scheme_id)
        else:
            result.errors.append(f"{context} field 'scheme_id' must be a non-empty string.")

        for field_name in ("display_name", "organism", "version"):
            if field_name in scheme and (not isinstance(scheme[field_name], str) or not scheme[field_name].strip()):
                result.errors.append(f"{context} field '{field_name}' must be a non-empty string.")

        primers = scheme.get("primers")
        if not isinstance(primers, list) or not primers:
            result.errors.append(f"{context} field 'primers' must be a non-empty list.")
            continue

        primer_ids_seen = set()
        sequences_seen: dict[str, str] = {}
        for primer_index, primer in enumerate(primers, start=1):
            primer_context = f"{scheme_id or context}/primer #{primer_index}"
            if not isinstance(primer, dict):
                result.errors.append(f"{primer_context} must be an object.")
                continue

            for field_name in ("id", "name", "sequence", "role", "segment"):
                if field_name not in primer:
                    result.errors.append(f"{primer_context} is missing required field '{field_name}'.")

            primer_id = primer.get("id", "")
            primer_name = primer.get("name", primer_id)
            sequence = primer.get("sequence", "")
            if isinstance(primer_id, str) and primer_id.strip():
                if primer_id in primer_ids_seen:
                    result.errors.append(f"{primer_context} has duplicate primer id '{primer_id}'.")
                primer_ids_seen.add(primer_id)
            else:
                result.errors.append(f"{primer_context} field 'id' must be a non-empty string.")

            for field_name in ("name", "role"):
                if field_name in primer and (not isinstance(primer[field_name], str) or not primer[field_name].strip()):
                    result.errors.append(f"{primer_context} field '{field_name}' must be a non-empty string.")

            if not isinstance(sequence, str) or not sequence.strip():
                result.errors.append(f"{primer_context} field 'sequence' must be a non-empty string.")
                continue

            normalized_sequence = sequence.strip().upper()
            invalid_codes = sorted(set(normalized_sequence) - ALLOWED_SEQUENCE_CODES)
            if invalid_codes:
                result.errors.append(f"{primer_context} sequence contains unsupported IUPAC code(s): {', '.join(invalid_codes)}.")

            previous_name = sequences_seen.get(normalized_sequence)
            if previous_name:
                result.warnings.append(f"{primer_context} has the same sequence as {previous_name}.")
            else:
                sequences_seen[normalized_sequence] = str(primer_name)

            segment = primer.get("segment", "")
            if segment is None:
                segment = ""
            if not isinstance(segment, str):
                result.errors.append(f"{primer_context} field 'segment' must be a string.")
            elif organism.upper().startswith("INFLUENZA") and segment.upper() not in INFLUENZA_SEGMENTS:
                result.errors.append(f"{primer_context} has unsupported influenza segment '{segment}'.")

    return result


def normalized_library_to_records(primer_library: dict) -> dict[str, list[PrimerRecord]]:
    records_by_organism: dict[str, list[PrimerRecord]] = {}
    database_version = primer_library.get("database_version", "")
    for scheme in primer_library.get("schemes", []):
        organism = scheme["organism"]
        scheme_id = scheme["scheme_id"]
        scheme_version = scheme["version"]
        records = records_by_organism.setdefault(organism, [])
        for primer in scheme["primers"]:
            subtype_tags = tuple(primer.get("subtype_tags") or infer_subtype_tags(primer.get("name", "")))
            records.append(
                PrimerRecord(
                    organism=organism,
                    name=primer["name"].strip(),
                    sequence=primer["sequence"].strip().upper(),
                    segment=(primer.get("segment") or "").strip().upper(),
                    role=primer.get("role", "primer").strip(),
                    subtype_tags=subtype_tags,
                    scheme_id=scheme_id,
                    scheme_version=scheme_version,
                    database_version=database_version,
                    pool=(primer.get("pool") or "").strip(),
                    strand=(primer.get("strand") or "").strip(),
                    assay_type="pcr",
                    assay_name=scheme["display_name"],
                )
            )
    return records_by_organism


def legacy_library_to_normalized_database(
    primer_library: dict,
    database_version: str,
    source: str = "legacy",
) -> dict:
    """Convert legacy {organism: {primer_name: sequence}} data to normalized schema v1."""
    normalized = {
        "schema_version": "1.0",
        "database_version": database_version,
        "schemes": [],
    }
    for organism, primers in primer_library.items():
        scheme_id = f"{slugify(source)}-{slugify(organism)}"
        scheme = {
            "scheme_id": scheme_id,
            "display_name": f"{source} {organism} primers",
            "organism": organism,
            "version": database_version,
            "status": "current",
            "source": source,
            "references": [],
            "primers": [],
        }
        for primer_name, sequence in primers.items():
            segment = infer_primer_segment(organism, primer_name)
            subtype_tags = list(infer_subtype_tags(primer_name))
            scheme["primers"].append(
                {
                    "id": primer_name,
                    "name": primer_name,
                    "sequence": sequence.strip().upper(),
                    "role": infer_normalized_role(primer_name),
                    "segment": segment,
                    "gene": segment,
                    "pool": infer_pool(primer_name),
                    "strand": "",
                    "subtype_tags": subtype_tags,
                    "notes": "Converted from legacy primer JSON; metadata inferred from primer name where possible.",
                }
            )
        normalized["schemes"].append(scheme)
    return normalized


def load_primer_records(library_file: str) -> tuple[dict[str, list[PrimerRecord]], ValidationResult]:
    primer_library = load_primer_library(library_file)
    was_virus_organized = is_virus_organized_primer_library(primer_library)
    if was_virus_organized:
        validation = validate_virus_organized_primer_library(primer_library, library_file=library_file)
        if validation.errors:
            error_text = "\n".join(f"- {error}" for error in validation.errors)
            sys.exit(f"Primer library validation failed for {library_file}:\n{error_text}")
        primer_library = flatten_virus_organized_library(primer_library)
    else:
        validation = ValidationResult()

    has_schemes = is_normalized_primer_library(primer_library)
    has_panels = is_panel_primer_library(primer_library)

    if has_schemes or has_panels:
        record_maps = []

        if not was_virus_organized:
            if has_schemes:
                validation = merge_validation_results(validation, validate_normalized_primer_library(primer_library))
            if has_panels:
                validation = merge_validation_results(validation, validate_panel_primer_library(primer_library, library_file=library_file))
        if validation.errors:
            error_text = "\n".join(f"- {error}" for error in validation.errors)
            sys.exit(f"Primer library validation failed for {library_file}:\n{error_text}")

        if has_schemes:
            record_maps.append(normalized_library_to_records(primer_library))
        if has_panels:
            record_maps.append(panel_library_to_records(primer_library, library_file, validation=validation))
        if validation.errors:
            error_text = "\n".join(f"- {error}" for error in validation.errors)
            sys.exit(f"Primer library validation failed for {library_file}:\n{error_text}")
        return merge_primer_record_maps(*record_maps), validation

    validation = validate_legacy_primer_library(primer_library)
    if validation.errors:
        error_text = "\n".join(f"- {error}" for error in validation.errors)
        sys.exit(f"Primer library validation failed for {library_file}:\n{error_text}")
    return legacy_library_to_records(primer_library), validation


def print_validation_messages(validation: ValidationResult):
    for warning in validation.warnings:
        print(f"Primer library warning: {warning}", file=sys.stderr)


def print_metadata_validation_messages(validation: ValidationResult):
    for warning in validation.warnings:
        print(f"Metadata warning: {warning}", file=sys.stderr)


def load_primer_library(library_file: str) -> dict:
    """
    Loads the primer library from a JSON file.
    The file should contain a dictionary mapping virus types to their primer dictionaries.
    """
    try:
        with open(library_file, "r") as f:
            return json.load(f)
    except Exception as e:
        sys.exit(f"Error loading primer library from {library_file}: {e}")

def build_influenza_subtype_primers(primer_library: dict, subtype: str) -> dict:
    """
    Legacy dict helper for selecting Influenza-A H1/H3 primers.

    H1/H3 use primer names tagged for that subtype plus untagged Influenza-A
    primer names, and exclude names tagged for the other subtype.
    """
    generic = primer_library.get("Influenza-A")
    if generic is None:
        sys.exit("Primer JSON lacks the 'Influenza-A' section.")

    subtype = subtype.upper()
    if subtype == "A":
        return generic  # full A panel

    if subtype in {"H1", "H3"}:
        subset = {}
        for primer_name, sequence in generic.items():
            tags = infer_subtype_tags(primer_name)
            if subtype in tags or not tags:
                subset[primer_name] = sequence
        if not subset:
            sys.exit(f"No {subtype} primers found inside 'Influenza-A'.")
        return subset

    sys.exit(f"Unsupported Influenza-A subtype '{subtype}'.")


def build_influenza_subtype_records(primer_records: dict[str, list[PrimerRecord]], subtype: str) -> list[PrimerRecord]:
    """
    Select Influenza-A records for A, H1, or H3. H1/H3 include subtype-specific
    records plus records without subtype tags, and exclude records tagged for the
    other subtype.
    """
    generic = primer_records.get("Influenza-A")
    if generic is None:
        sys.exit("Primer JSON lacks the 'Influenza-A' section.")

    subtype = subtype.upper()
    if subtype == "A":
        return generic

    if subtype in {"H1", "H3"}:
        subset = [
            record for record in generic
            if subtype in record.subtype_tags or not record.subtype_tags
        ]
        if not subset:
            sys.exit(f"No {subtype} primers found inside 'Influenza-A'.")
        return subset

    sys.exit(f"Unsupported Influenza-A subtype '{subtype}'.")


def _filter_assay_records(
    records: list[PrimerRecord],
    assay_type: str = "all",
    assay_id: str | None = None,
) -> list[PrimerRecord]:
    """Filter records by PCR/NGS source and optional scheme or panel ID."""
    selected = records
    if assay_type != "all":
        selected = [record for record in selected if record.assay_type == assay_type]
    if assay_id:
        selected = [record for record in selected if record.scheme_id == assay_id]
    if selected:
        return selected

    available = sorted({record.scheme_id for record in records if record.scheme_id})
    criteria = f"assay type '{assay_type}'"
    if assay_id:
        criteria += f" and assay ID '{assay_id}'"
    suffix = f" Available assay IDs: {', '.join(available)}." if available else ""
    sys.exit(f"No primers match {criteria}.{suffix}")


def select_primer_records(
    primer_records: dict[str, list[PrimerRecord]],
    virus_type: str,
    flu_type: str | None,
    assay_type: str = "all",
    assay_id: str | None = None,
) -> tuple[str, list[PrimerRecord]]:
    organism_by_casefold = {organism.casefold(): organism for organism in primer_records}
    if virus_type.lower() != "influenza":
        canonical_virus = organism_by_casefold.get(virus_type.casefold(), virus_type)
        selected = primer_records.get(canonical_virus)
        if not selected:
            sys.exit(f"No primer information available for virus type '{virus_type}'.")
        return canonical_virus, _filter_assay_records(selected, assay_type, assay_id)

    if not flu_type:
        sys.exit("Error: for influenza please supply --flu-type A, H1, H3 or B.")

    subtype = flu_type.upper()
    if subtype == "B":
        selected = primer_records.get("Influenza-B")
        if not selected:
            sys.exit("Primer JSON lacks the 'Influenza-B' section.")
        return "Influenza-B", _filter_assay_records(selected, assay_type, assay_id)

    selected = build_influenza_subtype_records(primer_records, subtype)
    return f"Influenza-{subtype}", _filter_assay_records(selected, assay_type, assay_id)

def get_segment(subject_id: str) -> str:
    """
    Given a subject FASTA header (without the leading '>'),
    extract the segment code from supported tokens such as:
      "03-M|252500127", "03-HA|252500127", or "contig1|03-M|INFL16-2025"
    Returns the segment code (e.g., "M" or "HA") or an empty string if not found.
    """
    for token in subject_id.split("|"):
        match = re.fullmatch(r"\d{1,2}-([A-Za-z0-9]+)", token.strip())
        if match:
            return match.group(1).upper()
    return ""

def get_subject_ids(fasta_file: str) -> list:
    """
    Extract subject sequence IDs from a FASTA file.
    Reads each header (lines starting with '>') and takes the first word as the ID.
    """
    subject_ids = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                subject_id = line[1:].strip().split()[0]
                subject_ids.append(subject_id)
    return subject_ids


def ensure_blastn_available():
    try:
        return resolve_blastn()
    except BlastUnavailableError as exc:
        sys.exit(f"Error: {exc}")


def run_blastn(
    primer_name: str,
    primer_seq: str,
    subject_file: str,
    subject_sequences: dict[str, str] | None = None,
    execution: BlastExecution | None = None,
) -> dict:
    """
    Run BLASTn with the primer (query) against the subject FASTA file.
    Accepts partial alignments. For each hit, the aligned query (qseq) and subject (sseq)
    sequences are obtained. The custom mismatch function (accounting for ambiguous bases)
    is applied to the aligned region, and positions of mismatches are recorded.
    Bases in the primer not included in the alignment are also counted as mismatches.
    Percent identity is recalculated over the full primer length.
    Returns a dictionary mapping each subject sequence ID (sseqid) to the best alignment.
    """
    execution = execution or BlastExecution()
    try:
        executable = execution.executable or resolve_blastn()
    except BlastUnavailableError:
        if execution.strict_errors:
            raise
        sys.exit("Error: blastn command not found. Please install BLAST+ and ensure blastn is in your PATH.")
    timeout = None
    if execution.deadline is not None:
        timeout = execution.deadline - time.monotonic()
        if timeout <= 0:
            raise AnalysisTimeoutError("Analysis timed out. Try fewer files or select a single assay.")
    query_filename = ""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_query:
        tmp_query.write(f">{primer_name}\n{primer_seq}\n")
        query_filename = tmp_query.name

    blast_command = [
        executable,
        "-query", query_filename,
        "-subject", subject_file,
        "-reward", "2",
        "-penalty", "-3",
        "-word_size", "4",
        "-dust", "yes",
        "-max_target_seqs", BLAST_MAX_TARGET_SEQS,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
    ]

    try:
        try:
            result = subprocess.run(
                blast_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                timeout=timeout,
                env=blast_environment(executable),
            )
        except FileNotFoundError:
            if execution.strict_errors:
                raise BlastUnavailableError("BLAST executable unavailable.") from None
            sys.exit("Error: blastn command not found. Please install BLAST+ and ensure blastn is in your PATH.")
        except subprocess.TimeoutExpired:
            raise AnalysisTimeoutError("Analysis timed out. Try fewer files or select a single assay.") from None
        except OSError:
            if execution.strict_errors:
                raise BlastUnavailableError("BLAST executable could not start.") from None
            raise
    finally:
        if query_filename and os.path.exists(query_filename):
            os.remove(query_filename)

    if result.returncode != 0:
        if execution.strict_errors:
            raise BlastError("BLAST analysis failed. Check the input or try a smaller analysis.")
        print(f"BLASTn error for primer '{primer_name}' against {subject_file}:\n{result.stderr}", file=sys.stderr)
        return {}

    if subject_sequences is None:
        subject_sequences = read_fasta_sequences(subject_file)

    hits_by_subject = {}
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 14:
            continue
        (qseqid, sseqid, pident, align_length, mismatches, gapopen,
         qstart, qend, sstart, send, evalue, bitscore, qseq, sseq) = parts

        try:
            pident = float(pident)
            align_length = int(align_length)
            gapopen = int(gapopen)
            qstart = int(qstart)
            qend = int(qend)
            sstart = int(sstart)
            send = int(send)
            evalue = float(evalue)
            bitscore = float(bitscore)
        except ValueError:
            continue

        query_alignment, subject_alignment = reconstruct_full_primer_alignment(
            primer_seq,
            qseq,
            sseq,
            qstart,
            qend,
            sstart,
            send,
            subject_sequences.get(sseqid, ""),
        )

        adjusted_mismatches = count_mismatches(query_alignment, subject_alignment)
        adjusted_alignment_length = len(primer_seq)
        adjusted_percent_identity = ((adjusted_alignment_length - adjusted_mismatches) / adjusted_alignment_length) * 100

        mismatch_positions = get_mismatch_positions(query_alignment, subject_alignment, 1, len(primer_seq), primer_seq)
        mismatch_details = get_mismatch_details(query_alignment, subject_alignment, 1, len(primer_seq), primer_seq)

        current_hit = {
            "qseqid": qseqid,
            "sseqid": sseqid,
            "pident": adjusted_percent_identity,
            "alignment_length": adjusted_alignment_length,
            "mismatches": adjusted_mismatches,
            "gapopen": gapopen,
            "qstart": qstart,
            "qend": qend,
            "sstart": sstart,
            "send": send,
            "evalue": evalue,
            "bitscore": bitscore,
            "mismatch_positions": mismatch_positions,
            "mismatch_details": mismatch_details,
            "query_alignment": query_alignment,
            "subject_alignment": subject_alignment,
        }

        if sseqid not in hits_by_subject:
            hits_by_subject[sseqid] = current_hit
        else:
            stored_hit = hits_by_subject[sseqid]
            if (current_hit["mismatches"] < stored_hit["mismatches"] or
                (current_hit["mismatches"] == stored_hit["mismatches"] and current_hit["bitscore"] > stored_hit["bitscore"])):
                hits_by_subject[sseqid] = current_hit
    return hits_by_subject

def filter_subject_ids_for_primer(subject_ids: list[str], primer: PrimerRecord, virus_type: str, fasta_file: str, quiet: bool = False) -> list[str]:
    if not virus_type.upper().startswith("INFLUENZA") or not primer.segment:
        return subject_ids

    filtered_subject_ids = []
    unparseable_count = 0
    for subject_id in subject_ids:
        subject_segment = get_segment(subject_id)
        if subject_segment == primer.segment:
            filtered_subject_ids.append(subject_id)
        elif not subject_segment:
            unparseable_count += 1

    if unparseable_count and not quiet:
        print(
            f"Warning: {unparseable_count} subject header(s) in {fasta_file} had no parseable influenza segment "
            f"while filtering primer '{primer.name}' for segment {primer.segment}.",
            file=sys.stderr,
        )
    return filtered_subject_ids


def process_fasta_file(
    fasta_file: str,
    virus_type: str,
    primers: list[PrimerRecord],
    metadata_records: list[MetadataRecord] | None = None,
    *,
    execution: BlastExecution | None = None,
) -> list:
    """
    For a given FASTA file and virus type, run BLASTn for each primer.
    Returns a list of dictionaries (one per subject per primer) with the alignment results.
    For Influenza, if a primer record has an intended segment, only subject
    sequences with that segment are reported.
    """
    execution = execution or BlastExecution()
    results = []
    if not primers:
        print(f"No primer information available for virus type '{virus_type}'.", file=sys.stderr)
        return results

    subject_sequences = read_fasta_sequences(fasta_file)
    subject_ids = get_subject_ids(fasta_file)
    metadata_matches = build_metadata_matches(subject_ids, metadata_records or [], fasta_file, quiet=execution.quiet)

    for primer in primers:
        primer_seq = primer.sequence.upper()
        if not execution.quiet:
            print(f"Running BLASTn for primer '{primer.name}' on file '{fasta_file}' ...")

        filtered_subject_ids = filter_subject_ids_for_primer(subject_ids, primer, virus_type, fasta_file, quiet=execution.quiet)
        hits_by_subject = run_blastn(primer.name, primer_seq, fasta_file, subject_sequences=subject_sequences, execution=execution)

        for subject in filtered_subject_ids:
            subject_segment = get_segment(subject)
            assay_columns = {
                "Assay_Type": primer.assay_type,
                "Assay_ID": primer.scheme_id,
                "Assay_Name": primer.assay_name,
                "Primer_Pool": primer.pool,
                "Primer_Start": primer.start,
                "Primer_End": primer.end,
            }
            if subject in hits_by_subject:
                hit = hits_by_subject[subject]
                result_row = {
                    "Fasta_File": os.path.basename(fasta_file),
                    "Virus_Type": virus_type,
                    "Primer_Name": primer.name,
                    "Primer_Role": primer.role,
                    "Primer_Sequence": primer_seq,
                    "Primer_Segment": primer.segment,
                    "Subject_Sequence_ID": hit["sseqid"],
                    "Subject_Segment": subject_segment,
                    "Hit_Status": "hit",
                    "Percent_Identity": hit["pident"],
                    "Alignment_Length": hit["alignment_length"],
                    "Mismatches": hit["mismatches"],
                    "Gap_Openings": hit["gapopen"],
                    "Query_Start": hit["qstart"],
                    "Query_End": hit["qend"],
                    "Subject_Start": hit["sstart"],
                    "Subject_End": hit["send"],
                    "E_value": hit["evalue"],
                    "Bitscore": hit["bitscore"],
                    "Mismatch_Positions": hit["mismatch_positions"],
                    "Mismatch_Details": hit["mismatch_details"],
                    "Query_Alignment": hit["query_alignment"],
                    "Subject_Alignment": hit["subject_alignment"],
                }
            else:
                result_row = {
                    "Fasta_File": os.path.basename(fasta_file),
                    "Virus_Type": virus_type,
                    "Primer_Name": primer.name,
                    "Primer_Role": primer.role,
                    "Primer_Sequence": primer_seq,
                    "Primer_Segment": primer.segment,
                    "Subject_Sequence_ID": subject,
                    "Subject_Segment": subject_segment,
                    "Hit_Status": "no_hit",
                    "Percent_Identity": "No hit",
                    "Alignment_Length": "",
                    "Mismatches": "",
                    "Gap_Openings": "",
                    "Query_Start": "",
                    "Query_End": "",
                    "Subject_Start": "",
                    "Subject_End": "",
                    "E_value": "",
                    "Bitscore": "",
                    "Mismatch_Positions": "",
                    "Mismatch_Details": "",
                    "Query_Alignment": "",
                    "Subject_Alignment": "",
                }
            result_row.update(assay_columns)
            result_row.update(metadata_columns_for_subject(subject, metadata_matches, virus_type, primer))
            results.append(result_row)
    return results

def write_csv_rows(results: list, stream):
    """Serialize the canonical CSV columns to a file or in-memory text stream."""
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)


def write_csv_report(results: list, output_file: str):
    """Write the results to a CSV file."""
    if not results:
        print("No results to write.", file=sys.stderr)
        return

    try:
        with open(output_file, "w", newline="") as csvfile:
            write_csv_rows(results, csvfile)
        print(f"Report successfully written to {output_file}")
    except Exception as e:
        print(f"Error writing CSV file: {e}", file=sys.stderr)
