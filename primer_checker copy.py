#!/usr/bin/env python3
"""
primer_blast.py

This script:
  - Reads one or more FASTA files (each containing one or more subject sequences).
  - Uses a virus type (passed as a command-line argument) to select a set of primers loaded from an external JSON file.
  - For Influenza, an additional parameter (--flu-type) specifies whether to use the Influenza-A or Influenza-B primer set.
  - For each primer and each FASTA file, runs a BLASTn search (via the BLAST+ CLI)
    where the primer is the query and the FASTA file is the subject.
  - Accepts partial alignments. For partial alignments, bases at the start and/or end of the primer
    that are not aligned are penalized (counted as mismatches) and mismatches are recalculated
    over the full primer length. In addition, the script records the positions of all mismatches.
  - For each subject (sample) in the FASTA file, reports the best hit (if any) for that primer.
    For Influenza, only subject sequences matching the primer’s intended segment (e.g., HA or M)
    are processed.
  - Writes a CSV report with the results for visualization.

Usage:
    python3 primer_blast.py --primers primers.json --virus influenza --flu-type A --fasta file1.fasta file2.fasta --output primer_report.csv
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

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

# --- End Ambiguity Functions ---

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

def get_segment(subject_id: str) -> str:
    """
    Given a subject FASTA header (without the leading '>'),
    extract the segment code assuming the header format is like:
      "03-M|252500127" or "03-HA|252500127"
    Returns the segment code (e.g., "M" or "HA") or an empty string if not found.
    """
    try:
        header_part = subject_id.split("|")[0]
        parts = header_part.split("-")
        if len(parts) >= 2:
            return parts[1].upper()
    except Exception:
        return ""
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

def run_blastn(primer_name: str, primer_seq: str, subject_file: str) -> dict:
    """
    Run BLASTn with the primer (query) against the subject FASTA file.
    Accepts partial alignments. For each hit, the aligned query (qseq) and subject (sseq)
    sequences are obtained. The custom mismatch function (accounting for ambiguous bases)
    is applied to the aligned region, and positions of mismatches are recorded.
    Bases in the primer not included in the alignment are also counted as mismatches.
    Percent identity is recalculated over the full primer length.
    Returns a dictionary mapping each subject sequence ID (sseqid) to the best alignment.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_query:
        tmp_query.write(f">{primer_name}\n{primer_seq}\n")
        query_filename = tmp_query.name

    blast_command = [
        "blastn",
        "-query", query_filename,
        "-subject", subject_file,
        "-reward", "2",
        "-penalty", "-3",
        "-word_size", "4",
        "-dust", "yes",
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
    ]

    try:
        result = subprocess.run(
            blast_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
    except FileNotFoundError:
        sys.exit("Error: blastn command not found. Please ensure BLAST+ is installed and in your PATH.")

    os.remove(query_filename)

    if result.returncode != 0:
        print(f"BLASTn error for primer '{primer_name}' against {subject_file}:\n{result.stderr}", file=sys.stderr)
        return {}

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

        missing_start = qstart - 1
        missing_end = len(primer_seq) - qend

        custom_mismatches = count_mismatches(qseq, sseq)
        adjusted_mismatches = custom_mismatches + missing_start + missing_end
        adjusted_alignment_length = len(primer_seq)
        adjusted_percent_identity = ((adjusted_alignment_length - adjusted_mismatches) / adjusted_alignment_length) * 100

        mismatch_positions = get_mismatch_positions(qseq, sseq, qstart, qend, primer_seq)

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
            "mismatch_positions": mismatch_positions
        }

        if sseqid not in hits_by_subject:
            hits_by_subject[sseqid] = current_hit
        else:
            stored_hit = hits_by_subject[sseqid]
            if (current_hit["mismatches"] < stored_hit["mismatches"] or
                (current_hit["mismatches"] == stored_hit["mismatches"] and current_hit["bitscore"] > stored_hit["bitscore"])):
                hits_by_subject[sseqid] = current_hit
    return hits_by_subject

def process_fasta_file(fasta_file: str, virus_type: str) -> list:
    """
    For a given FASTA file and virus type, run BLASTn for each primer.
    Returns a list of dictionaries (one per subject per primer) with the alignment results.
    For Influenza, if a primer name indicates an intended segment (HA or M),
    only subject sequences with that segment are processed.
    """
    results = []
    # Use the globally loaded VIRUS_PRIMERS dictionary.
    primers = VIRUS_PRIMERS.get(virus_type)
    if not primers:
        print(f"No primer information available for virus type '{virus_type}'.", file=sys.stderr)
        return results

    subject_ids = get_subject_ids(fasta_file)

    for primer_name, primer_seq in primers.items():
        primer_seq = primer_seq.upper()
        print(f"Running BLASTn for primer '{primer_name}' on file '{fasta_file}' ...")
        
        filtered_subject_ids = subject_ids
        if virus_type.upper().startswith("INFLUENZA"):
            intended_segment = None
            if "HA" in primer_name.upper():
                intended_segment = "HA"
            elif "M" in primer_name.upper():
                intended_segment = "M"
            if intended_segment:
                filtered_subject_ids = [sid for sid in subject_ids if get_segment(sid) == intended_segment]

        hits_by_subject = run_blastn(primer_name, primer_seq, fasta_file)

        for subject in filtered_subject_ids:
            if subject in hits_by_subject:
                hit = hits_by_subject[subject]
                result_row = {
                    "Fasta_File": os.path.basename(fasta_file),
                    "Virus_Type": virus_type,
                    "Primer_Name": primer_name,
                    "Primer_Sequence": primer_seq,
                    "Subject_Sequence_ID": hit["sseqid"],
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
                    "Mismatch_Positions": hit["mismatch_positions"]
                }
            else:
                result_row = {
                    "Fasta_File": os.path.basename(fasta_file),
                    "Virus_Type": virus_type,
                    "Primer_Name": primer_name,
                    "Primer_Sequence": primer_seq,
                    "Subject_Sequence_ID": subject,
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
                    "Mismatch_Positions": ""
                }
            results.append(result_row)
    return results

def write_csv_report(results: list, output_file: str):
    """Write the results to a CSV file."""
    if not results:
        print("No results to write.", file=sys.stderr)
        return

    fieldnames = [
        "Fasta_File",
        "Virus_Type",
        "Primer_Name",
        "Primer_Sequence",
        "Subject_Sequence_ID",
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
        "Mismatch_Positions"
    ]
    try:
        with open(output_file, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"Report successfully written to {output_file}")
    except Exception as e:
        print(f"Error writing CSV file: {e}", file=sys.stderr)

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
        required=True,
        help="Virus type (e.g., 'SARS-CoV-2', 'Influenza', 'RSV-A', 'RSV-B')."
    )
    parser.add_argument(
        "--flu-type",
        type=str,
        choices=["H1", "H3", "B"],
        help="For influenza, specify the type: A or B."
    )
    parser.add_argument(
        "--fasta",
        type=str,
        nargs="+",
        required=True,
        help="One or more FASTA file paths (subject sequences) to process."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="primer_report.csv",
        help="Output CSV filename (default: primer_report.csv)."
    )
    args = parser.parse_args()

    # Load the primer library from the specified JSON file.
    global VIRUS_PRIMERS
    VIRUS_PRIMERS = load_primer_library(args.primers)

    virus_type = args.virus
    if virus_type.lower() == "influenza":
        if not args.flu_type:
            sys.exit("Error: For influenza, please specify --flu-type A or B.")
        virus_type = "Influenza-" + args.flu_type.upper()

    fasta_files = args.fasta
    output_file = args.output

    all_results = []
    for fasta_file in fasta_files:
        if not os.path.exists(fasta_file):
            print(f"FASTA file '{fasta_file}' not found. Skipping.", file=sys.stderr)
            continue
        print(f"Processing FASTA file: {fasta_file}")
        file_results = process_fasta_file(fasta_file, virus_type)
        all_results.extend(file_results)

    write_csv_report(all_results, output_file)

if __name__ == "__main__":
    main()
