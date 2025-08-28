# Primer BLAST Script

## Overview

`primer_blast.py` is a Python script that uses BLASTn to analyze primer alignments in provided FASTA files. It is designed to:

- **Load virus-specific primer sets** from an external JSON file.
- **Run BLASTn** (via the BLAST+ CLI) to align primers (as queries) against subject sequences.
- **Accept partial alignments** and penalize missing bases.
- **Handle ambiguous nucleotide codes** (IUPAC) when counting mismatches.
- **Report detailed alignment metrics**, including percent identity, mismatch count, and the exact positions of mismatches.
- **Filter subject sequences** (e.g., for Influenza, only processing sequences matching a specific segment like HA or M).

The final output is a CSV report that can be visualized with tools like PowerBI.

## Features

- **External Primer Library:** Primer sets are loaded from a JSON file (default: `primers.json`), allowing easy updates without modifying the script.
- **Ambiguous Nucleotide Handling:** Custom functions account for IUPAC ambiguity codes so that bases like `Y` or `R` are compared correctly.
- **Detailed Mismatch Reporting:** In addition to counting mismatches, the script records the positions (1-indexed relative to the primer) where mismatches occur.
- **Flexible Input:** Processes one or more FASTA files containing subject sequences.
- **Influenza-Specific Filtering:** For Influenza virus, an additional parameter (`--flu-type`) specifies whether to use the Influenza-A or Influenza-B primer set, and sequences are filtered by segment (e.g., HA or M) based on FASTA header formatting.

## Prerequisites

- **Python 3**
- **BLAST+ Tools:** Ensure `blastn` is installed and available in your system's PATH.
- **Primer Library File:** A JSON file (e.g., `primers.json`) containing your virus-specific primer sets.

## Sample Primer Library (`primers.json`)

```json
{
  "SARS-CoV-2": {
    "Primer_1": "CTGCAGATTTGGATGATTTCTCC"
  },
  "Influenza-A": {
    "Primer1": "CAAGACCAATCYTGTCACCTCTGAC"
  },
  "Influenza-B": {
    "Primer1": "AGACCAGAGGGAAACTATGCCC"
  },
  "RSV-A": {
    "Primer1": "GACCRATCCTGTCACCTCTGAC"
  },
  "RSV-B": {
    "Primer1": "GACCRATCCTGTCACCTCTGAC"
  }
}
```

## Usage
Run the script from the command line with the required parameters. For example:

```bash
python3 primer_checker.py --primers primers.json --virus influenza --flu-type A --fasta file1.fasta file2.fasta --output primer_report.csv
```

## Command-Line Arguments

- **--primers:** Path to the primer library file in JSON format (default: `primers.json`).
- **--virus:** The virus type to process (e.g., `SARS-CoV-2`, `influenza`, `RSV-A`, `RSV-B`).
- **--flu-type:** For Influenza, specify the subtype: `A` or `B`.
- **--fasta:** One or more FASTA files containing the subject sequences.
- **--output:** The output CSV file for the report (default: `primer_report.csv`).

## How It Works

### Loading the Primer Library
- The script loads the primer sets from the provided JSON file.

### Parsing Subject Sequences
- It reads each FASTA file, extracting subject sequence IDs.
- For Influenza, the script filters sequences by segment (e.g., `HA` or `M`) using a predefined header format.

### BLASTn Alignment
- For each primer, BLASTn is executed with the primer as the query and the subject sequences as the database.
- Partial alignments are accepted.

### Mismatch Analysis
- **Mismatches:** Mismatches in the aligned region are recalculated using a custom function that handles ambiguous bases.
- **Missing Bases:** Missing bases (due to partial alignment) are penalized.
- **Mismatch Positions:** The script records the exact positions of mismatches relative to the full primer.

### Report Generation
- A CSV file is generated containing detailed results, including:
  - FASTA file name
  - Virus type and primer name/sequence
  - Subject sequence ID
  - Percent identity, alignment length, mismatches, gap openings, query/subject start-end positions, e-value, bitscore, and mismatch positions

## Notes

- **BLAST+ Requirement:** Ensure that the `blastn` command is available in your system's PATH.
- **FASTA Header Format:** The script expects a specific format in FASTA headers (e.g., `"03-M|252500127"`) for segment extraction, which is crucial for Influenza processing.
- **Customization:** You can update the primer library JSON file and adjust BLASTn parameters in the script as needed.

## License

This script is provided as-is without any warranty. You are free to modify and distribute it as necessary.

## Contact

For questions or feedback, please open an issue or contact the maintainer.

