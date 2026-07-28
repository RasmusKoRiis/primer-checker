# Primer Checker

## Overview

`primer_checker.py` is a Python script that uses BLASTn to analyze primer alignments in provided FASTA files. It is designed to:

- **Load virus-specific primer sets** from an external JSON file.
- **Validate legacy primer JSON files** before analysis.
- **Run BLASTn** (via the BLAST+ CLI) to align primers (as queries) against subject sequences.
- **Accept partial alignments** and penalize missing bases.
- **Handle ambiguous nucleotide codes** (IUPAC) when counting mismatches.
- **Report detailed alignment metrics**, including percent identity, mismatch count, and the exact positions of mismatches.
- **Filter subject sequences** (e.g., for Influenza, only processing sequences matching a specific segment like HA, M, or NS).
- **Write an optional local HTML report** for browser-based investigation of organisms, primers, samples, hit status, percent identity, and mismatches.

The canonical output is a CSV report that can be visualized with tools like PowerBI. The optional HTML report is intended for local review after a run.

## Project Layout

- `primer_analysis.py`: primer loading, metadata matching, FASTA parsing, BLAST analysis, and CSV output.
- `primer_report.py`: self-contained HTML report generation.
- `report_text/`: editable English and Norwegian wording used in the HTML report.
- `primer_checker.py`: command-line entrypoint and compatibility import surface.
- `scripts/run_primer_checker_batch.py`: batch runner for mixed FASTA folders.
- `primer_db/`: primer databases. `primer_db/fhi_primers.unified.json` is the preferred FHI database, organized by virus with PCR schemes and NGS panels separated.
- `fixtures/`: small test data and metadata examples.
- `docs/`: maintenance notes and implementation plans.
- `result/`: local generated CSV/HTML outputs.


## Editing the HTML report text

The visible wording is stored separately from the Python and JavaScript code:

- Edit `report_text/english.json` for English text.
- Edit `report_text/norwegian.json` for Norwegian text.

Change only the values on the right-hand side of the keys. Keep the key names, JSON commas and quotation marks, and placeholders such as `{sample}`, `{percent}`, `{count}`, `{visible}`, and `{total}` unchanged. Both files must contain the same keys.

After saving, rerun the normal primer-checker command to generate a new HTML report. Existing HTML files do not update automatically. See `report_text/README.md` for the short editing guide.

## Features

- **External Primer Library:** Primer sets are loaded from a JSON file (default: `primers.json`), allowing easy updates without modifying the script.
- **Primer Library Validation:** `--validate-primers` checks legacy and normalized primer JSON structure and sequence alphabet without requiring BLAST or FASTA input.
- **Normalized Primer Database:** The tool accepts the legacy `{organism: {primer_name: sequence}}` format, flat `schemes[]` / `panels[]` databases, and the preferred virus-organized `viruses[]` format with separate `pcr.schemes[]` and `ngs.panels[]` sections.
- **Ambiguous Nucleotide Handling:** Custom functions account for IUPAC ambiguity codes so that bases like `Y` or `R` are compared correctly.
- **Detailed Mismatch Reporting:** In addition to counting mismatches, the script records the positions (1-indexed relative to the primer) where mismatches occur.
- **Flexible Input:** Processes one or more FASTA files containing subject sequences.
- **Influenza-Specific Filtering:** For Influenza virus, an additional parameter (`--flu-type`) specifies whether to use the Influenza-A, H1, H3, or Influenza-B primer set. H1/H3 runs use primers tagged for that subtype plus untagged Influenza-A primers, and exclude primers tagged for the other subtype. Sequences are filtered by segment (e.g., HA, M, or NS) based on FASTA header formatting.
- **Batch Folder Wrapper:** `scripts/run_primer_checker_batch.py` can scan a folder of FASTA files, infer the correct virus/subtype from each filename, and write one combined CSV/HTML report.
- **Optional Sample Metadata:** `--metadata-csv` attaches sample date and assay-specific Ct values to result rows when metadata sample IDs match FASTA headers.
- **Visual Investigation Report:** `--html-report` writes a self-contained HTML/CSS/JS file that opens directly in a browser and provides filters, primer risk badges, stacked sample-nucleotide percentage charts, mismatch count distributions, detailed documentation, English/Norwegian text toggle, and clickable sample alignment popups.

## Prerequisites

- **Python 3**
- **BLAST+ Tools:** Ensure `blastn` is installed and available in your system's PATH.
- **Primer Library File:** A JSON file (e.g., `primers.json`) containing your virus-specific primer sets.
- **pytest:** Required only for running the automated tests.

## Sample Primer Library (`primers.json`)

Legacy format:

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

Normalized format:

```json
{
  "schema_version": "1.0",
  "database_version": "2026-05-06",
  "schemes": [
    {
      "scheme_id": "fhi-influenza-b",
      "display_name": "FHI Influenza-B primers",
      "organism": "Influenza-B",
      "version": "2026-05-06",
      "status": "current",
      "source": "FHI",
      "references": [],
      "primers": [
        {
          "id": "triplex_InfB_F_NS",
          "name": "triplex_InfB_F_NS",
          "sequence": "TCCTCAAYTCACTCTTCGAGCG",
          "role": "forward_primer",
          "segment": "NS",
          "gene": "NS",
          "pool": "triplex",
          "strand": "",
          "subtype_tags": [],
          "notes": ""
        }
      ]
    }
  ]
}
```

Preferred virus-organized format with PCR schemes and BED/FASTA-backed NGS panels:

```json
{
  "schema_version": "3.0",
  "database_version": "2026-05-21",
  "viruses": [
    {
      "organism": "SARS-CoV-2",
      "pcr": {
        "schemes": [
          {
            "scheme_id": "fhi-sars-cov-2",
            "display_name": "FHI SARS-CoV-2 primers",
            "organism": "SARS-CoV-2",
            "version": "2026-05-06",
            "primers": []
          }
        ]
      },
      "ngs": {
        "panels": [
          {
            "panel_id": "sars2-ngs-vmidt-2.2",
            "display_name": "SARS-CoV-2 VMIDT 2.2 NGS Panel",
            "organism": "SARS-CoV-2",
            "technology": "amplicon_ngs",
            "panel_version": "VMIDT.2.2",
            "reference": {
              "name": "NC_045512.2",
              "coordinate_system": "0-based BED"
            },
            "files": {
              "bed": "assets/SARS-CoV-2/VMIDT.2.2/SARS-CoV-2.scheme.bed",
              "primers_fasta": "assets/SARS-CoV-2/VMIDT.2.2/SARS-CoV-2.primers.fasta"
            },
            "mapping": {
              "bed_name_field": "name",
              "pool_from_bed_field": 4
            }
          }
        ]
      }
    }
  ]
}
```

The same database JSON can contain PCR primer `schemes[]` and NGS `panels[]`. The `files` paths are resolved relative to the JSON database file. A panel can load primer sequences from an ARTIC-style primer BED with sequence in column 7, from `primers_fasta`, or from both. When both are present, the FASTA sequence is used for matching and matching BED rows provide pool, strand, coordinate, and reference metadata. BED field numbers in `mapping` are zero-based, so `pool_from_bed_field: 4` reads the fifth BED column.

Use `--assay-type pcr` or `--assay-type ngs` to select a technology. Use `--assay-id` with a `scheme_id` or `panel_id` when the virus has multiple assays. The default `--assay-type all` preserves the previous behavior.

```bash
python3 primer_checker.py \
  --primers primer_db/fhi_primers.unified.json \
  --virus SARS-CoV-2 \
  --assay-type ngs \
  --assay-id sars2-ngs-vmidt-2.2 \
  --fasta sample.fasta \
  --output primer_report.csv
```

## Usage
Run the script from the command line with the required parameters. For example:

```bash
python3 primer_checker.py --primers primers.json --virus influenza --flu-type A --fasta file1.fasta file2.fasta --output primer_report.csv
```

Generate both CSV and local HTML output:

```bash
python3 primer_checker.py --primers primers.json --virus influenza --flu-type B --fasta file1.fasta --output primer_report.csv --html-report primer_report.html
```

Attach one or more previous CSV reports to the HTML report:

```bash
python3 primer_checker.py --primers primers.json --virus influenza --flu-type B --fasta file1.fasta --output primer_report.csv --html-report primer_report.html
```

Attach sample metadata to the CSV and HTML report:

```bash
python3 primer_checker.py \
  --primers primer_db/fhi_primers.normalized.json \
  --virus RSV-A \
  --fasta file1.fasta \
  --metadata-csv fixtures/fake_metadata.csv \
  --output primer_report.csv \
  --html-report primer_report.html
```

Validate the primer database without running BLAST:

```bash
python3 primer_checker.py --primers primers.json --validate-primers
```

Convert a legacy primer database to normalized JSON:

```bash
python3 scripts/convert_legacy_primers.py \
  --input /path/to/fhi_primers.json \
  --output primer_db/fhi_primers.normalized.json \
  --database-version 2026-05-06 \
  --source FHI
```

Run a folder of mixed FASTA files by filename and write combined CSV and HTML reports:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder /path/to/fasta_folder \
  --primers primer_db/fhi_primers.normalized.json \
  --metadata-csv fixtures/fake_metadata.csv \
  --output batch_primer_report.csv \
  --html-report batch_primer_report.html
```

Only run primer analysis and skip the HTML report:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder /path/to/fasta_folder \
  --primers primer_db/fhi_primers.normalized.json \
  --output batch_primer_report.csv \
  --analysis-only
```

Quick test with the small included fixture dataset:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder fixtures/simple_test_data \
  --primers fixtures/simple_test_data/simple_primers.json \
  --metadata-csv fixtures/simple_test_data/metadata.csv \
  --output result/simple_test_report.csv \
  --html-report result/simple_test_report.html
```

Check how filenames will be routed before running BLAST:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder /path/to/fasta_folder \
  --primers primer_db/fhi_primers.normalized.json \
  --dry-run
```

## Command-Line Arguments

- **--primers:** Path to the primer library file in JSON format (default: `primers.json`).
- **--virus:** The virus type to process (e.g., `SARS-CoV-2`, `influenza`, `RSV-A`, `RSV-B`).
- **--assay-type:** Analyze `pcr`, `ngs`, or `all` records for the selected virus (default: `all`).
- **--assay-id:** Restrict analysis to one exact PCR `scheme_id` or NGS `panel_id`.
- **--flu-type:** For Influenza, specify the subtype: `A`, `H1`, `H3`, or `B`.
- **--fasta:** One or more FASTA files containing the subject sequences.
- **--output:** The output CSV file for the report (default: `primer_report.csv`).
- **--html-report:** Optional self-contained HTML report file.
- **--metadata-csv:** Optional sample metadata CSV. The file must contain a sample ID column such as `SampleID` or `Sample_ID`, and may contain `Sample_Date` plus one or more Ct assay columns such as `CT_H3`, `CT_H1`, `Triplex-InfA_CT`, `Triplex-InfB_CT`, `Triplex-SC2_CT`, `CT_RSVA`, or `CT_RSVB`.
- **--validate-primers:** Validate the primer JSON and exit without requiring `--virus`, `--fasta`, or BLAST.

## Optional Metadata CSV

Metadata is attached to each result row when a metadata sample ID matches the FASTA subject header. Matching is strict enough to avoid unsafe substring matches:

- `SampleID` `454511` matches FASTA header `Genome|454511`.
- `SampleID` `4545` does not match FASTA header `Genome|454511`.
- `SampleID` `Genome|4545` matches FASTA header `Genome|4545`, but not `Genome|4546`.

Recognized metadata columns:

- sample ID: `SampleID`, `Sample_ID`, `Sample`, or `ID`
- date: `Sample_Date`, `SampleDate`, `Date`, `Collection_Date`, or `Sampling_Date`
- Ct: generic `Ct`, `CT`, `Ct_Value`, `Cq`, `Cp`, or `Cycle_Threshold`
- assay-specific Ct: any column starting or ending with `CT`, `CQ`, or `CP`, for example `CT_H3`, `CT_H1`, `CT_INFA`, `CT_INFB`, `Triplex-InfA_CT`, `Triplex-InfB_CT`, `CT_RSVA`, `CT_RSVB`, or `Triplex-SC2_CT`

Attached metadata is written to the CSV columns `Metadata_Sample_ID`, `Sample_Date`, `Ct_Value`, and `Ct_Source`. `Ct_Value` is selected per result row from the best matching assay-specific metadata column:

- H3 rows prefer `CT_H3`.
- H1 rows prefer `CT_H1`.
- Influenza-A triplex primer rows prefer `Triplex-InfA_CT`.
- Influenza-B triplex primer rows prefer `Triplex-InfB_CT`.
- SARS-CoV-2 rows prefer `Triplex-SC2_CT`, `CT_SC2`, or similar SC2/SARS-CoV-2 names.
- RSV-A and RSV-B rows prefer `CT_RSVA` and `CT_RSVB`.
- Influenza-A and Influenza-B rows prefer `CT_INFA` and `CT_INFB`.
- Generic `Ct_Value` is used as a fallback.

ISO dates such as `2025-05-07` are preferred for CSV consumers. The HTML report currently omits sample metadata fields while the surveillance-focused report layout is being revised.

A small test metadata file based on the external FASTA sample IDs is included at `fixtures/fake_metadata.csv`.

### Metadata Export From SQL/GISAID R Scripts

The R scripts in `scripts/` now write primer-checker metadata CSVs with the same schema as `fixtures/fake_metadata.csv`:

- `scripts/influenza_gisaid.R` writes `INFLUENZA PRIMER CHECKER METADATA - <week-year>.csv`
- `scripts/rsv_gisaid.R` writes `RSV PRIMER CHECKER METADATA - <week-year>.csv`
- `scripts/sc2_gisaid.R` writes `SC2 PRIMER CHECKER METADATA - <week-year>.csv`

These scripts no longer make GISAID CSV/Excel/FASTA outputs, and they no longer filter by SID/run ID. They export metadata for the loaded database rows. If a source database does not contain one of the expected Ct columns, that output column is left blank.

Each script accepts an optional first argument for the output directory:

```bash
Rscript scripts/influenza_gisaid.R /path/to/output
Rscript scripts/rsv_gisaid.R /path/to/output
Rscript scripts/sc2_gisaid.R /path/to/output
```

If no output directory is supplied, each script uses its existing default result directory.

For the influenza SQL-derived `fludb` table, the metadata exporter maps:

- `pcr_h1_ct` -> `CT_H1`
- `pcr_h3_ct` -> `CT_H3`
- `pcr_bvic_ct` / `pcr_byam_ct` -> `CT_INFB`
- `prove_tatt` -> `Sample_Date`

Influenza triplex Ct columns are left blank unless triplex Ct fields are added to the source table later.

## Batch Folder Wrapper

`scripts/run_primer_checker_batch.py` is intended for a folder containing several FASTA files from different organisms or subtypes. It supports `.fa`, `.fasta`, `.fna`, and `.fas` files.

Filename routing rules:

- `H3`, `H3N2` -> `--virus influenza --flu-type H3`
- `H1`, `H1N1` -> `--virus influenza --flu-type H1`
- `INFA`, `FLUA` -> `--virus influenza --flu-type A`
- `INFB`, `FLUB` -> `--virus influenza --flu-type B`
- `SC2`, `SARS`, `SARS2`, `COV2`, `COVID`, `COVID19` -> `--virus SARS-CoV-2`
- `RSVA` -> `--virus RSV-A`
- `RSVB` -> `--virus RSV-B`
- Other non-influenza organisms can be routed when the FASTA filename contains the organism name from the loaded primer database, for example `Norovirus-GII_run1.fasta` for organism `Norovirus-GII`.

By default, the batch runner writes both the combined CSV and `batch_primer_report.html`. Use `--analysis-only` to skip HTML report generation. Use `--dry-run` first to verify routing. Use `--fail-on-unclassified` if the run should stop when a filename does not match any rule. Ambiguous `RSV` filenames without `A` or `B` are not classified because the primer database stores RSV-A and RSV-B separately.

For H1/H3 filename-based influenza runs, subtype-specific primer selection is based on primer metadata:

- H3 FASTA files use Influenza-A primers tagged `H3` plus untagged Influenza-A primers.
- H1 FASTA files use Influenza-A primers tagged `H1` plus untagged Influenza-A primers.
- H1 FASTA files do not use H3-tagged primers, and H3 FASTA files do not use H1-tagged primers.

## How It Works

### Loading the Primer Library
- The script loads the primer sets from the provided JSON file.
- Supported database formats are the legacy shape `{organism: {primer_name: sequence}}`, flat `schemes[]` / `panels[]` databases, and the preferred virus-organized `viruses[]` database used by `primer_db/fhi_primers.unified.json`.
- Legacy primer records are converted internally to structured primer records with inferred segment and role metadata where possible.
- Normalized records use explicit metadata where available.
- Panel records load primer sequences from BED or FASTA assets and keep panel ID, panel version, pool, strand, reference, and coordinate metadata internally. In the virus-organized database, PCR records live under `pcr.schemes[]` and NGS panels live under `ngs.panels[]` for each organism.

### Parsing Subject Sequences
- It reads each FASTA file, extracting subject sequence IDs.
- For Influenza, the script filters sequences by segment (e.g., `HA`, `M`, or `NS`) using supported header tokens such as `"03-M|252500127"` or `"contig1|03-M|INFL16-2025"`.

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
  - Primer and subject segment where available
  - Primer/sample alignment strings and mismatch base-change details where BLAST hits are available
  - Attached metadata sample ID, sample date, selected Ct value, and Ct source where available
  - Hit status
  - Subject sequence ID
  - Percent identity, alignment length, mismatches, gap openings, query/subject start-end positions, e-value, bitscore, and mismatch positions
- If `--html-report` is supplied, a local visual report is also generated from the same result rows.

## Testing

Run the focused unit tests with:

```bash
python3 -m pytest -q
```

The unit tests do not require BLAST.

## Notes

- **BLAST+ Requirement:** Ensure that the `blastn` command is available in your system's PATH.
- **FASTA Header Format:** The script supports segmented influenza tokens such as `"03-M|252500127"` and `"contig1|03-M|INFL16-2025"` for segment extraction, which is crucial for Influenza processing.
- **Customization:** You can update the primer library JSON file and adjust BLASTn parameters in the script as needed.

## License

This script is provided as-is without any warranty. You are free to modify and distribute it as necessary.

## Contact

For questions or feedback, please open an issue or contact the maintainer.
