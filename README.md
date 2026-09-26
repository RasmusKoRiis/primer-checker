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

## CLI usage
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


## Web application

The browser application uses **Next.js, React, and TypeScript** for uploads,
configuration, native result tables, and alignment inspection. A **FastAPI**
endpoint calls the same functions in `primer_analysis.py` as the CLI. The CLI,
metadata matching, IUPAC mismatch definitions, influenza selection, and existing
standalone report generator remain supported.

### Local development

Use Python 3.12 and Node 22 (see `.python-version` and `.nvmrc`). Install native
NCBI BLAST+ and check that `blastn -version` works first. The CLI itself still
requires only Python 3.10+ and BLAST; the dependencies below are for web/testing.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-dev.txt
npm ci
npm run dev
```

Open [localhost:3000](http://localhost:3000). The single command starts Next.js
and FastAPI on port 8000; Ctrl-C stops both. Upload FASTA files or choose **Use a
synthetic example**, select a virus and assay, and analyze. All files in one
request use the same virus/subtype/assay selection.
For mixed-virus folder routing, use the existing batch CLI.

The UI provides primer/sample tables, sorting, primer/sample/segment/assay/
mismatch/status filters, individual alignment inspection, CSV and standalone
HTML downloads, and a JSON provenance download. Tables aggregate the filtered
comparisons; no-hit results are separate from measured mismatches. The new UI
uses descriptive counts, while the existing standalone HTML retains its
existing risk interpretation. No new scientific risk thresholds are introduced.
The website accepts sequence files and an optional primer database; sample
metadata uploads and Ct/date fields are omitted from its interface. The CLI/API
still support metadata for existing integrations.

Provenance records the analysis time, application/database versions, loaded
primer-record fingerprint, selections, input file hashes, BLAST
version and parameters, and Git revision where supplied. Web CSVs preserve the
canonical analysis columns and append provenance columns; standalone CLI CSV
columns are unchanged. Web CSV text cells that could execute spreadsheet
formulas are prefixed with an apostrophe. Raw JSON/HTML analytical values remain
unchanged.

### Documentation and language

The **Documentation** tab opens the built-in analysis and report guide without
leaving the application. It explains primer selection, BLAST matching, IUPAC
compatibility, full-primer identity and mismatch calculations, no-hit results,
website tables, and the downloaded HTML report's charts and review levels.
Links such as `/#documentation-method` and `/#documentation-levels` open a
specific section directly. Returning to Analysis preserves files, database
editor contents, results, and filters.

**FASTA header rules and examples** appears beside the sequence upload. It shows
ordinary and influenza headers, explains that only the first word after `>` is
the identifier, and includes downloadable synthetic FASTA templates. Influenza
requires a pipe-separated segment tag such as `01-HA|sample_001` before any
spaces; `sample_001 HA` and `sample_001|HA` do not identify the segment.

**View primer database format and download templates** is available below the
database source buttons. It previews complete editable JSON examples for PCR/NGS,
influenza, and the legacy dictionary format, with a field guide and downloads.
Replace the synthetic primers with your own before use; the PCR/NGS template
defaults to `assay_type: "pcr"`, which can be changed to `"ngs"`. The same guides
are in Documentation at `/#documentation-fasta` and `/#documentation-database`.
Preview and download content share `app/lib/input-examples.json`, whose examples
are checked against the real upload validator and preflight endpoint in tests.

The **English / Norsk** buttons switch the analysis form, database builder,
results, and documentation between English and Norwegian Bokmål. The preference
is remembered in browser local storage; no sequence data is stored there.
Language changes preserve the current analysis and only affect presentation.
Uploaded names, sequence data, CSV column names, and calculations stay unchanged.
The downloaded HTML report retains its own English/Norsk switch.

Report guidance is imported directly from `report_text/english.json` and
`report_text/norwegian.json`, so the website and generated report share that
wording. Website-specific text is localized in `app/lib/norwegian.json`, using
English phrases as keys and preserving `{placeholder}` names. Unrecognized
technical API diagnostics fall back to their original wording.

### Your own primer database

Choose **Upload JSON** to use a self-contained primer database for this analysis.
The web accepts legacy `{organism: {primer_name: sequence}}` dictionaries and
normalized schema `1.0` databases with inline `schemes[].primers[]` sequences.
Existing files that reference BED/FASTA assets (`panels[]` or `viruses[]`) remain
supported by the CLI and the installed reference library, not by public uploads.

Choose **Build a database**, enter a name, organism, version, PCR/NGS type, and
primer names/sequences (5′ → 3′), then **Create database**. Add primers with
**Add primer**. Role, segment, pool, and influenza A subtype tags are supported.
Influenza uses organism `Influenza-A` or `Influenza-B` and requires HA/M/NS
segments, matching the existing engine. The validated database becomes active
immediately. **Download database JSON** saves it for later uploads or CLI use:

```bash
python primer_checker.py --primers custom-primers.json --virus SARS-CoV-2 \
  --fasta sample.fasta --output results.csv
```

Exported schemes use optional `assay_type: "pcr" | "ngs"`; omitted values retain
the original PCR behavior. This supports self-contained NGS primer lists without
external BED assets. Database uploads are limited to 250,000 bytes, 500 primers,
and 200 bases per primer. They never replace the reference database and are not
saved on the server. The result provenance records the uploaded filename, its
file hash, and the fingerprint of the primer records actually used. Download
both the custom database and results if you need to reproduce an analysis later.

### API

- `GET /api/catalog`: available viruses, influenza subtypes, assays, database
  fingerprint/version, and upload limits, derived from the installed database.
- `GET /api/health`: verifies the database and that the BLAST executable runs.
- `POST /api/database`: multipart `database` JSON; validates an uploaded library
  and returns its catalog without saving it or running BLAST.
- `POST /api/preflight`: the same multipart fields as analysis; returns exact
  record/comparison/search counts and warnings without running BLAST. Rejects
  invalid or oversized workloads before they can start.
- `POST /api/analyze`: multipart `files` (repeat for multiple FASTAs), optional
  `metadata` CSV, optional `database` JSON, `virus`, optional `flu_type`, `assay_type` (`pcr`, `ngs`, `all`;
  web default `pcr`), and optional exact `assay_id`.

```bash
curl --fail-with-body http://localhost:8000/api/analyze \
  -F 'files=@public/example.fasta' \
  -F 'virus=SARS-CoV-2' \
  -F 'assay_type=pcr'
```

The stateless JSON response contains `rows`, `summary`, `manifest`, `warnings`,
and `downloads` (`csv`, `html`, `manifest` strings). The browser creates downloads
locally; there is no persistent analysis ID or download store. Errors have the
shape `{"error":{"code":"invalid_fasta","message":"..."}}`, with appropriate
422, 413, 502, 503, or 504 status codes. Avoid saving real biological request
payloads or responses in shared logs.

### Limits and privacy

Use public, synthetic, or anonymized consensus sequences. **Do not upload
confidential, identifiable, or otherwise restricted data to this public
deployment.** No suitability for confidential NIPH/surveillance data is claimed.
Uploads use temporary storage and are removed after processing. Browser results
are lost on reload unless downloaded.

The web supports up to 10 FASTAs, 3 MB combined upload data (FASTA and database
files in the UI; also metadata for direct API callers),
200 sequence records, 2,000 primer/record comparisons, 300 BLAST searches, and
50 million total sequence bases × selected primers. A 240-second analysis
deadline and a separate response-size limit also apply. The browser blocks
oversized files before upload, automatically checks valid-sized batches with
preflight, and enables Analyze only for the current validated selection. The
analysis endpoint repeats every check, including for direct API callers.

These conservative per-analysis limits target small workloads on Vercel Hobby.
They do not enforce an account-wide monthly quota or protect against repeated
requests; watch project usage in Vercel. Larger jobs belong in the CLI.
The web does not accept ZIP,
FASTQ, gapped sequences, or duplicate identifiers within a FASTA. Metadata sent
directly to the API uses the existing column matching rules above. Large
NGS/surveillance workloads may need fewer files/one panel at a time, or the CLI.

### Verification

```bash
source .venv/bin/activate
python -m pytest -q
python -m ruff check api web_service tests/web scripts/package_blast.py
npm test
npm run lint
npm run typecheck
npm run build
npx playwright install chromium
npm run test:e2e
```

The real-BLAST equivalence tests are skipped if BLAST is unavailable; install
BLAST for full validation. The end-to-end suite requires it, starts the Python
API and production Next.js server, uploads the synthetic fixture, checks a
mismatch, downloads reports, builds/downloads/re-uploads a custom database,
checks rejection of oversized workloads, and verifies mobile layout. GitHub Actions uses
Python 3.12, Node 22, and the packaged Linux BLAST binary. It also verifies
BLAST execution inside Amazon Linux 2023.

To inspect the production build manually, run the two processes separately:

```bash
# Terminal 1, with the virtual environment active
python -m uvicorn api.index:app --host 127.0.0.1 --port 8000 --no-access-log
# Terminal 2
npm run build
npm run start
```

## Deployment

The repository root is prepared for a Vercel Next.js project with a Python ASGI
function and pinned Linux x86_64 BLAST bundle. No database, permanent upload
storage, or paid external service is required. The production branch must remain
`main`; `feat/webapp-vercel` is for previews and review.

See [the deployment guide](docs/deployment.md) for the exact GitHub/Vercel
connection steps, account limitations, build commands, optional environment
variables, current Vercel limits, BLAST checksums/libraries, remaining preview
verification, fallback compute design, and future custom domain setup. The
implementation environment has no Vercel login, so a hosted preview and project
connection still require the account owner. No production or DNS changes have
been made.

The discovery baseline and design decisions are recorded in
[the implementation note](docs/web-implementation.md).
