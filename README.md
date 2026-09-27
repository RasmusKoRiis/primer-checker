# Primer Checker

## Overview

`primer_checker.py` is a Python script that uses BLASTn to analyze primer alignments in provided FASTA files. It is designed to:

- **Load virus-specific primer sets** from an external JSON file.
- **Validate legacy primer JSON files** before analysis.
- **Run BLASTn** (via the BLAST+ CLI) to align primers (as queries) against subject sequences.
- **Accept partial alignments** and penalize missing bases.
- **Handle ambiguous nucleotide codes** (IUPAC) when counting mismatches.
- **Report detailed alignment metrics**, including percent identity, mismatch count, and the exact positions of mismatches.
- **Filter subject sequences** (for influenza, compare each primer only with FASTA records carrying its database-defined segment label).
- **Write an optional local HTML report** for browser-based investigation of organisms, primers, samples, hit status, percent identity, and mismatches.

The canonical output is a CSV report that can be visualized with tools like PowerBI. The optional HTML report is intended for local review after a run.

## Project Layout

- `primer_analysis.py`: primer loading, metadata matching, FASTA parsing, BLAST analysis, and CSV output.
- `primer_report.py`: self-contained HTML report generation.
- `report_text/`: editable English and Norwegian wording used in the HTML report.
- `primer_checker.py`: command-line entrypoint and compatibility import surface.
- `scripts/run_primer_checker_batch.py`: batch runner for mixed FASTA folders.
- `primer_db/`: clearly labelled synthetic test databases. `primer_db/dummy_primers.json` is the website default; no real assay library is bundled.
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
- **Influenza-Specific Filtering:** `--flu-type` selects a type or subtype actually present in the database, such as `A`, `H5N1`, or `B/VICTORIA`. A type selects all its primers; a subtype selects primers with that exact tag plus untagged primers of the same type. Segment labels are defined by the database and matched to FASTA headers, without a fixed segment list.
- **Batch Folder Wrapper:** `scripts/run_primer_checker_batch.py` can scan a folder of FASTA files, infer the correct virus/subtype from each filename, and write one combined CSV/HTML report.
- **Optional Sample Metadata:** `--metadata-csv` attaches sample date and assay-specific Ct values to result rows when metadata sample IDs match FASTA headers.
- **Visual Investigation Report:** `--html-report` writes a self-contained HTML/CSS/JS file that opens directly in a browser and provides filters, primer risk badges, stacked sample-nucleotide percentage charts, mismatch count distributions, detailed documentation, English/Norwegian text toggle, and clickable sample alignment popups.

## Prerequisites

- **Python 3**
- **BLAST+ Tools:** Ensure `blastn` is installed and available in your system's PATH.
- **Primer Library File:** A JSON file (e.g., `primers.json`) containing your virus-specific primer sets.
- **pytest:** Required only for running the automated tests.

## Sample Primer Library (`primers.json`)

All examples shipped here are **dummy data for software testing only**, not validated assays.
The self-contained [dummy database](primer_db/dummy_primers.json) can be downloaded
from the website and uploaded again. It contains three synthetic PCR oligos,
two NGS oligos, and eight synthetic influenza segment examples. The influenza
organism and subtype labels demonstrate routing; the sequences are invented.

Legacy format:

```json
{"Demo-virus": {"DUMMY_F": "GTCAGACATCGATGCTACGTCAGGATCGTACCTAGCTGAC"}}
```

Normalized format (accepted by the website and CLI):

```json
{
  "schema_version": "1.0",
  "database_version": "dummy-1.0",
  "purpose": "synthetic-test-only",
  "schemes": [{
    "scheme_id": "dummy-pcr",
    "display_name": "DUMMY — PCR example",
    "organism": "Demo-virus",
    "version": "dummy-1.0",
    "assay_type": "pcr",
    "primers": [{
      "id": "DUMMY_F", "name": "DUMMY_F",
      "sequence": "GTCAGACATCGATGCTACGTCAGGATCGTACCTAGCTGAC",
      "role": "forward_primer", "segment": "", "subtype_tags": []
    }]
  }]
}
```

The CLI also accepts virus-organized schema `3.0` databases and external
BED/FASTA panel assets. See the entirely synthetic
[panel example](primer_db/assets/panel_database.example.json).

The same database JSON can contain PCR primer `schemes[]` and NGS `panels[]`. The `files` paths are resolved relative to the JSON database file. A panel can load primer sequences from an ARTIC-style primer BED with sequence in column 7, from `primers_fasta`, or from both. When both are present, the FASTA sequence is used for matching and matching BED rows provide pool, strand, coordinate, and reference metadata. BED field numbers in `mapping` are zero-based, so `pool_from_bed_field: 4` reads the fifth BED column.

Use `--assay-type pcr` or `--assay-type ngs` to select a technology. Use `--assay-id` with a `scheme_id` or `panel_id` when the virus has multiple assays. The default `--assay-type all` preserves the previous behavior.

```bash
python3 primer_checker.py \
  --primers primer_db/dummy_primers.json \
  --virus Demo-virus \
  --assay-type ngs \
  --assay-id dummy-ngs \
  --fasta public/example.fasta \
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
  --primers /path/to/my-primers.json \
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
  --input /path/to/legacy-primers.json \
  --output /path/to/my-primers.json \
  --database-version 2026-05-06 \
  --source "Your laboratory"
```

Run a folder of mixed FASTA files by filename and write combined CSV and HTML reports:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder /path/to/fasta_folder \
  --primers /path/to/my-primers.json \
  --metadata-csv fixtures/fake_metadata.csv \
  --output batch_primer_report.csv \
  --html-report batch_primer_report.html
```

Only run primer analysis and skip the HTML report:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder /path/to/fasta_folder \
  --primers /path/to/my-primers.json \
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
  --primers /path/to/my-primers.json \
  --dry-run
```

## Command-Line Arguments

- **--primers:** Path to the primer library file in JSON format (default: `primers.json`).
- **--virus:** The virus type to process (e.g., `SARS-CoV-2`, `influenza`, `RSV-A`, `RSV-B`).
- **--assay-type:** Analyze `pcr`, `ngs`, or `all` records for the selected virus (default: `all`).
- **--assay-id:** Restrict analysis to one exact PCR `scheme_id` or NGS `panel_id`.
- **--flu-type:** For influenza, specify a type or exact subtype label from the database (for example `A`, `H1`, `H5N1`, or `B/VICTORIA`).
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

Filename routing rules (influenza choices must exist in the loaded database):

- `H3`, `H3N2` -> `--virus influenza --flu-type H3`
- `H1`, `H1N1` -> `--virus influenza --flu-type H1`
- `INFA`, `FLUA` -> `--virus influenza --flu-type A`
- `INFB`, `FLUB` -> `--virus influenza --flu-type B`
- Other exact subtype labels, such as `H5N1` or `H7N9`, are read from the database. Exact full tags take priority over H-family fallback; `H1N1` routes to `H1` only if the database has `H1` and no `H1N1` tag. Conflicting subtype labels remain unclassified. Other influenza types can use markers such as `INFC` when that type is in the database.
- `SC2`, `SARS`, `SARS2`, `COV2`, `COVID`, `COVID19` -> `--virus SARS-CoV-2`
- `RSVA` -> `--virus RSV-A`
- `RSVB` -> `--virus RSV-B`
- Other non-influenza organisms can be routed when the FASTA filename contains the organism name from the loaded primer database, for example `Norovirus-GII_run1.fasta` for organism `Norovirus-GII`.

By default, the batch runner writes both the combined CSV and `batch_primer_report.html`. Use `--analysis-only` to skip HTML report generation. Use `--dry-run` first to verify routing. Use `--fail-on-unclassified` if the run should stop when a filename does not match any rule. Ambiguous `RSV` filenames without `A` or `B` are not classified because the primer database stores RSV-A and RSV-B separately.

For all influenza subtype runs, primer selection is based on database metadata. For example:

- H3 FASTA files use Influenza-A primers tagged `H3` plus untagged Influenza-A primers.
- H1 FASTA files use Influenza-A primers tagged `H1` plus untagged Influenza-A primers.
- H1 FASTA files do not use H3-tagged primers, and H3 FASTA files do not use H1-tagged primers.

## How It Works

### Loading the Primer Library
- The script loads the primer sets from the provided JSON file.
- Supported database formats are the legacy shape `{organism: {primer_name: sequence}}`, flat `schemes[]` / `panels[]` databases, and virus-organized `viruses[]` databases. The bundled dummy database uses self-contained `schemes[]`.
- Legacy primer records are converted internally to structured primer records with inferred segment and role metadata where possible.
- Normalized records use explicit metadata where available.
- Panel records load primer sequences from BED or FASTA assets and keep panel ID, panel version, pool, strand, reference, and coordinate metadata internally. In the virus-organized database, PCR records live under `pcr.schemes[]` and NGS panels live under `ngs.panels[]` for each organism.

### Parsing Subject Sequences
- It reads each FASTA file, extracting subject sequence IDs.
- For influenza, the script matches database segment labels to header tokens such as `"01-PB2|sample"`, `"06-NA|sample"`, or `"contig1|03-M|sample"`. Segment names contain 1–32 letters or digits; matching is case-insensitive. The numeric prefix is not used to infer a segment.

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

### Bundled dummy database

The default is **dummy-1.0**, a synthetic software test set, independent of any
real primer database. Choose **Use a synthetic example** to load two artificial
FASTA records: PCR gives six hits and one mismatch (`9:T>A`). The reverse primer
binds the opposite strand. **Download dummy database** saves the matching JSON.
NGS mode uses two invented oligos against the same example (four hits, one mismatch).
For influenza routing, use `primer_db/dummy-influenza.fasta` with the dummy
influenza scheme. No biological performance is implied by these examples.

The site labels dummy data in English/Norwegian before analysis and in results.
Web CSV exports include `Database_Purpose=synthetic-test-only`; provenance records
`is_dummy: true`, and HTML downloads carry a visible bilingual dummy-data notice.
The JSON's `purpose` marker preserves that label when the dummy database is uploaded.
Remove this marker when replacing all example sequences with your own real library.
For your own analyses, upload or build your own primer database. There is no
runtime dependency on an external GitHub primer library.

The former library, panel assets, archive and saved report were removed from
this branch's current files. Earlier Git commits and other branches/deployments
may still contain them; this change does not rewrite Git history.

### Your own primer database

Choose **Upload a database** to use a self-contained primer database for this analysis.
The web accepts legacy `{organism: {primer_name: sequence}}` dictionaries and
normalized schema `1.0` databases with inline `schemes[].primers[]` sequences.
Existing files that reference BED/FASTA assets (`panels[]` or `viruses[]`) remain
supported by the CLI and an operator-installed database, not by public uploads.

Choose **Build a database**, enter a name, organism, version, PCR/NGS type, and
primer names/sequences (5′ → 3′), then **Create database**. Add primers with
**Add primer**. Role, segment, pool, and influenza subtype tags are supported.
Supply reverse primers as the actual oligo in **5′ → 3′** direction too; do not
reverse-complement them before entry. BLAST searches both subject strands,
independently of role metadata. Alignments show both rows in primer orientation,
reverse-complementing reverse-strand subject hits. Position 1 is always the
primer's 5′ end, and the last five primer bases remain its 3′ terminal region.
The displayed subject coordinates describe the original local BLAST hit
(1-based, inclusive), before extension to the full primer; reverse hits have
decreasing coordinates.

Alignment gaps are preserved. An insertion is anchored to the preceding primer
base without shifting later positions (`19:->T` means T inserted after position
19; a leading insertion uses anchor 1). Each gap column counts as a difference,
but position charts count each affected anchor once per hit. Changes sharing
an anchor are combined in the base-difference field. Regression fixtures in
`fixtures/reverse_primer` cover both strands, terminal mismatches, and indels.

For influenza, the database specifies the organism type (for example
`Influenza-A`, `Influenza-B`, `Influenza-C`, or `Influenza-D`) and any segment
label matching the FASTA headers. The validated database becomes active
immediately. **Download database JSON** saves it for later uploads or CLI use:

```bash
python primer_checker.py --primers custom-primers.json --virus SARS-CoV-2 \
  --fasta sample.fasta --output results.csv
```

Exported schemes use optional `assay_type: "pcr" | "ngs"`; omitted values retain
the original PCR behavior. This supports self-contained NGS primer lists without
external BED assets. Database uploads are limited to 250,000 bytes, 500 primers,
and 200 bases per primer. They never replace the installed database and are not
saved on the server. The result provenance records the uploaded filename, its
file hash, and the fingerprint of the primer records actually used. Download
both the custom database and results if you need to reproduce an analysis later.

Influenza selection is database-driven in both the website and CLI:

- Segment suggestions include `PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `M`, and
  `NS`, but any label with 1–32 letters or digits is valid. A primer with
  `segment: "PB2"` is compared only with headers such as `>01-PB2|sample`.
  Header numbers do not determine the segment, and no gene-name aliases are
  assumed: use matching labels in the database and FASTA.
- Enter comma-separated subtype tags in the builder, or a JSON list such as
  `"subtype_tags": ["H5N1", "H7N9"]`. Labels are normalized to uppercase.
  Matching is exact: `H5` and `H5N1` are distinct selections unless the primer
  explicitly lists both. Up to 32 tags per primer are accepted; each tag is
  1–64 letters/digits/dots/underscores/hyphens, starting with a letter or digit.
- An explicit empty list `[]` (a blank builder field) means the primer is
  shared by all subtypes of that influenza type. If the field is omitted,
  legacy name tokens such as `H5` or `H5N1` may supply tags. Legacy segment
  inference recognizes the usual segment names; use explicit metadata for
  other labels.
- The menu lists only types and tags present in the loaded database. `A` or
  `B` selects every primer for that type. A subtype selects its tagged primers
  plus untagged primers of the same type. A-subtype names retain short IDs
  such as `H1`; other types use qualified IDs such as `B/VICTORIA` to avoid
  mixing organisms. Assay choices and counts reflect the selected subtype.

The downloadable influenza example contains all eight usual segment labels
and both `H1` and `H5N1` primers. It contains synthetic sequences for checking
file format and selection behavior, not biological assay validation.

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

The [live application](https://primer-checker.vercel.app) is deployed on the
**Hobby** plan under `rasmus-projects1` and can be used without a Vercel login.
The tested feature branch was explicitly published for live testing; the draft
PR is still open and no branches were merged. The existing fork,
`RasmusKRiis/primer-checker`, is connected for Git deployments, with `main` as
the production branch and other branches creating protected previews.
PCR/NGS analysis, custom databases, report downloads, and workload rejection
were verified on Vercel. The public browser flow also completes a real analysis.
No custom domain or paid service was added.

See [the deployment guide](docs/deployment.md) for verified results, direct
preview and production updates, the GitHub connection, Hobby limits,
BLAST checksums/libraries, fallback compute design, and future domain setup.

The discovery baseline and design decisions are recorded in
[the implementation note](docs/web-implementation.md).
