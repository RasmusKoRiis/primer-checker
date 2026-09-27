# Primer Checker Improvement Plan

> Historical planning notes from before the web application. References to the
> old external primer library describe that earlier state. The current repository
> ships only synthetic dummy data; see `primer_db/README.md`.


This is a living investigation plan. It records observations before implementation work. The external data and primer database folders were inspected read-only.

## 1. Current understanding

- The repository is a small CLI analysis tool centered on `primer_checker.py`.
- The tool loads a JSON primer library, runs `blastn` for each primer against one or more subject FASTA files, recalculates mismatch metrics with IUPAC ambiguity support, keeps one best hit per subject, and writes a CSV report.
- Main input files:
  - Primer library JSON, default `primers.json`, passed with `--primers`.
  - FASTA files passed with `--fasta`.
  - Virus selector passed with `--virus`.
  - Influenza subtype selector passed with `--flu-type A|H1|H3|B` when `--virus influenza`.
- Main output file:
  - CSV report, default `primer_report.csv`, with one row per primer/subject sequence after any influenza segment filtering.
- Main command:
  - `python3 primer_checker.py --primers <primers.json> --virus <virus> [--flu-type <A|H1|H3|B>] --fasta <files...> --output <report.csv>`
- External realistic data is in `/home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/data`.
- The available primer database is `/home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json`.
- The local environment has Python 3.12.10, but `blastn` and `pytest` were not available during investigation.

## 2. Repository map

- `README.md`: usage and high-level behavior. It refers to `primer_blast.py` in the overview, while the actual script is `primer_checker.py`.
- `primer_checker.py`: complete CLI, primer loading, FASTA header parsing, BLAST invocation, mismatch recalculation, influenza subtype handling, and CSV writing.
- `primer_checker copy.py`: tracked in git but already deleted in the current worktree before this investigation. The deleted tracked version is very similar to the current script and appears to be an older backup.
- `CODEX_IMPROVEMENT_PLAN.md`: this plan.
- No tests, package metadata, dependency lockfile, config file, sample primer JSON, or sample FASTA files are present in the repository root.

## 3. Data observations

- Data folder structure is flat.
- File types found:
  - FASTA: `H1.fasta`, `H3.fasta`, `INFA_test.fasta`, `INFB_test.fasta`, `RSVA_primercheck.fasta`, `RSVB.fasta`, `SC2.fasta`.
  - Python backup: `primer_checker_backup.py`.
- FASTA sequence counts:
  - `H1.fasta`: 386 sequences, HA headers only in sampled/summary parsing.
  - `H3.fasta`: 265 sequences, HA headers only in sampled/summary parsing.
  - `INFA_test.fasta`: 1016 sequences, mostly HA and M plus two contig-style headers.
  - `INFB_test.fasta`: 325 sequences, HA, M, and NS segments.
  - `RSVA_primercheck.fasta`: 78 genome sequences.
  - `RSVB.fasta`: 54 genome sequences.
  - `SC2.fasta`: 1427 nanopore SARS-CoV-2 sequences.
- Naming/header conventions observed:
  - Influenza segmented headers commonly look like `01-HA|252506719`, `03-M|252502109`, or `08-NS|252502109`.
  - Some influenza accessions include strain-like identifiers after the pipe, for example `A/Netherlands/10685/2024`.
  - `INFA_test.fasta` includes `>contig1|03-M|INFL16-2025` and `>contig2|03-M|INFL16-2025`, which the current `get_segment()` will not recognize as M because it only examines text before the first pipe.
  - RSV headers are `Genome|...`.
  - SARS-CoV-2 headers are `Nanopore|...`.
- Sequence observations:
  - FASTA records are multiline.
  - Lengths vary substantially by file, as expected for segments versus whole genomes.
  - One invalid alphabet example was found in `INFA_test.fasta`, associated with the contig header, containing trailing non-IUPAC text. A validator should report this clearly.
- Suspicious edge cases:
  - Influenza B has NS segment records and NS primers, but the current analysis only recognizes HA and M from primer names.
  - The current subject ID extraction keeps only the first whitespace-delimited token, which is fine for these files but should be documented.
  - The current segment parser cannot handle `contig1|03-M|...` even though that header contains usable segment information.

## 4. Primer database observations

- Current primer database structure:
  - A single JSON file: `/home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json`.
  - Top-level keys are organism/group names: `SARS-CoV-2`, `Influenza-A`, `Influenza-B`, `RSV-A`, `RSV-B`.
  - Each value is a simple object mapping primer name to primer sequence.
- File format:
  - Valid JSON.
  - No manifest, schema, version, provenance, source references, release date, or checksum metadata.
- Primer counts:
  - `SARS-CoV-2`: 3 primers/probe entries.
  - `Influenza-A`: 11 primers/probe entries.
  - `Influenza-B`: 7 primers/probe entries.
  - `RSV-A`: 3 primers/probe entries.
  - `RSV-B`: 3 primers/probe entries.
- Fields represented clearly:
  - Primer name and sequence are represented.
- Fields not represented explicitly:
  - Primer role/type, pool, gene/segment, organism/subtype, scheme name, scheme version, target coordinates, reference accession, citation/source, intended assay, strand/orientation, status/deprecation, date added, and notes.
- Inferred metadata currently depends on naming conventions:
  - Segment is inferred by checking primer names for `HA` or `M`.
  - Influenza A H1/H3 subsets are inferred by checking primer names for `H1` or `H3`, plus names containing `_M`.
- Database quality issues found:
  - NS segment exists in `Influenza-B` primer names but current code does not filter NS records explicitly.
  - Probe entries are stored in the same shape as primers, with no explicit `type`.
  - Mixed naming conventions: `triplex_*`, `H3_*`, `FluSw-H1-*`, `BHA*`, `probe-*`, `RSVQA*`.
  - No uniqueness checks beyond JSON object keys.
  - No validation for allowed IUPAC bases, duplicate sequences, expected primer/probe roles, or required metadata.
- Adding a new scheme today requires manual JSON edits and code-compatible naming conventions. That is fragile because a biologically meaningful field must be encoded in the primer name.
- Old schemes cannot be retained cleanly for reproducibility unless top-level keys or primer names are manually duplicated.

## 5. Primer database improvement options

### Option A: Add a JSON Schema and read-only validator for the current legacy format

- Current problem: The current database can be syntactically valid JSON while still having invalid sequences, duplicate/confusing names, unsupported organism keys, or metadata encoded only in names.
- Proposed structure or workflow: Keep accepting the existing `{organism: {primer_name: sequence}}` format, but add a validator that checks top-level structure, unique names within each organism, allowed IUPAC bases, non-empty sequences, and inferred metadata warnings.
- Expected benefit: Highest immediate safety with minimal compatibility risk.
- Backward compatibility risk: Low. Does not require changing the current database file.
- Suggested implementation outline: Add a `validate_primer_db()` function and a CLI mode such as `--validate-primers`; optionally add `scripts/validate_primer_db.py` later.
- Suggested validation tests: Valid legacy JSON passes; malformed JSON fails; empty sequence fails; sequence with invalid characters fails; influenza NS primer emits a warning until explicit metadata exists.

### Option B: Introduce a normalized v1 schema while preserving legacy loading

- Current problem: The legacy format cannot represent segment, subtype, primer role, pool, references, or scheme version.
- Proposed structure or workflow: Support a new versioned object shape:
  - `schema_version`
  - `database_version`
  - `schemes[]`
  - scheme fields: `scheme_id`, `display_name`, `organism`, `subtype`, `version`, `source`, `references`, `primers[]`
  - primer fields: `id`, `name`, `sequence`, `role`, `segment`, `gene`, `pool`, `strand`, `reference_start`, `reference_end`, `notes`
- Expected benefit: Makes primer schemes self-describing, easier to validate, and safer to update.
- Backward compatibility risk: Medium if the code switches too quickly. Low if the loader supports both legacy and normalized formats.
- Suggested implementation outline: Add dataclasses or typed dictionaries internally; write a loader that converts legacy JSON into the same internal representation as normalized v1.
- Suggested validation tests: Legacy and normalized fixtures load to equivalent internal primer records; explicit segment metadata controls filtering; missing required v1 fields fail with clear messages.

### Option C: Split database by organism/scheme/version with a manifest

- Current problem: A single file does not scale well for multiple organisms, assay schemes, and retained historical versions.
- Proposed structure or workflow:
  - `primer_db/manifest.json` lists available schemes and current defaults.
  - `primer_db/schemes/<organism>/<scheme_id>/<version>/primers.json` stores versioned scheme data.
  - Historical versions remain read-only and selectable.
- Expected benefit: Reproducibility, clearer ownership, easier review of updates, and safer addition of new schemes.
- Backward compatibility risk: Medium. Requires new path resolution and migration plan.
- Suggested implementation outline: Start with manifest support in repo tests and keep `--primers` accepting a direct JSON file. Add `--scheme` and `--scheme-version` later.
- Suggested validation tests: Manifest references existing files; default scheme resolves; old version remains selectable; missing file or duplicate default fails.

### Option D: Add a conversion script from legacy JSON to normalized v1

- Current problem: Manually rewriting existing primers into a richer format invites transcription errors.
- Proposed structure or workflow: A conversion script reads `fhi_primers.json`, infers what it can from names, marks uncertain fields as `null` or `inferred: true`, and writes a normalized draft for human review.
- Expected benefit: Reduces manual update risk and gives a repeatable migration path.
- Backward compatibility risk: Low if generated output is not used automatically until reviewed.
- Suggested implementation outline: Add `scripts/convert_legacy_primers.py --input <legacy.json> --output <normalized.json>`.
- Suggested validation tests: Conversion preserves all names and sequences; inferred segments for HA/M/NS match expected examples; uncertain fields are explicit.

### Option E: Add database provenance and release workflow

- Current problem: There is no way to tell which primer database version produced a report.
- Proposed structure or workflow: Include database version and scheme metadata in the loaded primer records, print it at runtime, and add database columns to the CSV such as `Primer_Database_Version`, `Primer_Scheme`, and `Primer_Scheme_Version`.
- Expected benefit: Improves reproducibility and auditability with real data.
- Backward compatibility risk: Medium because CSV columns change. This can be mitigated by appending columns rather than renaming existing ones.
- Suggested implementation outline: First load metadata internally; then append optional metadata columns once tests cover output headers.
- Suggested validation tests: Reports include metadata for normalized DB; legacy DB reports `legacy` or empty metadata consistently.

## 6. Problems and improvement candidates

### README/script name mismatch

- Type: documentation / usability
- Evidence from code, data, or primer_db: `README.md` line 5 and `primer_checker.py` lines 3 and 20 refer to `primer_blast.py`, but the repository entry point is `primer_checker.py`.
- Expected benefit: Reduces confusion when running or validating the tool.
- Risk level: low
- Suggested implementation outline: Update README and module docstring to consistently use `primer_checker.py`.
- Suggested tests: Run `python3 primer_checker.py --help` and ensure README examples match the CLI.

### Missing dependency and test setup

- Type: testing / documentation
- Evidence from code, data, or primer_db: No tests or project metadata were found. `pytest` is not installed in the current environment.
- Expected benefit: Makes behavior testable and gives future contributors a clear setup path.
- Risk level: low
- Suggested implementation outline: Add focused tests and minimal test dependency documentation. Avoid requiring BLAST for unit tests.
- Suggested tests: `pytest` unit tests for pure functions; `pytest --collect-only` should discover tests once pytest is installed.

### Missing BLAST availability preflight

- Type: validation / usability
- Evidence from code, data, or primer_db: `blastn -version` failed with command not found; current code discovers this only after creating a temporary query file and entering analysis.
- Expected benefit: Clear, early error before processing many inputs.
- Risk level: low
- Suggested implementation outline: Add a startup check for `blastn` using `shutil.which()` and report installation guidance before processing files.
- Suggested tests: Mock missing and present `blastn`; ensure missing BLAST exits before FASTA processing.

### Temporary query file can leak when BLAST is missing

- Type: bug
- Evidence from code, data, or primer_db: `run_blastn()` creates a temp file at lines 175-177, catches `FileNotFoundError` at lines 198-199, and exits before `os.remove(query_filename)` at line 201.
- Expected benefit: Avoids leftover temp files on error paths.
- Risk level: low
- Suggested implementation outline: Use `try/finally` around temporary file cleanup, or use `TemporaryDirectory`.
- Suggested tests: Mock `subprocess.run` raising `FileNotFoundError` and assert cleanup is attempted.

### Global primer state

- Type: refactor / testing
- Evidence from code, data, or primer_db: `process_fasta_file()` reads global `VIRUS_PRIMERS` at lines 271-272; `main()` mutates it at lines 410-424.
- Expected benefit: Makes behavior easier to test, reuse, and reason about.
- Risk level: low
- Suggested implementation outline: Pass selected primer records into processing functions rather than reading global state.
- Suggested tests: Unit test `process_fasta_file()` with injected primers and mocked `run_blastn()`.

### Fragile influenza segment inference from primer names

- Type: validation / database
- Evidence from code, data, or primer_db: `process_fasta_file()` infers segment by checking if `"HA"` or `"M"` appears in the primer name at lines 284-291. The primer database contains `triplex_InfB_F_NS`, `triplex_InfB_R_NS`, and `triplex_InfB_Probe_NS`, but NS is not recognized.
- Expected benefit: Correct filtering for NS and future segments, fewer false matches from incidental letters.
- Risk level: medium
- Suggested implementation outline: Add explicit segment metadata via normalized primer records; as an interim fix, recognize `_NS` with boundary-aware parsing and warn when segment cannot be inferred.
- Suggested tests: NS primers only process NS subjects; HA/M behavior remains unchanged; primer names with incidental letters do not trigger segment filtering.

### Fragile FASTA header segment parsing

- Type: validation / usability
- Evidence from code, data, or primer_db: `get_segment()` assumes segment is in the first pipe-delimited token, split by hyphen, at lines 143-147. Real data includes `contig1|03-M|INFL16-2025`.
- Expected benefit: Clearer handling of realistic FASTA headers and fewer silent skips.
- Risk level: medium
- Suggested implementation outline: Define supported segment extraction patterns, return parse status, and warn/report unparseable segmented records.
- Suggested tests: Expected headers, contig-prefixed headers, strain-name headers, missing pipes, extra hyphens, lowercase segments, and non-segmented RSV/SARS-CoV-2 headers.

### No FASTA validation or summary

- Type: validation / usability
- Evidence from code, data, or primer_db: One invalid alphabet example was found in `INFA_test.fasta`; current code only extracts headers and delegates sequence content to BLAST.
- Expected benefit: Faster diagnosis of malformed inputs before BLAST failures or confusing no-hit output.
- Risk level: low
- Suggested implementation outline: Add an optional FASTA validation function that checks headers, empty records, duplicate IDs, allowed IUPAC alphabet, sequence counts, and segment summary.
- Suggested tests: Invalid characters fail or warn with file/header/line; duplicate IDs warn; empty sequence fails.

### Output CSV mixes numeric and string values

- Type: usability / validation
- Evidence from code, data, or primer_db: No-hit rows store `"No hit"` in `Percent_Identity` at line 323 and blanks in numeric columns.
- Expected benefit: Easier downstream analysis in PowerBI/R/Python and clearer status handling.
- Risk level: medium
- Suggested implementation outline: Add a separate `Hit_Status` column and use empty numeric fields for numeric columns; preserve current columns initially if backward compatibility matters.
- Suggested tests: No-hit row has `Hit_Status=no_hit`; hit row has numeric percent identity; CSV header is stable.

### No visual run report for investigation

- Type: feature / usability
- Evidence from code, data, or primer_db: The tool currently writes only a CSV report. The repository goal includes making results easier to validate with real data, and the realistic data spans multiple organisms, segments, primers/probes, and many samples.
- Expected benefit: Faster review after each run, easier primer-level investigation, and clearer communication of which primers/organisms/samples have poor matches or no hits.
- Risk level: medium
- Suggested implementation outline: Add an optional static HTML report generated after analysis, preferably with embedded JSON data and local CSS/JS so it can be opened directly in a browser. Include an overview dashboard by organism/virus, FASTA file, primer, hit status, mismatch counts, percent identity, and segment. Add controls to filter/select organisms, primers, FASTA files, hit statuses, and mismatch thresholds. Add per-primer visual panels showing summary counts, percent identity/mismatch distributions, mismatch-position maps, no-hit counts, risk badges, and a sortable sample table. Keep CSV as the canonical machine-readable output.
- Suggested tests: Generate HTML from small fixture results; assert expected sections and embedded data exist; test escaping of primer names/sample IDs; test filtering data model includes organism, primer, segment, hit status, mismatches, mismatch positions, and percent identity.

### Visual report lacks deeper primer risk and trend views

- Type: feature / usability
- Evidence from code, data, or primer_db: The first HTML report provided filters and per-primer tables, but did not highlight primer risk, mismatch locations, mismatch-count ratios, or previous-report context.
- Expected benefit: Makes risky primers easier to triage, helps distinguish terminal-end mismatches from middle-primer mismatches, and gives analysts a lightweight way to compare against previous runs.
- Risk level: medium
- Suggested implementation outline: Add risk thresholds in one JS config object; calculate primer stats from currently filtered rows; add Risk column and badges to the overview table; add overall and per-primer mismatch-count distribution plots; add SVG mismatch-position maps with 5' and 3' terminal highlighting; add optional `--previous-report-csv` attachments embedded as `previous-report-data`.
- Suggested tests: Fixture HTML contains risk thresholds, Risk column, mismatch distribution sections, sequence map labels, previous-report data, and safe JSON serialization; previous CSV loader handles valid, missing, and malformed files without crashing.

### Hard-coded BLAST parameters

- Type: feature / maintainability
- Evidence from code, data, or primer_db: `run_blastn()` hard-codes reward, penalty, word size, dust, and outfmt at lines 179-188.
- Expected benefit: Easier validation and reproducible tuning for different assays.
- Risk level: medium
- Suggested implementation outline: Move parameters to constants first, then optional CLI/config fields later.
- Suggested tests: Command construction tests verify defaults and any overrides.

### Silent skipping and weak exit semantics

- Type: validation / usability
- Evidence from code, data, or primer_db: Missing FASTA files are printed and skipped at lines 431-434; missing primer type returns no results at lines 272-275; no results only prints "No results to write" at lines 340-342.
- Expected benefit: Better automation behavior and less chance of empty successful runs.
- Risk level: medium
- Suggested implementation outline: Return non-zero exit status for invalid inputs or no analyzable data; collect warnings and errors.
- Suggested tests: Missing all FASTA files exits non-zero; unknown virus exits non-zero; some missing FASTA files produce warning but process valid files according to documented policy.

### Primer database lacks schema, metadata, and versioning

- Type: database
- Evidence from code, data, or primer_db: `fhi_primers.json` is a valid but minimal mapping of organism to primer name to sequence. It has no scheme metadata, references, pools, roles, segments, versions, or provenance.
- Expected benefit: Safer updates, reproducible reports, and support for multiple schemes.
- Risk level: medium
- Suggested implementation outline: Add legacy validation first, then normalized v1 loader with explicit metadata and compatibility conversion.
- Suggested tests: Legacy DB validates; normalized DB validates; version/provenance fields are exposed to reporting without breaking legacy input.

### Primer database update workflow is manual and fragile

- Type: database / usability
- Evidence from code, data, or primer_db: Updating requires hand-editing a single JSON file and preserving implicit naming conventions used by code.
- Expected benefit: Safer addition of new schemes and easier review of changes.
- Risk level: medium
- Suggested implementation outline: Add conversion and validation scripts; eventually add manifest-based scheme discovery.
- Suggested tests: New scheme fixture validates; duplicate primer IDs fail; old scheme remains selectable.

### Primer database should move to normalized versioned JSON

- Type: database / documentation
- Evidence from code, data, or primer_db: The current legacy JSON does not store explicit segment, role, source, scheme, version, or status fields.
- Expected benefit: Easier updates, safer validation, explicit primer metadata, and retained old versions for reproducibility.
- Risk level: medium
- Suggested implementation outline: Keep legacy loading, add support for a normalized `schema_version`/`database_version`/`schemes[]` JSON format, add a validator, and add a converter from current legacy JSON or CSV into normalized scheme files.
- Suggested tests: Normalized fixture validates and loads; missing required fields fail clearly; legacy and normalized records produce equivalent analysis behavior for the current FHI primers.

## 7. Recommended implementation order

General code improvements:

1. Add unit tests for pure functions: IUPAC matching, mismatch positions, subject ID parsing, segment parsing, influenza subtype selection.
2. Add BLAST preflight and temp-file cleanup.
3. Remove global primer state by passing selected primer records/functions explicitly.
4. Improve FASTA header parsing and add warnings/summaries for skipped segmented records.
5. Add optional FASTA validation.
6. Improve output status semantics with a separate hit status field.
7. Add an optional static HTML visual report from the same structured results used for CSV output.
8. Move hard-coded BLAST settings to named defaults, then optional CLI/config later.

Primer database improvements:

1. Add a legacy primer database validator for the current JSON shape.
2. Add an internal normalized primer record model while preserving legacy JSON loading.
3. Add explicit segment support, including NS, using metadata where available and conservative inference for legacy DB.
4. Add a normalized v1 schema and fixtures in the repository.
5. Add a conversion script from legacy JSON to normalized v1 draft.
6. Add manifest/scheme/version layout once the loader and validator are tested.
7. Add database version/scheme metadata columns to CSV reports.

## 8. First implementation batch

Recommended first batch:

1. Add tests for existing pure behavior and edge cases discovered here.
2. Add a legacy primer DB validation function and CLI path that validates `fhi_primers.json` without modifying it.
3. Add a small internal `PrimerRecord` representation and a backward-compatible loader that converts the legacy JSON shape into records.
4. Fix BLAST preflight/temp-file cleanup.
5. Add explicit legacy segment inference support for HA, M, and NS with warnings for unparseable influenza headers.
6. Add a first static HTML visual report output that can be opened after each run and filtered by organism/virus, primer, FASTA file, hit status, and mismatch threshold.
7. Update README naming/setup notes.

Done-when criteria for this batch:

- Existing CLI remains backward compatible for `--primers <legacy.json>`.
- Unit tests cover core parsing/mismatch/database validation behavior without requiring BLAST.
- Primer DB validation reports clear errors/warnings and passes the existing `fhi_primers.json` except for documented warnings.
- NS influenza B primers no longer run against every segment.
- Missing `blastn` fails early with a clear message and no temp-file leak.
- A generated HTML report can be opened locally and shows overview plus per-primer investigation views using the same result rows as the CSV.

Implementation status:

- Implemented in `primer_checker.py`: `PrimerRecord`, legacy primer DB validation, backward-compatible record loading, HA/M/NS segment inference, improved FASTA segment parsing for standard and contig-prefixed influenza headers, BLAST preflight, temp-file cleanup, `--validate-primers`, and optional `--html-report`.
- Implemented in `primer_checker.py` second report pass: primer risk categories, sortable Risk overview column, mismatch-count ratio plots, per-primer mismatch-position SVG maps, 5'/3' terminal mismatch risk contribution, and optional repeated `--previous-report-csv` attachment support.
- Implemented in `primer_checker.py` third report pass: default highest-mismatch/no-hit ordering in per-primer sample tables, clickable per-primer table column sorting, and 10-row collapsed tables with show-all/show-top-10 toggles. The top-level "Primers to Review" and "Mismatch Count Distribution" modules were removed after review; the Primer Overview risk column and per-primer mismatch distribution plots remain.
- Implemented in `tests/test_primer_checker.py`: focused tests for IUPAC matching, mismatch positions, primer DB validation/loading, influenza subtype and NS segment selection, FASTA segment parsing, HTML report generation/escaping, BLAST missing-path cleanup, and validation-only CLI behavior.
- Implemented in `HTML_REPORT_IMPROVEMENTS.md`: short user/developer note describing report changes, risk categories, previous CSV attachment, assumptions, and next improvements.
- Implemented normalized primer DB support: `primer_checker.py` now validates and loads both legacy and normalized JSON formats; `scripts/convert_legacy_primers.py` converts old JSON to normalized schema v1; `primer_db/fhi_primers.normalized.json` is the converted FHI database copy in this repo; `PRIMER_DATABASE_MAINTENANCE.md` documents the workflow.
- Implemented filename-based batch routing: H1/H3 FASTA filenames now route to matching subtype-tagged Influenza-A primers plus untagged Influenza-A primers, excluding the opposite subtype; `scripts/run_primer_checker_batch.py` scans an input folder and runs each FASTA against the inferred SARS-CoV-2, Influenza-A/H1/H3/B, RSV-A, or RSV-B primer set.
- Implemented optional metadata attachment: `--metadata-csv` loads SampleID/Sample_Date plus generic or assay-specific Ct metadata, attaches matched metadata to CSV/HTML rows with exact-token FASTA-header matching, and the HTML report now includes a per-primer timeline for average mismatches and selected Ct over time. Supported assay-style columns include `CT_H3`, `CT_H1`, `Triplex-InfA_CT`, `Triplex-InfB_CT`, `Triplex-SC2_CT`, `CT_RSVA`, `CT_RSVB`, `CT_INFA`, and `CT_INFB`; result rows include `Ct_Source` so the selected assay column is visible. `fixtures/fake_metadata.csv` provides a small realistic test metadata file based on sample IDs from the external FASTA folder.
- Implemented R metadata export helpers: `scripts/primer_checker_metadata_helpers.R` defines the shared primer-checker metadata schema and builder functions; `scripts/influenza_gisaid.R`, `scripts/rsv_gisaid.R`, and `scripts/sc2_gisaid.R` now write primer-checker metadata CSVs only. They no longer create GISAID CSV/Excel/FASTA outputs and no longer filter by SID/run ID.
- Updated influenza metadata mapping for real `fludb` fields: `pcr_h1_ct` -> `CT_H1`, `pcr_h3_ct` -> `CT_H3`, `pcr_bvic_ct`/`pcr_byam_ct` -> `CT_INFB`, and `prove_tatt` -> `Sample_Date`. Influenza triplex Ct output columns remain blank until triplex Ct fields exist in the source data.
- Enhanced HTML report inspection: plot y-axes were added to the primer mismatch map and timeline; mismatch-position graphics now surface nucleotide-change details; per-primer detail tables keep no-hit rows at the bottom unless sorting by hit status; clicking a sample ID opens an offline alignment popup using embedded BLAST query/subject alignment strings.
- Implemented in `README.md`: script-name cleanup, validation mode documentation, supported influenza header patterns, pytest notes, and optional HTML report usage.
- Verification completed locally: syntax checks pass; `python3 primer_checker.py --help` works; the external `fhi_primers.json` validates without fatal errors; missing `blastn` fails early with a clear message; a fixture HTML report was generated under `/tmp`.
- Verification completed for enhanced report: pytest suite passes locally; generated enhanced fixture report contains risk thresholds, mismatch distribution sections, sequence-map code, and embedded previous report data; extracted generated JavaScript passes `node --check`.
- Verification not completed locally: full pytest execution, because `python3 -m pytest -q` fails in this environment with `No module named pytest`.

Remaining follow-up work:

- Install pytest in the execution environment and run the full test suite.
- Install BLAST+ and run a real end-to-end analysis against the external FASTA files.
- Consider whether the new CSV columns (`Primer_Segment`, `Subject_Segment`, `Hit_Status`) need downstream PowerBI adjustments.

## 9. Commands discovered

- `python3 primer_checker.py --help`
- `python3 primer_checker.py --primers /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json --virus influenza --flu-type A --fasta /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/data/INFA_test.fasta --output primer_report.csv`
- `python3 primer_checker.py --primers /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json --virus influenza --flu-type B --fasta /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/data/INFB_test.fasta --output primer_report.csv`
- `python3 primer_checker.py --primers /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json --virus SARS-CoV-2 --fasta /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/data/SC2.fasta --output primer_report.csv`
- Future report command shape to consider: `python3 primer_checker.py --primers <primers.json> --virus <virus> --fasta <files...> --output primer_report.csv --html-report primer_report.html`
- `python3 scripts/run_primer_checker_batch.py --input-folder <folder> --primers primer_db/fhi_primers.normalized.json --dry-run`
- `python3 scripts/run_primer_checker_batch.py --input-folder <folder> --primers primer_db/fhi_primers.normalized.json --output batch_primer_report.csv --html-report batch_primer_report.html`
- `python3 scripts/run_primer_checker_batch.py --input-folder <folder> --primers primer_db/fhi_primers.normalized.json --metadata-csv fixtures/fake_metadata.csv --output batch_primer_report.csv --html-report batch_primer_report.html`
- `jq empty /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json`
- `jq -r 'to_entries[] | [.key, (.value|type), (.value|length)] | @tsv' /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db/fhi_primers.json`
- `grep -Hc '^>' /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/data/*.fasta`
- `blastn -version` currently fails because `blastn` is not installed.
- `pytest --collect-only -q -p no:cacheprovider` currently fails because `pytest` is not installed.

## 10. Open questions

None that block the recommended first implementation batch.

Questions for later database design:

- Which primer scheme identifiers and version names should be authoritative for the FHI primer sets?
- Should probes be analyzed exactly like primers, or should probe rows have separate thresholds/semantics?
- Should CSV output append database metadata columns immediately, or wait until downstream PowerBI expectations are confirmed?

## 11. Suggested next /goal prompt

```text
/goal Implement the first low-risk reliability and primer database safety batch from CODEX_IMPROVEMENT_PLAN.md.

Scope:
- Keep the existing CLI backward compatible for legacy primer JSON files shaped as {organism: {primer_name: sequence}}.
- Do not modify files under /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/data or /home/rasmuskopperud.riis/Coding/run_folder/data/primer_checker/primer_db.
- Use the external fhi_primers.json and FASTA files only for read-only validation/smoke checks.
- Add focused tests that do not require BLAST.
- Add an optional local visual report, but keep CSV as the canonical machine-readable output.

Implement:
1. Add a small internal PrimerRecord representation and a backward-compatible primer library loader that converts the current legacy JSON into records.
2. Add primer database validation for legacy JSON: top-level/object shape, non-empty primer names, non-empty sequences, allowed IUPAC DNA/RNA bases, duplicate record detection where applicable, and useful warnings for inferred metadata.
3. Add explicit segment inference for legacy influenza primer names covering HA, M, and NS with boundary-aware matching, and use it so Influenza-B NS primers only apply to NS subjects.
4. Improve FASTA segment parsing enough to handle both existing `03-M|...` headers and discovered `contig1|03-M|...` style headers; report/warn on unparseable influenza headers instead of silently ignoring them.
5. Add BLAST availability preflight and ensure temporary query files are cleaned up on missing-BLAST/error paths.
6. Add an optional static HTML report output, for example `--html-report primer_report.html`, generated from the same structured result rows as the CSV. The report should open directly in a browser without a server or external network dependencies. It should include:
   - Overview cards/tables by organism or virus, FASTA file, primer, segment, hit status, no-hit count, mismatch count, and percent identity.
   - Per-primer investigation panels with mismatch/identity summaries and sortable sample rows.
   - Client-side controls to choose/filter organisms, primers, FASTA files, hit status, and mismatch thresholds.
   - Safe escaping/serialization for sample IDs, primer names, and sequences.
7. Add unit tests for IUPAC matching, mismatch positions, primer DB validation/loading, influenza subtype/segment selection, FASTA segment parsing, and HTML report generation from fixture rows.
8. Update README/docstring naming and setup notes, including the real script name, BLAST/pytest expectations, and the optional visual report.

Done when:
- `python3 primer_checker.py --help` still works.
- The test suite passes in an environment with pytest installed and does not require BLAST for unit tests.
- The existing external `fhi_primers.json` validates without fatal errors.
- Influenza-B NS primers are filtered to NS records in tests.
- Missing `blastn` produces an early clear error and does not leak temporary query files.
- `--html-report <path>` produces a local HTML file from fixture or mocked result rows, and tests confirm it contains the overview data, per-primer data, and filterable report data model.
- CODEX_IMPROVEMENT_PLAN.md is updated with what was implemented and any remaining follow-up work.
```
