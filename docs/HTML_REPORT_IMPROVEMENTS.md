# HTML Report Improvements

## What changed

- Added primer risk badges and row highlighting in the Primer Overview table.
- Added an overall mismatch-count distribution plot that updates with filters.
- Added per-primer mismatch-count distribution plots.
- Added per-primer mismatch-location SVG diagrams across the primer sequence.
- Highlighted terminal 5' and 3' primer regions in the mismatch-position risk logic and sequence map.
- Added optional attachment of previous CSV reports through `--previous-report-csv`.
- Made per-primer sample tables sortable by clicking table headers.
- Defaulted per-primer sample tables to show highest mismatch/no-hit rows first.
- Collapsed per-primer sample tables to 10 rows by default with a show-all/show-top-10 toggle.
- Removed the top-level "Primers to Review" and "Mismatch Count Distribution" modules after review; risk remains in Primer Overview and mismatch distributions remain inside primer panels.
- Kept the report static, self-contained, and offline-friendly with embedded JSON, CSS, SVG, and vanilla JavaScript.

## Risk categories

Risk is calculated from the currently filtered rows per primer. No-hit rows are counted separately from numeric mismatch averages.

Default thresholds are defined in the generated report JavaScript as `RISK_THRESHOLDS`:

- `Critical`: no-hit rate >= 20%, average mismatches >= 3, 3+ mismatch rate >= 25%, or terminal mismatch share >= 50%.
- `High`: no-hit rate >= 5%, average mismatches >= 2, 2+ mismatch rate >= 30%, or terminal mismatch share >= 25%.
- `Watch`: average mismatches >= 1, max mismatches >= 2, or any no-hit rows.
- `Low`: none of the above.

Terminal mismatch share means the proportion of observed mismatch-position events that fall within the first or last five primer bases. This is included because mismatches near the 5' and especially 3' primer ends can matter more than mismatches in the middle.

## Attaching previous reports

Attach one or more previous CSV reports when generating the HTML:

```bash
python3 primer_checker.py \
  --primers primers.json \
  --virus influenza \
  --flu-type B \
  --fasta current.fasta \
  --output primer_report.csv \
  --html-report primer_report.html \
  --previous-report-csv old_report_1.csv \
  --previous-report-csv old_report_2.csv
```

Attached reports are parsed by the CLI and embedded into the HTML in `previous-report-data`. If no previous CSV is provided, the report still works and shows an empty-state message. Invalid or unreadable CSVs are embedded with warnings instead of crashing report generation.

## Assumptions

- CSV output remains the canonical machine-readable output.
- The HTML report is for local investigation and does not load external scripts or send data anywhere.
- Existing result columns such as `Primer_Name`, `Mismatches`, `Mismatch_Positions`, `Hit_Status`, and `Percent_Identity` are enough for the first risk and plotting pass.
- Mismatch-position plots use hit rows only; no-hit rows are represented in mismatch-count distributions and risk calculations.

## Suggested next improvements

- Add user-configurable risk thresholds through a CLI option or small JSON config.
- Add optional side-by-side comparison plots for current versus previous reports.
- Add a primer-end weighted risk score once domain-specific thresholds are agreed.
- Add browser-level regression tests for sorting, expanding, and filtering interactions.
- Add an end-to-end browser smoke test with Playwright once the project has a test environment.
