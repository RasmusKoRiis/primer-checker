#!/usr/bin/env python3
"""HTML report generation for primer checker result rows."""

import csv
import html
import json
import os
import re
import sys
from pathlib import Path

from primer_analysis import CSV_FIELDNAMES, infer_primer_role

REPORT_TEXT_DIRECTORY = Path(__file__).resolve().parent / "report_text"
REPORT_TEXT_FILES = {
    "en": "english.json",
    "no": "norwegian.json",
}
REPORT_TEXT_PLACEHOLDER_PATTERN = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}")

METADATA_REPORT_FIELDNAMES = {
    "Metadata_Sample_ID",
    "Sample_Date",
    "Ct_Value",
    "Ct_Source",
}
REPORT_FIELDNAMES = [field for field in CSV_FIELDNAMES if field not in METADATA_REPORT_FIELDNAMES] + ["Primer_Role"]


def load_report_translations(
    text_directory: str | os.PathLike[str] | None = None,
) -> dict[str, dict[str, str]]:
    """Load and validate the editable report text files."""
    directory = Path(text_directory) if text_directory is not None else REPORT_TEXT_DIRECTORY
    translations: dict[str, dict[str, str]] = {}

    for language, filename in REPORT_TEXT_FILES.items():
        path = directory / filename
        try:
            with path.open(encoding="utf-8") as text_file:
                values = json.load(text_file)
        except FileNotFoundError as error:
            raise FileNotFoundError(f"Report text file not found: {path}") from error
        except json.JSONDecodeError as error:
            raise ValueError(
                f"Invalid JSON in report text file {path} at line {error.lineno}, "
                f"column {error.colno}: {error.msg}"
            ) from error

        if not isinstance(values, dict):
            raise ValueError(f"Report text file {path} must contain one JSON object.")

        invalid_keys = [key for key in values if not isinstance(key, str)]
        invalid_values = [key for key, value in values.items() if not isinstance(value, str)]
        if invalid_keys:
            raise ValueError(f"Report text file {path} contains a non-text key.")
        if invalid_values:
            joined_keys = ", ".join(sorted(invalid_values))
            raise ValueError(f"Report text values must be text in {path}: {joined_keys}")

        translations[language] = values

    english_keys = set(translations["en"])
    for language, values in translations.items():
        language_keys = set(values)
        missing_keys = sorted(english_keys - language_keys)
        extra_keys = sorted(language_keys - english_keys)
        if missing_keys or extra_keys:
            details = []
            if missing_keys:
                details.append("missing: " + ", ".join(missing_keys))
            if extra_keys:
                details.append("extra: " + ", ".join(extra_keys))
            raise ValueError(
                f"Report text keys in {REPORT_TEXT_FILES[language]} do not match english.json "
                f"({'; '.join(details)})."
            )

        for key in english_keys:
            english_placeholders = set(REPORT_TEXT_PLACEHOLDER_PATTERN.findall(translations["en"][key]))
            translated_placeholders = set(REPORT_TEXT_PLACEHOLDER_PATTERN.findall(values[key]))
            if translated_placeholders != english_placeholders:
                raise ValueError(
                    f"Placeholders for report text key {key!r} do not match between "
                    f"english.json and {REPORT_TEXT_FILES[language]}."
                )

    return translations


def _json_for_inline_script(value) -> str:
    """Serialize JSON without allowing values to terminate an inline script element."""
    return (
        json.dumps(value, ensure_ascii=True)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )


def prepare_report_row(row: dict) -> dict:
    prepared = {field: row.get(field, "") for field in REPORT_FIELDNAMES}
    prepared["Primer_Role"] = prepared["Primer_Role"] or infer_primer_role(prepared["Primer_Name"])
    return prepared


def safe_number(value):
    if value in ("", None, "No hit"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_previous_report_csv(path: str) -> dict:
    """Read a previous CSV report so its rows can be embedded in the static HTML report."""
    report = {
        "name": os.path.basename(path),
        "path": path,
        "rows": [],
        "warnings": [],
    }
    try:
        with open(path, newline="", encoding="utf-8") as csvfile:
            reader = csv.DictReader(csvfile)
            if not reader.fieldnames:
                report["warnings"].append("CSV has no header row.")
                return report
            missing_fields = [field for field in ("Primer_Name", "Mismatches", "Hit_Status") if field not in reader.fieldnames]
            if missing_fields:
                report["warnings"].append(f"Missing expected column(s): {', '.join(missing_fields)}.")
            for row in reader:
                report["rows"].append(prepare_report_row(row))
    except Exception as e:
        report["warnings"].append(f"Could not read CSV: {e}")
    return report


def load_previous_report_csvs(paths: list[str] | None) -> list[dict]:
    return [load_previous_report_csv(path) for path in paths or []]


def build_html_report(
    results: list[dict],
    title: str | None = None,
    previous_reports: list[dict] | None = None,
) -> str:
    """Build a self-contained HTML report with embedded result data."""
    report_rows = []
    for row in results:
        report_rows.append(prepare_report_row(row))

    previous_reports = [
        {
            **report,
            "rows": [prepare_report_row(row) for row in report.get("rows", [])],
        }
        for report in (previous_reports or [])
    ]
    translations = load_report_translations()
    data_json = _json_for_inline_script(report_rows)
    translations_json = _json_for_inline_script(translations)
    fields_json = json.dumps(REPORT_FIELDNAMES)
    escaped_title = html.escape(title or translations["en"]["report_title"])
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__REPORT_TITLE__</title>
  <style>
    :root {
      --bg: #f7f8fb;
      --panel: #ffffff;
      --ink: #17202a;
      --muted: #5d6878;
      --line: #d8dee9;
      --accent: #1f7a8c;
      --accent-weak: #e6f4f6;
      --bad: #b42318;
      --warn: #b7791f;
      --good: #287d3c;
      --critical-bg: #fde2e1;
      --high-bg: #fff0e6;
      --watch-bg: #fff8db;
      --low-bg: #e7f5eb;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: var(--bg);
    }
    header {
      padding: 24px 32px 16px;
      background: var(--panel);
      border-bottom: 1px solid var(--line);
    }
    h1 { margin: 0 0 6px; font-size: 28px; letter-spacing: 0; }
    h2 { margin: 24px 0 6px; font-size: 20px; letter-spacing: 0; }
    h3 { margin: 0 0 10px; font-size: 16px; letter-spacing: 0; }
    main { padding: 20px 32px 36px; }
    .meta, .section-note { color: var(--muted); margin: 0 0 12px; }
    .report-note {
      display: grid;
      gap: 8px;
      margin: 8px 0 14px;
      padding: 12px 14px;
      background: #f8fafc;
      border: 1px solid var(--line);
      border-left: 4px solid var(--accent);
      border-radius: 8px;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.55;
    }
    .report-note-compact { margin-top: 10px; }
    .note-row {
      display: grid;
      grid-template-columns: minmax(110px, 170px) 1fr;
      gap: 10px;
      align-items: start;
    }
    .note-row strong {
      color: var(--ink);
      font-weight: 700;
    }
    .report-note .table-toggle {
      justify-self: start;
      margin-top: 2px;
    }
    .filters {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      padding: 16px;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
    }
    label { display: grid; gap: 6px; color: var(--muted); font-size: 13px; }
    select, input {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 8px 10px;
      font: inherit;
      color: var(--ink);
      background: #fff;
    }
    .cards {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 12px;
      margin-top: 16px;
    }
    .card, .plot-card {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
    }
    .metric { font-size: 24px; font-weight: 700; margin-top: 4px; }
    .table-wrap {
      overflow-x: auto;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
    }
    table { width: 100%; border-collapse: collapse; min-width: 860px; }
    th, td {
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
      font-size: 13px;
    }
    th { background: #eef2f7; cursor: pointer; user-select: none; white-space: nowrap; }
    tr:last-child td { border-bottom: 0; }
    .badge, .risk-badge {
      display: inline-block;
      padding: 2px 8px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 600;
      white-space: nowrap;
    }
    .hit { color: var(--good); background: #e7f5eb; }
    .no_hit { color: var(--bad); background: #fde8e7; }
    .indeterminate, .risk-indeterminate { color: #475569; background: #e2e8f0; }
    .risk-low { color: var(--good); background: var(--low-bg); }
    .risk-watch { color: #8a5a00; background: var(--watch-bg); }
    .risk-high { color: #9a4b00; background: var(--high-bg); }
    .risk-critical { color: var(--bad); background: var(--critical-bg); }
    tr.risk-watch-row { background: #fffdf2; border-left: 4px solid #d69e2e; }
    tr.risk-high-row { background: #fff8f0; border-left: 4px solid #dd6b20; }
    tr.risk-critical-row { background: #fff5f5; border-left: 4px solid var(--bad); }
    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
      margin: 10px 0 12px;
      color: var(--muted);
      font-size: 13px;
    }
    .report-tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin: 0 0 16px;
    }
    .tab-button {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--panel);
      color: var(--ink);
      cursor: pointer;
      padding: 8px 12px;
      font: inherit;
      font-weight: 600;
    }
    .tab-button.active {
      background: var(--accent);
      border-color: var(--accent);
      color: #fff;
    }
    .tab-panel[hidden] { display: none; }
    .language-toggle {
      display: flex;
      justify-content: flex-end;
      gap: 8px;
      margin: 0 0 14px;
    }
    .language-button {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--panel);
      color: var(--ink);
      cursor: pointer;
      padding: 6px 10px;
      font: inherit;
      font-size: 13px;
      font-weight: 600;
    }
    .language-button.active {
      background: var(--accent-weak);
      border-color: var(--accent);
    }
    .primer-panel {
      margin-top: 12px;
      padding: 14px;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      scroll-margin-top: 18px;
    }
    .primer-panel.risk-high-panel { border-left: 5px solid #dd6b20; }
    .primer-panel.risk-critical-panel { border-left: 5px solid var(--bad); }
    .primer-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
      gap: 10px;
      margin: 10px 0;
      color: var(--muted);
    }
    .bar, .stacked-bar {
      height: 8px;
      background: #edf1f5;
      border-radius: 999px;
      overflow: hidden;
      margin-top: 6px;
    }
    .bar span { display: block; height: 100%; background: var(--accent); }
    .plot-card { margin: 14px 0; }
    .mismatch-distribution {
      display: grid;
      gap: 8px;
      margin-top: 8px;
    }
    .distribution-row {
      display: grid;
      grid-template-columns: 70px 1fr 110px;
      gap: 10px;
      align-items: center;
      font-size: 13px;
    }
    .distribution-track {
      height: 14px;
      background: #edf1f5;
      border-radius: 999px;
      overflow: hidden;
    }
    .distribution-fill {
      height: 100%;
      min-width: 2px;
      background: var(--accent);
    }
    .dist-0 { background: var(--good); }
    .dist-1 { background: #7aa6a9; }
    .dist-2 { background: var(--warn); }
    .dist-3, .dist-4plus { background: #dd6b20; }
    .dist-no_hit { background: var(--bad); }
    .dist-indeterminate { background: #94a3b8; }
    .primer-sequence-map {
      width: 100%;
      overflow-x: auto;
      padding-bottom: 4px;
    }
    .alignment-button {
      border: 0;
      background: transparent;
      color: var(--accent);
      cursor: pointer;
      font: inherit;
      padding: 0;
      text-decoration: underline;
      text-underline-offset: 2px;
    }
    .primer-link {
      color: var(--accent);
      font-weight: 700;
      text-decoration: underline;
      text-underline-offset: 2px;
    }
    .modal-backdrop {
      position: fixed;
      inset: 0;
      display: none;
      align-items: center;
      justify-content: center;
      padding: 20px;
      background: rgba(15, 23, 42, 0.48);
      z-index: 1000;
    }
    .modal-backdrop.open { display: flex; }
    .modal {
      width: min(980px, 96vw);
      max-height: 88vh;
      overflow: auto;
      background: var(--panel);
      border-radius: 8px;
      border: 1px solid var(--line);
      box-shadow: 0 18px 60px rgba(15, 23, 42, 0.28);
      padding: 18px;
    }
    .modal-header {
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: start;
      margin-bottom: 12px;
    }
    .modal-close {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      cursor: pointer;
      padding: 6px 10px;
      font: inherit;
    }
    .alignment-view {
      display: grid;
      gap: 6px;
      overflow-x: auto;
      padding: 12px;
      background: #f8fafc;
      border: 1px solid var(--line);
      border-radius: 6px;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      font-size: 13px;
      white-space: pre;
    }
    .alignment-row {
      display: grid;
      grid-template-columns: 72px max-content;
      gap: 10px;
    }
    .alignment-matchline { color: var(--muted); }
    .sequence-svg {
      max-width: 100%;
      min-width: 520px;
      height: auto;
      display: block;
    }
    .axis-label { fill: var(--muted); font-size: 11px; }
    .axis-title { fill: var(--muted); font-size: 12px; font-weight: 600; }
    .base-label { fill: var(--ink); font-size: 11px; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    .three-prime { fill: #b42318; font-weight: 700; }
    .five-prime { fill: #1f7a8c; font-weight: 700; }
    .review-reason {
      color: var(--muted);
      font-size: 13px;
    }
    .table-actions {
      display: flex;
      justify-content: space-between;
      gap: 10px;
      align-items: center;
      margin: 8px 0;
      color: var(--muted);
      font-size: 13px;
    }
    .table-toggle {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--accent-weak);
      color: var(--ink);
      cursor: pointer;
      padding: 6px 10px;
      font: inherit;
    }
    .primer-detail-table th {
      cursor: pointer;
    }
    .comparison-table td.delta-up { color: var(--bad); font-weight: 700; }
    .comparison-table td.delta-down { color: var(--good); font-weight: 700; }
    .amplicon-map {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(110px, 1fr));
      gap: 8px;
      margin: 12px 0;
    }
    .amplicon-tile {
      min-height: 78px;
      padding: 9px;
      border: 1px solid var(--line);
      border-left-width: 5px;
      border-radius: 7px;
      background: var(--panel);
      font-size: 12px;
    }
    .amplicon-tile strong { display: block; margin-bottom: 4px; }
    .amplicon-tile.risk-low { border-left-color: var(--good); }
    .amplicon-tile.risk-watch { border-left-color: #d69e2e; }
    .amplicon-tile.risk-high { border-left-color: #dd6b20; }
    .amplicon-tile.risk-critical { border-left-color: var(--bad); }
    .change-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px 12px;
      margin: 10px 0 0;
      color: var(--muted);
      font-size: 12px;
    }
    .change-swatch {
      display: inline-block;
      width: 11px;
      height: 11px;
      border-radius: 2px;
      margin-right: 5px;
      vertical-align: -1px;
      border: 1px solid rgba(23, 32, 42, 0.18);
    }
    .documentation {
      max-width: 1160px;
      display: grid;
      gap: 16px;
    }
    .documentation h2 {
      margin: 0;
      font-size: 24px;
    }
    .documentation p,
    .documentation li,
    .documentation span {
      color: var(--muted);
      line-height: 1.62;
    }
    .documentation p { margin: 0; }
    .documentation ul,
    .documentation ol {
      margin: 0;
      padding-left: 20px;
    }
    .doc-lead {
      max-width: 880px;
      padding: 16px 18px;
      background: #eef7f8;
      border: 1px solid #c8e4e8;
      border-left: 5px solid var(--accent);
      border-radius: 8px;
      color: #334155;
      font-size: 15px;
      line-height: 1.65;
    }
    .doc-section {
      min-width: 0;
      background: var(--panel);
      border: 1px solid var(--line);
      border-left: 4px solid #9fb3c8;
      border-radius: 8px;
      overflow: hidden;
    }
    .doc-section-priority { border-left-color: var(--accent); }
    .doc-section > h3 {
      margin: 0;
      padding: 13px 18px;
      background: #f8fafc;
      border-bottom: 1px solid var(--line);
      color: var(--ink);
      font-size: 16px;
    }
    .doc-section > p,
    .doc-section > ul,
    .doc-section > ol {
      margin: 16px 18px;
    }
    .doc-section > p + ul,
    .doc-section > p + ol {
      margin-top: -6px;
    }
    .doc-section li + li { margin-top: 8px; }
    .doc-faq-list,
    .doc-risk-list {
      list-style: none;
      padding-left: 0;
    }
    .doc-faq-list {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 10px;
    }
    .doc-faq-list li {
      display: grid;
      gap: 6px;
      min-height: 100%;
      padding: 12px 14px;
      background: #fbfcfe;
      border: 1px solid var(--line);
      border-left: 4px solid #c8e4e8;
      border-radius: 7px;
    }
    .doc-faq-list strong {
      color: var(--ink);
      font-size: 14px;
      line-height: 1.35;
    }
    .doc-faq-list span {
      font-size: 14px;
    }
    .doc-risk-list {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 10px;
    }
    .doc-risk-list li {
      padding: 12px 14px;
      background: #fbfcfe;
      border: 1px solid var(--line);
      border-left-width: 5px;
      border-radius: 7px;
      line-height: 1.55;
    }
    .doc-risk-critical { border-left-color: var(--bad); }
    .doc-risk-high { border-left-color: #dd6b20; }
    .doc-risk-watch { border-left-color: #d69e2e; }
    .doc-risk-low { border-left-color: var(--good); }
    .doc-two-column {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(330px, 1fr));
      gap: 14px;
      align-items: stretch;
    }
    .doc-step-list li::marker {
      color: var(--accent);
      font-weight: 700;
    }
    .empty { padding: 18px; color: var(--muted); }
    @media (max-width: 700px) {
      header, main { padding-left: 16px; padding-right: 16px; }
      .distribution-row { grid-template-columns: 52px 1fr; }
      .distribution-row span:last-child { grid-column: 2; }
      .note-row { grid-template-columns: 1fr; gap: 3px; }
      .doc-section { padding: 14px; }
    }
  </style>
</head>
<body>
  <!-- Visible wording is loaded from report_text/english.json and norwegian.json. -->
  <header>
    <h1 data-i18n="report_title">Primer binding-site report</h1>
    <p class="meta" data-i18n="report_meta">Shows how well each primer matches the uploaded sequences. Use the report to find primers that may need closer review.</p>
  </header>
  <main>
    <div class="language-toggle" data-i18n-aria-label="language_label" aria-label="Report language">
      <button class="language-button active" type="button" data-language="en">English</button>
      <button class="language-button" type="button" data-language="no">Norsk</button>
    </div>
    <nav class="report-tabs" data-i18n-aria-label="report_tabs_label" aria-label="Report sections">
      <button class="tab-button active" type="button" data-tab="overview" data-i18n="tab_overview">Results</button>
      <button class="tab-button" type="button" data-tab="documentation" data-i18n="tab_documentation">Help and documentation</button>
    </nav>
    <section class="tab-panel" id="overview-tab">
      <section class="filters" data-i18n-aria-label="report_filters_label" aria-label="Filter report results">
        <label><span data-i18n="filter_virus">Organism or virus</span><select id="filter-virus"></select></label>
        <label><span data-i18n="filter_primer">Primer</span><select id="filter-primer"></select></label>
        <label><span data-i18n="filter_fasta">Input FASTA file</span><select id="filter-fasta"></select></label>
        <label><span data-i18n="filter_risk">Review level</span><select id="filter-risk"></select></label>
      </section>
      <section class="cards" id="overview-cards"></section>
      <section id="ngs-panel-overview"></section>
      <h2 data-i18n="primer_overview_title">Primer overview</h2>
      <div class="report-note">
        <div class="note-row"><strong data-i18n="overview_note_label">How to read</strong><span data-i18n="primer_overview_note">First review primers marked Watch, High or Critical. Click a primer name to open its details. Check which sequences have no hit or several mismatches, and whether mismatches occur near either primer end.</span></div>
        <div class="note-row"><strong data-i18n="risk_note_label">Important</strong><span data-i18n="risk_criteria_note">The review level is an automatic screening signal based on the filtered data. It does not by itself mean that the assay has failed. Confirm important findings with assay performance data and laboratory context.</span></div>
      </div>
      <div class="legend" id="risk-legend"></div>
      <div class="table-wrap"><table id="summary-table"></table></div>
      <h2 data-i18n="primer_investigation_title">Primer details</h2>
      <div class="report-note"><div class="note-row"><strong data-i18n="primer_investigation_note_label">How to use this section</strong><span data-i18n="primer_investigation_note">For each primer, first check the summary values and the mismatch-position chart. Then review the mismatch distribution and the individual sequences.</span></div></div>
      <section id="primer-panels"></section>
    </section>
    <section class="tab-panel documentation" id="documentation-tab" hidden>
      <h2 data-i18n="doc_title">Help and documentation</h2>
      <p class="doc-lead" data-i18n="doc_intro"></p>

      <section class="doc-section doc-section-priority">
        <h3 data-i18n="doc_use_title"></h3>
        <ol class="doc-step-list">
          <li data-i18n="doc_use_1"></li>
          <li data-i18n="doc_use_2"></li>
          <li data-i18n="doc_use_3"></li>
          <li data-i18n="doc_use_4"></li>
          <li data-i18n="doc_use_5"></li>
          <li data-i18n="doc_use_6"></li>
        </ol>
      </section>

      <section class="doc-section doc-section-priority">
        <h3 data-i18n="doc_faq_title"></h3>
        <ul class="doc-faq-list">
          <li><strong data-i18n="doc_faq_1_q"></strong><span data-i18n="doc_faq_1"></span></li>
          <li><strong data-i18n="doc_faq_2_q"></strong><span data-i18n="doc_faq_2"></span></li>
          <li><strong data-i18n="doc_faq_3_q"></strong><span data-i18n="doc_faq_3"></span></li>
          <li><strong data-i18n="doc_faq_4_q"></strong><span data-i18n="doc_faq_4"></span></li>
          <li><strong data-i18n="doc_faq_5_q"></strong><span data-i18n="doc_faq_5"></span></li>
          <li><strong data-i18n="doc_faq_6_q"></strong><span data-i18n="doc_faq_6"></span></li>
          <li><strong data-i18n="doc_faq_7_q"></strong><span data-i18n="doc_faq_7"></span></li>
        </ul>
      </section>

      <section class="doc-section">
        <h3 data-i18n="doc_glossary_title"></h3>
        <ul>
          <li data-i18n="doc_term_hit"></li>
          <li data-i18n="doc_term_no_hit"></li>
          <li data-i18n="doc_term_mismatch"></li>
          <li data-i18n="doc_term_identity"></li>
          <li data-i18n="doc_term_terminal"></li>
        </ul>
      </section>

      <section class="doc-section">
        <h3 data-i18n="doc_cards_title"></h3>
        <p data-i18n="doc_filters_text"></p>
        <ul>
          <li data-i18n="doc_cards_1"></li>
          <li data-i18n="doc_cards_2"></li>
          <li data-i18n="doc_cards_3"></li>
          <li data-i18n="doc_cards_4"></li>
          <li data-i18n="doc_cards_5"></li>
        </ul>
      </section>

      <section class="doc-section">
        <h3 data-i18n="doc_overview_title"></h3>
        <p data-i18n="doc_overview_text"></p>
        <ul>
          <li data-i18n="doc_overview_1"></li>
          <li data-i18n="doc_overview_2"></li>
          <li data-i18n="doc_overview_3"></li>
          <li data-i18n="doc_overview_4"></li>
          <li data-i18n="doc_overview_5"></li>
          <li data-i18n="doc_overview_6"></li>
          <li data-i18n="doc_overview_7"></li>
        </ul>
      </section>

      <section class="doc-section">
        <h3 data-i18n="doc_risk_title"></h3>
        <p data-i18n="doc_risk_text"></p>
        <ul class="doc-risk-list">
          <li class="doc-risk-critical" data-i18n="doc_risk_critical"></li>
          <li class="doc-risk-high" data-i18n="doc_risk_high"></li>
          <li class="doc-risk-watch" data-i18n="doc_risk_watch"></li>
          <li class="doc-risk-low" data-i18n="doc_risk_low"></li>
        </ul>
      </section>

      <section class="doc-section">
        <h3 data-i18n="doc_panels_title"></h3>
        <p data-i18n="doc_panels_text"></p>
      </section>

      <div class="doc-two-column">
        <section class="doc-section">
          <h3 data-i18n="doc_chart_title"></h3>
          <p data-i18n="doc_chart_text"></p>
        </section>
        <section class="doc-section">
          <h3 data-i18n="doc_distribution_title"></h3>
          <p data-i18n="doc_distribution_text"></p>
        </section>
      </div>

      <div class="doc-two-column">
        <section class="doc-section">
          <h3 data-i18n="doc_detail_title"></h3>
          <p data-i18n="doc_detail_text"></p>
        </section>
        <section class="doc-section">
          <h3 data-i18n="doc_alignment_title"></h3>
          <p data-i18n="doc_alignment_text"></p>
        </section>
      </div>

      <section class="doc-section">
        <h3 data-i18n="doc_workflow_title"></h3>
        <ol class="doc-step-list">
          <li data-i18n="doc_workflow_1"></li>
          <li data-i18n="doc_workflow_2"></li>
          <li data-i18n="doc_workflow_3"></li>
          <li data-i18n="doc_workflow_4"></li>
          <li data-i18n="doc_workflow_5"></li>
        </ol>
      </section>

      <section class="doc-section">
        <h3 data-i18n="doc_ngs_title"></h3>
        <p data-i18n="doc_ngs_text"></p>
      </section>

      <section class="doc-section doc-section-priority">
        <h3 data-i18n="doc_limits_title"></h3>
        <ul>
          <li data-i18n="doc_limits_1"></li>
          <li data-i18n="doc_limits_2"></li>
          <li data-i18n="doc_limits_3"></li>
          <li data-i18n="doc_limits_4"></li>
          <li data-i18n="doc_limits_5"></li>
        </ul>
      </section>
    </section>
  </main>
  <div class="modal-backdrop" id="alignment-modal" role="dialog" aria-modal="true" aria-labelledby="alignment-modal-title">
    <div class="modal">
      <div class="modal-header">
        <div>
          <h2 id="alignment-modal-title" data-i18n="modal_title">Primer alignment</h2>
          <p class="section-note" id="alignment-modal-meta"></p>
        </div>
        <button class="modal-close" type="button" id="alignment-modal-close" data-i18n="close">Close</button>
      </div>
      <div id="alignment-modal-body"></div>
    </div>
  </div>
  <script id="report-data" type="application/json">__REPORT_DATA__</script>
  <script>
    const rows = JSON.parse(document.getElementById('report-data').textContent);
    const fields = __FIELDS_JSON__;
    const TRANSLATIONS = __TRANSLATIONS_JSON__;
    let currentLanguage = 'en';
    function t(key, replacements = {}) {
      let value = (TRANSLATIONS[currentLanguage] && TRANSLATIONS[currentLanguage][key]) || TRANSLATIONS.en[key] || key;
      for (const [name, replacement] of Object.entries(replacements)) {
        value = value.replaceAll('{' + name + '}', replacement);
      }
      return value;
    }
    function applyTranslations() {
      document.documentElement.lang = currentLanguage;
      document.title = t('report_title');
      document.querySelectorAll('[data-i18n]').forEach(element => { element.textContent = t(element.dataset.i18n); });
      document.querySelectorAll('[data-i18n-aria-label]').forEach(element => {
        element.setAttribute('aria-label', t(element.dataset.i18nAriaLabel));
      });
      document.querySelectorAll('[data-language]').forEach(button => button.classList.toggle('active', button.dataset.language === currentLanguage));
    }
    const RISK_THRESHOLDS = {
      highNoHitRate: 0.05,
      criticalNoHitRate: 0.20,
      highTwoPlusMismatchRate: 0.30,
      criticalThreePlusMismatchRate: 0.25,
      highTerminalMismatchRate: 0.05,
      criticalTerminalMismatchRate: 0.50,
      terminalBases: 5
    };
    const riskRank = { Low: 0, Watch: 1, High: 2, Critical: 3, Indeterminate: 4 };
    const distCategories = [
      { key: '0', labelKey: 'dist_0' },
      { key: '1', labelKey: 'dist_1' },
      { key: '2', labelKey: 'dist_2' },
      { key: '3', labelKey: 'dist_3' },
      { key: '4plus', labelKey: 'dist_4plus' },
      { key: 'no_hit', labelKey: 'dist_no_hit' },
      { key: 'indeterminate', labelKey: 'status_indeterminate' }
    ];
    const mismatchChangePalette = ['#1f7a8c', '#b42318', '#287d3c', '#b7791f', '#6f42c1', '#c2410c', '#0f766e', '#be185d', '#4d7c0f', '#0369a1', '#92400e', '#475569'];
    const mismatchChangeColorMap = new Map();
    const filters = {
      virus: document.getElementById('filter-virus'),
      primer: document.getElementById('filter-primer'),
      fasta: document.getElementById('filter-fasta'),
      risk: document.getElementById('filter-risk')
    };
    let sortState = { key: 'Risk', direction: -1 };
    const primerTableState = {};
    const detailColumns = [
      { key: 'Fasta_File', labelKey: 'detail_fasta' },
      { key: 'Virus_Type', labelKey: 'detail_virus' },
      { key: 'Subject_Sequence_ID', labelKey: 'detail_sample' },
      { key: 'Subject_Segment', labelKey: 'detail_segment' },
      { key: 'Hit_Status', labelKey: 'detail_status' },
      { key: 'Percent_Identity', labelKey: 'detail_identity' },
      { key: 'Mismatches', labelKey: 'detail_mismatches' },
      { key: 'Mismatch_Positions', labelKey: 'detail_positions' },
      { key: 'Mismatch_Details', labelKey: 'detail_bases' }
    ];

    function parseNumber(value) {
      if (value === '' || value === null || value === undefined || value === 'No hit') return null;
      const parsed = Number(String(value).replace(',', '.'));
      return Number.isFinite(parsed) ? parsed : null;
    }
    const numeric = parseNumber;
    function parseMismatchPositions(value) {
      if (!value) return [];
      return String(value).split(',')
        .map(part => Number(part.trim()))
        .filter(value => Number.isInteger(value) && value > 0);
    }
    function parseMismatchDetails(value) {
      const details = new Map();
      if (!value) return details;
      for (const part of String(value).split(',')) {
        const match = part.trim().match(/^(\\d+):(.+)$/);
        if (!match) continue;
        const pos = Number(match[1]);
        const change = match[2];
        if (!Number.isInteger(pos)) continue;
        if (!details.has(pos)) details.set(pos, new Map());
        details.get(pos).set(change, (details.get(pos).get(change) || 0) + 1);
      }
      return details;
    }
    function mergeMismatchDetailMaps(target, source) {
      for (const [pos, changes] of source.entries()) {
        if (!target.has(pos)) target.set(pos, new Map());
        for (const [change, count] of changes.entries()) {
          target.get(pos).set(change, (target.get(pos).get(change) || 0) + count);
        }
      }
    }
    function topMismatchChanges(changes) {
      if (!changes || !changes.size) return '';
      return [...changes.entries()]
        .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
        .slice(0, 3)
        .map(([change, count]) => change + ' x' + count)
        .join(', ');
    }
    function uniqueValues(key) {
      return [...new Set(rows.map(row => row[key] || '').filter(Boolean))].sort();
    }
    function fillSelect(select, values, labelKey) {
      select.innerHTML = '<option value="">' + t(labelKey) + '</option>' + values.map(value => '<option>' + escapeHtml(value) + '</option>').join('');
    }
    function fillRiskSelect() {
      const riskValues = ['Low', 'Watch', 'High', 'Critical', 'Indeterminate'];
      filters.risk.innerHTML = '<option value="">' + t('all_risks') + '</option>' +
        riskValues.map(risk => '<option value="' + risk + '">' + escapeHtml(t(risk.toLowerCase())) + '</option>').join('');
    }
    function escapeHtml(value) {
      return String(value).replace(/[&<>"']/g, character => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[character]));
    }
    function formatPercent(value) {
      return Number.isFinite(value) ? (value * 100).toFixed(1) + '%' : t('not_available');
    }
    function isNoHit(row) {
      return row.Hit_Status === 'no_hit' || row.Percent_Identity === 'No hit';
    }
    function statusLabel(row) {
      if (row.Hit_Status === 'indeterminate') return t('status_indeterminate');
      return isNoHit(row) ? t('status_no_hit') : t('status_hit');
    }
    function displayValue(value) {
      return value === null || value === undefined || value === '' || value === 'n/a' ? t('not_available') : value;
    }
    function formatIdentityValue(value) {
      const parsed = parseNumber(value);
      return parsed === null ? escapeHtml(displayValue(value)) : parsed.toFixed(2);
    }
    function sampleBaseFromChange(change) {
      const text = String(change || '').trim();
      const parts = text.split('>');
      return (parts.length > 1 ? parts[parts.length - 1] : text).toUpperCase() || 'other';
    }
    function colorForMismatchChange(change) {
      const key = String(change || 'other');
      if (!mismatchChangeColorMap.has(key)) {
        mismatchChangeColorMap.set(key, mismatchChangePalette[mismatchChangeColorMap.size % mismatchChangePalette.length]);
      }
      return mismatchChangeColorMap.get(key);
    }
    function currentRows() {
      const baseRows = rows.filter(row => {
        if (filters.virus.value && row.Virus_Type !== filters.virus.value) return false;
        if (filters.primer.value && row.Primer_Name !== filters.primer.value) return false;
        if (filters.fasta.value && row.Fasta_File !== filters.fasta.value) return false;
        return true;
      });
      if (!filters.risk.value) return baseRows;
      const primersForRisk = new Set(
        groupByPrimer(baseRows)
          .filter(([, group]) => calculateRisk(calculatePrimerStats(group)) === filters.risk.value)
          .map(([primer]) => primer)
      );
      return baseRows.filter(row => primersForRisk.has(row.Primer_Name));
    }
    function average(values) {
      const numericValues = values.map(numeric).filter(value => value !== null);
      if (!numericValues.length) return '';
      return (numericValues.reduce((sum, value) => sum + value, 0) / numericValues.length).toFixed(2);
    }
    function groupByPrimer(data) {
      const grouped = new Map();
      for (const row of data) {
        if (!grouped.has(row.Primer_Name)) grouped.set(row.Primer_Name, []);
        grouped.get(row.Primer_Name).push(row);
      }
      return [...grouped.entries()].sort((a, b) => a[0].localeCompare(b[0]));
    }
    function roleForRow(row) {
      const explicitRole = String(row.Primer_Role || '').toLowerCase();
      if (explicitRole === 'probe') return 'probe';
      const tokens = String(row.Primer_Name || '').toUpperCase().split(/[^A-Z0-9]+/).filter(Boolean);
      return tokens.some(token => token.includes('PROBE') || ['TM', 'VIC2', 'YAM2'].includes(token)) ? 'probe' : 'primer';
    }
    function terminalMismatchPosition(position, sequenceLength, role) {
      const terminalBases = RISK_THRESHOLDS.terminalBases;
      return role === 'probe' ? position <= terminalBases : position > sequenceLength - terminalBases;
    }
    function terminalLabelKey(role) {
      return role === 'probe' ? 'terminal_share_probe' : 'terminal_share_primer';
    }
    function chartNoteKey(role) {
      return role === 'probe' ? 'chart_note_use_probe' : 'chart_note_use_primer';
    }
    function calculateMismatchDistribution(group) {
      const distribution = { '0': 0, '1': 0, '2': 0, '3': 0, '4plus': 0, no_hit: 0, indeterminate: 0 };
      for (const row of group) {
        if (row.Hit_Status === 'indeterminate') {
          distribution.indeterminate += 1;
          continue;
        }
        if (isNoHit(row)) {
          distribution.no_hit += 1;
          continue;
        }
        const mismatches = parseNumber(row.Mismatches);
        if (mismatches === null) continue;
        if (mismatches >= 4) distribution['4plus'] += 1;
        else distribution[String(Math.max(0, Math.trunc(mismatches)))] += 1;
      }
      return distribution;
    }
    function calculateMismatchPositionCounts(group) {
      const hits = group.filter(row => !isNoHit(row) && row.Hit_Status !== 'indeterminate');
      const primerSequence = group.find(row => row.Primer_Sequence)?.Primer_Sequence || '';
      const counts = Array.from({ length: primerSequence.length }, () => 0);
      const detailsByPosition = new Map();
      for (const row of hits) {
        for (const pos of parseMismatchPositions(row.Mismatch_Positions)) {
          if (pos >= 1 && pos <= counts.length) counts[pos - 1] += 1;
        }
        mergeMismatchDetailMaps(detailsByPosition, parseMismatchDetails(row.Mismatch_Details));
      }
      return { primerSequence, counts, detailsByPosition, hitRows: hits.length };
    }
    function calculatePrimerStats(group) {
      const totalRows = group.length;
      const role = roleForRow(group[0] || {});
      const hitRows = group.filter(row => !isNoHit(row) && row.Hit_Status !== 'indeterminate');
      const noHits = group.filter(isNoHit).length;
      const indeterminate = group.filter(row => row.Hit_Status === 'indeterminate').length;
      const mismatchValues = hitRows.map(row => parseNumber(row.Mismatches)).filter(value => value !== null);
      const identityValues = hitRows.map(row => parseNumber(row.Percent_Identity)).filter(value => value !== null);
      const avgMismatches = mismatchValues.length ? mismatchValues.reduce((a, b) => a + b, 0) / mismatchValues.length : 0;
      const maxMismatches = mismatchValues.length ? Math.max(...mismatchValues) : 0;
      const avgPercentIdentity = identityValues.length ? identityValues.reduce((a, b) => a + b, 0) / identityValues.length : null;
      const twoPlus = mismatchValues.filter(value => value >= 2).length;
      const threePlus = mismatchValues.filter(value => value >= 3).length;
      const positionData = calculateMismatchPositionCounts(group);
      const terminalBases = RISK_THRESHOLDS.terminalBases;
      const terminalMismatchHits = hitRows.filter(row => {
        return parseMismatchPositions(row.Mismatch_Positions).some(pos => {
          return terminalMismatchPosition(pos, positionData.counts.length, role);
        });
      }).length;
      return {
        totalRows,
        indeterminate,
        hits: hitRows.length,
        noHits,
        noHitRate: hitRows.length + noHits ? noHits / (hitRows.length + noHits) : 0,
        avgMismatches,
        maxMismatches,
        avgPercentIdentity,
        twoPlusMismatchRate: hitRows.length ? twoPlus / hitRows.length : 0,
        threePlusMismatchRate: hitRows.length ? threePlus / hitRows.length : 0,
        terminalMismatchRate: hitRows.length ? terminalMismatchHits / hitRows.length : 0,
        terminalRole: role,
        distribution: calculateMismatchDistribution(group)
      };
    }
    function calculateRisk(stats) {
      if (stats.indeterminate && !stats.hits && !stats.noHits) return 'Indeterminate';
      if (
        stats.noHitRate >= RISK_THRESHOLDS.criticalNoHitRate ||
        stats.threePlusMismatchRate >= RISK_THRESHOLDS.criticalThreePlusMismatchRate ||
        stats.terminalMismatchRate >= RISK_THRESHOLDS.criticalTerminalMismatchRate
      ) return 'Critical';
      if (
        stats.noHitRate >= RISK_THRESHOLDS.highNoHitRate ||
        stats.twoPlusMismatchRate >= RISK_THRESHOLDS.highTwoPlusMismatchRate ||
        stats.terminalMismatchRate >= RISK_THRESHOLDS.highTerminalMismatchRate
      ) return 'High';
      if (stats.maxMismatches >= 2 || stats.noHits > 0) return 'Watch';
      return 'Low';
    }
    function riskExplanation(stats) {
      return [
        t('status_indeterminate') + ': ' + stats.indeterminate,
        t('no_hit_rate') + ': ' + formatPercent(stats.noHitRate),
        t('max_mismatches') + ': ' + stats.maxMismatches,
        t('two_plus_rate') + ': ' + formatPercent(stats.twoPlusMismatchRate),
        t('three_plus_rate') + ': ' + formatPercent(stats.threePlusMismatchRate),
        t(terminalLabelKey(stats.terminalRole)) + ': ' + formatPercent(stats.terminalMismatchRate)
      ].join(' | ');
    }
    function riskBadge(risk, title) {
      return '<span class="risk-badge risk-' + risk.toLowerCase() + '" title="' + escapeHtml(title) + '">' + escapeHtml(t(risk.toLowerCase())) + '</span>';
    }
    function renderCards(data) {
      const hits = data.filter(row => row.Hit_Status === 'hit').length;
      const noHits = data.filter(row => row.Hit_Status === 'no_hit').length;
      const primers = new Set(data.map(row => row.Primer_Name)).size;
      const samples = new Set(data.map(row => row.Subject_Sequence_ID)).size;
      const avgIdentity = average(data.map(row => row.Percent_Identity)) || t('not_available');
      document.getElementById('overview-cards').innerHTML = [
        [t('rows'), data.length],
        [t('primers'), primers],
        [t('samples'), samples],
        [t('hits'), hits],
        [t('no_hits'), noHits],
        [t('status_indeterminate'), data.filter(row => row.Hit_Status === 'indeterminate').length],
        [t('avg_identity'), avgIdentity]
      ].map(([label, value]) => '<div class="card"><h3>' + label + '</h3><div class="metric">' + value + '</div></div>').join('');
    }
    function renderRiskLegend() {
      document.getElementById('risk-legend').innerHTML =
        '<span>' + t('risk') + ':</span>' +
        riskBadge('Low', t('risk_low_title')) +
        riskBadge('Watch', t('risk_watch_title')) +
        riskBadge('High', t('risk_high_title')) +
        riskBadge('Critical', t('risk_critical_title'));
    }
    function ampliconName(primerName) {
      const match = String(primerName || '').match(/^(.*)_(LEFT|RIGHT)(?:_\\d+)?$/i);
      return match ? match[1] : '';
    }
    function primerDirection(primerName) {
      const match = String(primerName || '').match(/_(LEFT|RIGHT)(?:_\\d+)?$/i);
      return match ? match[1].toUpperCase() : '';
    }
    function ampliconOrder(name) {
      const numbers = String(name).match(/\\d+/g);
      return numbers && numbers.length ? Number(numbers[numbers.length - 1]) : Number.MAX_SAFE_INTEGER;
    }
    function renderNgsPanelOverview(data) {
      const target = document.getElementById('ngs-panel-overview');
      const ngsRows = data.filter(row => row.Assay_Type === 'ngs' || ampliconName(row.Primer_Name));
      if (!ngsRows.length) {
        target.innerHTML = '';
        return;
      }
      const grouped = new Map();
      for (const row of ngsRows) {
        const amplicon = ampliconName(row.Primer_Name);
        if (!amplicon) continue;
        const key = (row.Assay_ID || '') + '|' + amplicon;
        if (!grouped.has(key)) grouped.set(key, { amplicon, assayId: row.Assay_ID || '', assayName: row.Assay_Name || row.Assay_ID || t('unnamed_ngs'), rows: [], primers: new Set(), pools: new Set() });
        const group = grouped.get(key);
        group.rows.push(row);
        group.primers.add(row.Primer_Name);
        if (row.Primer_Pool) group.pools.add(row.Primer_Pool);
      }
      const amplicons = [...grouped.values()].sort((a, b) =>
        a.assayId.localeCompare(b.assayId) || ampliconOrder(a.amplicon) - ampliconOrder(b.amplicon) || a.amplicon.localeCompare(b.amplicon)
      );
      if (!amplicons.length) {
        target.innerHTML = '';
        return;
      }
      function renderAmpliconTile(group) {
        const primerGroups = groupByPrimer(group.rows);
        const directionRisks = { LEFT: [], RIGHT: [] };
        for (const [primer, primerRows] of primerGroups) {
          const direction = primerDirection(primer);
          if (directionRisks[direction]) directionRisks[direction].push(calculateRisk(calculatePrimerStats(primerRows)));
        }
        const bestDirectionRisk = direction => directionRisks[direction].sort((a, b) => riskRank[a] - riskRank[b])[0];
        const bestLeft = bestDirectionRisk('LEFT');
        const bestRight = bestDirectionRisk('RIGHT');
        const leftViable = bestLeft !== undefined && riskRank[bestLeft] < riskRank.High;
        const rightViable = bestRight !== undefined && riskRank[bestRight] < riskRank.High;
        const risk = leftViable && rightViable
          ? (riskRank[bestLeft] >= riskRank[bestRight] ? bestLeft : bestRight)
          : (bestLeft === 'Indeterminate' || bestRight === 'Indeterminate' ? 'Indeterminate' : 'Critical');
        const stats = calculatePrimerStats(group.rows);
        const directionStatus = t('direction_forward') + ': ' + (leftViable ? t('viable') + ' (' + t(bestLeft.toLowerCase()) + ')' : t('no_viable_primer')) +
          ' | ' + t('direction_reverse') + ': ' + (rightViable ? t('viable') + ' (' + t(bestRight.toLowerCase()) + ')' : t('no_viable_primer'));
        const title = group.assayName + ' | ' + group.amplicon + ' | ' + directionStatus + ' | ' + t('primers_label') + ': ' + [...group.primers].join(', ');
        return '<div class="amplicon-tile risk-' + risk.toLowerCase() + '" title="' + escapeHtml(title) + '">' +
          '<strong>' + escapeHtml(group.amplicon) + '</strong>' +
          riskBadge(risk, title) +
          '<div>' + group.primers.size + ' ' + (group.primers.size === 1 ? t('primer_count_singular') : t('primer_count_plural')) +
          (group.pools.size ? ' · ' + t('pool') + ' ' + escapeHtml([...group.pools].join(', ')) : '') + '</div>' +
          '<div>' + t('direction_forward') + ': ' + (leftViable ? t(bestLeft.toLowerCase()) : t('missing')) + '</div>' +
          '<div>' + t('direction_reverse') + ': ' + (rightViable ? t(bestRight.toLowerCase()) : t('missing')) + '</div></div>';
      }
      const byPanel = new Map();
      for (const group of amplicons) {
        const key = group.assayId || group.assayName;
        if (!byPanel.has(key)) byPanel.set(key, { name: group.assayName, id: group.assayId, amplicons: [] });
        byPanel.get(key).amplicons.push(group);
      }
      const panels = [...byPanel.values()].map(panel =>
        '<section class="primer-panel"><h3>' + escapeHtml(panel.name) +
        (panel.id && panel.id !== panel.name ? ' <span class="section-note">(' + escapeHtml(panel.id) + ')</span>' : '') + '</h3>' +
        '<div class="amplicon-map">' + panel.amplicons.map(renderAmpliconTile).join('') + '</div></section>'
      ).join('');
      target.innerHTML = '<h2>' + t('ngs_title') + '</h2>' +
        '<div class="report-note"><div class="note-row"><strong>' + t('ngs_note_label') + '</strong><span>' + t('ngs_note') + '</span></div></div>' +
        '<div class="legend">' + riskBadge('Low', t('no_current_warning')) + riskBadge('Watch', t('review')) + riskBadge('High', t('elevated_risk')) + riskBadge('Critical', t('strong_risk')) + '</div>' +
        panels;
    }
    function distributionChart(distribution) {
      const total = Object.values(distribution).reduce((sum, value) => sum + value, 0);
      if (!total) return '<div class="empty">' + t('no_data_mismatch_count') + '</div>';
      return '<div class="mismatch-distribution">' + distCategories.map(category => {
        const count = distribution[category.key] || 0;
        const pct = total ? count / total : 0;
        const label = t(category.labelKey);
        return '<div class="distribution-row" title="' + escapeHtml(label + ': ' + count + ' ' + t('rows_lower') + ', ' + formatPercent(pct)) + '">' +
          '<strong>' + escapeHtml(label) + '</strong>' +
          '<div class="distribution-track"><div class="distribution-fill dist-' + category.key + '" style="width:' + Math.max(0, pct * 100) + '%"></div></div>' +
          '<span>' + count + ' (' + formatPercent(pct) + ')</span>' +
          '</div>';
      }).join('') + '</div>';
    }
    function renderSequenceMap(group) {
      const { primerSequence, counts, detailsByPosition, hitRows } = calculateMismatchPositionCounts(group);
      const role = roleForRow(group[0] || {});
      if (!primerSequence || !counts.length || !hitRows || counts.every(count => count === 0)) {
        return '<div class="plot-card"><h3>' + t('chart_title') + '</h3><div class="empty">' + t('chart_empty') + '</div></div>';
      }
      const width = Math.max(620, counts.length * 28 + 110);
      const height = 164;
      const plotHeight = 68;
      const baseline = 92;
      const left = 76;
      const right = 20;
      const cellWidth = (width - left - right) / counts.length;
      const baseCountsByPosition = new Map();
      const observedSampleBases = new Set();
      let maxPct = 0.01;
      for (const [pos, changes] of detailsByPosition.entries()) {
        const baseCounts = new Map();
        for (const [change, changeCount] of changes.entries()) {
          const sampleBase = sampleBaseFromChange(change);
          baseCounts.set(sampleBase, (baseCounts.get(sampleBase) || 0) + changeCount);
          observedSampleBases.add(sampleBase);
        }
        baseCountsByPosition.set(pos, baseCounts);
        const positionTotal = [...baseCounts.values()].reduce((sum, value) => sum + value, 0);
        maxPct = Math.max(maxPct, hitRows ? positionTotal / hitRows : 0);
      }
      const terminalBases = RISK_THRESHOLDS.terminalBases;
      const yTicks = [0, maxPct / 2, maxPct].filter((value, index, arr) => arr.findIndex(other => Math.abs(other - value) < 0.000001) === index);
      const tickMarks = yTicks.map(value => {
        const y = baseline - (value / maxPct) * plotHeight;
        return '<g><line x1="' + left + '" y1="' + y.toFixed(1) + '" x2="' + (width - right) + '" y2="' + y.toFixed(1) + '" stroke="#e2e8f0"></line>' +
          '<text class="axis-label" x="' + (left - 8) + '" y="' + (y + 3).toFixed(1) + '" text-anchor="end">' + formatPercent(value) + '</text></g>';
      }).join('');
      const cells = counts.map((count, index) => {
        const pos = index + 1;
        const base = primerSequence[index] || '';
        const baseCounts = baseCountsByPosition.get(pos) || new Map();
        const totalCount = [...baseCounts.values()].reduce((sum, value) => sum + value, 0);
        const totalPct = hitRows ? totalCount / hitRows : 0;
        const x = left + index * cellWidth;
        const terminal = terminalMismatchPosition(pos, counts.length, role);
        let yCursor = baseline;
        const segments = [...baseCounts.entries()]
          .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
          .map(([sampleBase, sampleBaseCount]) => {
            const segmentPct = hitRows ? sampleBaseCount / hitRows : 0;
            const segmentHeight = segmentPct ? Math.max(2, (segmentPct / maxPct) * plotHeight) : 0;
            yCursor -= segmentHeight;
            return '<rect x="' + x.toFixed(1) + '" y="' + yCursor.toFixed(1) + '" width="' + Math.max(3, cellWidth - 5).toFixed(1) + '" height="' + segmentHeight.toFixed(1) + '" fill="' + colorForMismatchChange(sampleBase) + '"><title>' +
              escapeHtml(t('position') + ' ' + pos + ', ' + t('sample_base') + ' ' + sampleBase + ': ' + formatPercent(segmentPct) + ' (' + sampleBaseCount + '/' + hitRows + ' ' + t('hit_rows') + ')') +
              '</title></rect>';
          }).join('');
        const detailText = [...baseCounts.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0])).map(([sampleBase, sampleBaseCount]) => sampleBase + ' x' + sampleBaseCount).join(', ');
        const title = t('position') + ' ' + pos + ', ' + t('base') + ' ' + base + ': ' + formatPercent(totalPct) + ' ' + t('with_any_mismatch') + ' (' + totalCount + '/' + hitRows + ' ' + t('hit_rows') + ')' + (detailText ? '; ' + detailText : '');
        return '<g><title>' + escapeHtml(title) + '</title>' +
          (segments || '<rect x="' + x.toFixed(1) + '" y="' + (baseline - 1) + '" width="' + Math.max(3, cellWidth - 5).toFixed(1) + '" height="1" fill="#d8dee9"></rect>') +
          (terminal ? '<rect x="' + x.toFixed(1) + '" y="' + (baseline + 3) + '" width="' + Math.max(3, cellWidth - 5).toFixed(1) + '" height="3" fill="#b42318" opacity="0.5"></rect>' : '') +
          '<text class="base-label" x="' + (x + cellWidth / 2).toFixed(1) + '" y="116" text-anchor="middle">' + escapeHtml(base) + '</text>' +
          (pos === 1 || pos === counts.length || pos % 5 === 0 ? '<text class="axis-label" x="' + (x + cellWidth / 2).toFixed(1) + '" y="138" text-anchor="middle">' + pos + '</text>' : '') +
          '</g>';
      }).join('');
      const legend = [...observedSampleBases].sort().map(sampleBase => '<span><span class="change-swatch" style="background:' + colorForMismatchChange(sampleBase) + '"></span>' + escapeHtml(t('sample_base') + ' ' + sampleBase) + '</span>').join('');
      return '<div class="plot-card"><h3>' + t('chart_title') + '</h3>' +
        '<div class="report-note report-note-compact">' +
        '<div class="note-row"><strong>' + t('chart_note_what_label') + '</strong><span>' + t('chart_note_what') + '</span></div>' +
        '<div class="note-row"><strong>' + t('chart_note_use_label') + '</strong><span>' + t(chartNoteKey(role)) + '</span></div>' +
        '</div>' +
        '<div class="primer-sequence-map"><svg class="sequence-svg" viewBox="0 0 ' + width + ' ' + height + '" role="img" aria-label="' + escapeHtml(t('chart_title')) + '">' +
        '<text class="axis-title" transform="translate(14 55) rotate(-90)" text-anchor="middle">' + escapeHtml(t('percent_axis')) + '</text>' +
        tickMarks +
        '<text class="five-prime" x="' + (left - 25) + '" y="116">5&apos;</text><text class="three-prime" x="' + (width - 18) + '" y="116">3&apos;</text>' +
        '<line x1="' + left + '" y1="' + baseline + '" x2="' + (width - right) + '" y2="' + baseline + '" stroke="#9aa6b2"></line>' +
        cells +
        '</svg></div>' +
        (legend ? '<div class="change-legend">' + legend + '</div>' : '') +
        '</div>';
    }
    function summaryRows(data) {
      return groupByPrimer(data).map(([primer, group]) => {
        const stats = calculatePrimerStats(group);
        const risk = calculateRisk(stats);
        return {
          Primer_Name: primer,
          Virus_Type: [...new Set(group.map(row => row.Virus_Type))].join(', '),
          Primer_Segment: [...new Set(group.map(row => row.Primer_Segment).filter(Boolean))].join(', '),
          Risk: risk,
          Risk_Title: riskExplanation(stats),
          Risk_Rank: riskRank[risk],
          Rows: stats.totalRows,
          Hits: stats.hits,
          No_Hits: stats.noHits,
          No_Hit_Rate: formatPercent(stats.noHitRate),
          Max_Mismatches: stats.maxMismatches,
          Avg_Percent_Identity: stats.avgPercentIdentity === null ? t('not_available') : stats.avgPercentIdentity.toFixed(2),
          TwoPlus_Rate: formatPercent(stats.twoPlusMismatchRate),
          ThreePlus_Rate: formatPercent(stats.threePlusMismatchRate),
          Terminal_Mismatch_Share: formatPercent(stats.terminalMismatchRate)
        };
      });
    }
    function safeDomId(value) {
      return String(value).replace(/[^A-Za-z0-9_-]/g, '_');
    }
    function primerPanelId(primer) {
      return 'primer-panel-' + safeDomId(primer);
    }
    function primerTableId(primer) {
      return 'primer-table-' + safeDomId(primer);
    }
    function primerState(primer) {
      if (!primerTableState[primer]) {
        primerTableState[primer] = { expanded: false, sortKey: null, direction: -1 };
      }
      return primerTableState[primer];
    }
    function valueForSort(row, key) {
      if (key === 'Mismatches') return parseNumber(row.Mismatches) ?? Number.NEGATIVE_INFINITY;
      if (key === 'Percent_Identity') return parseNumber(row.Percent_Identity) ?? Number.NEGATIVE_INFINITY;
      if (key === 'Hit_Status') return statusLabel(row);
      return row[key] || '';
    }
    function sortedDetailRows(group, state) {
      if (!state.sortKey) {
        return [...group].sort((a, b) => {
          const aNoHit = isNoHit(a);
          const bNoHit = isNoHit(b);
          if (aNoHit !== bNoHit) return aNoHit ? -1 : 1;
          if (!aNoHit && !bNoHit) {
            const mismatchDelta = (parseNumber(b.Mismatches) ?? Number.NEGATIVE_INFINITY) - (parseNumber(a.Mismatches) ?? Number.NEGATIVE_INFINITY);
            if (mismatchDelta) return mismatchDelta;
          }
          return String(a.Subject_Sequence_ID || '').localeCompare(String(b.Subject_Sequence_ID || ''));
        });
      }
      return [...group].sort((a, b) => {
        const left = valueForSort(a, state.sortKey);
        const right = valueForSort(b, state.sortKey);
        if (typeof left === 'number' && typeof right === 'number') return (left - right) * state.direction;
        return String(left).localeCompare(String(right)) * state.direction;
      });
    }
    function rowKey(row) {
      return [
        row.Fasta_File || '',
        row.Primer_Name || '',
        row.Subject_Sequence_ID || '',
        row.Query_Start || '',
        row.Query_End || ''
      ].join('||');
    }
    function alignmentMatchLine(query, subject) {
      const length = Math.max(query.length, subject.length);
      let line = '';
      for (let i = 0; i < length; i += 1) {
        const q = query[i] || ' ';
        const s = subject[i] || ' ';
        line += q && s && q !== ' ' && s !== ' ' && basesMatchForReport(q, s) ? '|' : ' ';
      }
      return line;
    }
    const IUPAC_BASES_FOR_REPORT = {
      A: ['A'],
      C: ['C'],
      G: ['G'],
      T: ['T'],
      U: ['T'],
      R: ['A', 'G'],
      Y: ['C', 'T'],
      S: ['G', 'C'],
      W: ['A', 'T'],
      K: ['G', 'T'],
      M: ['A', 'C'],
      B: ['C', 'G', 'T'],
      D: ['A', 'G', 'T'],
      H: ['A', 'C', 'T'],
      V: ['A', 'C', 'G'],
      N: ['A', 'C', 'G', 'T']
    };
    function basesMatchForReport(a, b) {
      const leftBase = String(a).toUpperCase();
      const rightBase = String(b).toUpperCase();
      if (leftBase === '-' || rightBase === '-') return false;
      const leftOptions = IUPAC_BASES_FOR_REPORT[leftBase] || [leftBase];
      const rightOptions = IUPAC_BASES_FOR_REPORT[rightBase] || [rightBase];
      return leftOptions.some(base => rightOptions.includes(base));
    }
    function showAlignmentModal(row) {
      const modal = document.getElementById('alignment-modal');
      const title = document.getElementById('alignment-modal-title');
      const meta = document.getElementById('alignment-modal-meta');
      const body = document.getElementById('alignment-modal-body');
      const query = row.Query_Alignment || row.Primer_Sequence || '';
      const subject = row.Subject_Alignment || '';
      const strand = Number(row.Subject_Start) <= Number(row.Subject_End) ? t('strand_forward') : t('strand_reverse');
      const hitCoordinates = isNoHit(row) ? '' : t('local_hit_coordinates', { strand, start: row.Subject_Start, end: row.Subject_End });
      title.textContent = row.Subject_Sequence_ID || t('modal_title');
      meta.textContent = [row.Primer_Name, row.Fasta_File, statusLabel(row), hitCoordinates, row.Mismatches ? row.Mismatches + ' ' + t('mismatches_word') : '', row.Mismatch_Details || ''].filter(Boolean).join(' | ');
      if (isNoHit(row) || !subject) {
        body.innerHTML = '<div class="empty">' + t('no_alignment') + '</div>';
      } else {
        const matchLine = alignmentMatchLine(query, subject);
        body.innerHTML = '<div class="alignment-view">' +
          '<div class="alignment-row"><strong>' + t('primer') + '</strong><span>' + escapeHtml(query) + '</span></div>' +
          '<div class="alignment-row alignment-matchline"><strong></strong><span>' + escapeHtml(matchLine) + '</span></div>' +
          '<div class="alignment-row"><strong>' + t('sample') + '</strong><span>' + escapeHtml(subject) + '</span></div>' +
          '</div><p class="section-note">' + t('alignment_orientation') + '</p>' +
          '<p class="section-note">' + t('mismatch_details') + ': ' + escapeHtml(row.Mismatch_Details || t('none')) + '</p>';
      }
      modal.classList.add('open');
    }
    function closeAlignmentModal() {
      document.getElementById('alignment-modal').classList.remove('open');
    }
    function renderDetailTable(primer, group) {
      const state = primerState(primer);
      const sortedRows = sortedDetailRows(group, state);
      const visibleRows = state.expanded ? sortedRows : sortedRows.slice(0, 10);
      const tableId = primerTableId(primer);
      const hiddenCount = Math.max(0, sortedRows.length - visibleRows.length);
      const rowsHtml = visibleRows.map(row => '<tr>' +
        detailColumns.map(col => {
          if (col.key === 'Hit_Status') return '<td><span class="badge ' + escapeHtml(row[col.key]) + '">' + escapeHtml(statusLabel(row)) + '</span></td>';
          if (col.key === 'Subject_Sequence_ID') return '<td><button class="alignment-button" type="button" data-alignment-key="' + escapeHtml(rowKey(row)) + '" title="' + escapeHtml(t('open_alignment', { sample: row[col.key] || t('not_available') })) + '">' + escapeHtml(displayValue(row[col.key])) + '</button></td>';
          if (col.key === 'Percent_Identity') return '<td>' + formatIdentityValue(row[col.key]) + '</td>';
          return '<td>' + escapeHtml(displayValue(row[col.key])) + '</td>';
        }).join('') +
        '</tr>').join('');
      return '<div class="report-note report-note-compact">' +
          '<div class="note-row"><strong>' + t('detail_table_note_label') + '</strong><span>' + t('sample_rows_sorted', { visible: visibleRows.length, total: sortedRows.length }) + '</span></div>' +
          (sortedRows.length > 10 ? '<button class="table-toggle" type="button" data-table-toggle="' + escapeHtml(primer) + '">' + (state.expanded ? t('show_top_10') : t('show_all', { total: sortedRows.length })) + '</button>' : '') +
        '</div>' +
        '<div class="table-wrap" style="margin-top:12px"><table class="primer-detail-table" id="' + tableId + '"><thead><tr>' +
          detailColumns.map(col => '<th data-primer="' + escapeHtml(primer) + '" data-detail-sort="' + col.key + '">' + escapeHtml(t(col.labelKey)) + (state.sortKey === col.key ? (state.direction > 0 ? ' ▲' : ' ▼') : '') + '</th>').join('') +
          '</tr></thead><tbody>' + rowsHtml + '</tbody></table></div>' +
        (hiddenCount ? '<div class="review-reason">' + t('lower_priority_hidden', { count: hiddenCount }) + '</div>' : '');
    }
    function renderSummary(data) {
      const groups = summaryRows(data);
      groups.sort((a, b) => {
        const left = a[sortState.key];
        const right = b[sortState.key];
        if (sortState.key === 'Risk') return (a.Risk_Rank - b.Risk_Rank) * sortState.direction;
        const leftNum = Number(left);
        const rightNum = Number(right);
        if (Number.isFinite(leftNum) && Number.isFinite(rightNum)) return (leftNum - rightNum) * sortState.direction;
        return String(left).localeCompare(String(right)) * sortState.direction;
      });
      const columns = ['Risk', 'Primer_Name', 'Virus_Type', 'Primer_Segment', 'Rows', 'Hits', 'No_Hits', 'No_Hit_Rate', 'Max_Mismatches', 'Avg_Percent_Identity', 'TwoPlus_Rate', 'ThreePlus_Rate', 'Terminal_Mismatch_Share'];
      function summaryCellHtml(row, col) {
        if (col === 'Risk') return riskBadge(row.Risk, row.Risk_Title);
        if (col === 'Primer_Name') {
          return '<a class="primer-link" href="#' + primerPanelId(row.Primer_Name) + '" title="' + escapeHtml(t('jump_to_primer')) + '">' + escapeHtml(row.Primer_Name) + '</a>';
        }
        return escapeHtml(row[col]);
      }
      document.getElementById('summary-table').innerHTML =
        '<thead><tr>' + columns.map(col => '<th data-sort="' + col + '">' + t('col_' + col) + '</th>').join('') + '</tr></thead>' +
        '<tbody>' + groups.map(row => '<tr class="risk-' + row.Risk.toLowerCase() + '-row">' + columns.map(col => '<td>' + summaryCellHtml(row, col) + '</td>').join('') + '</tr>').join('') + '</tbody>';
      document.querySelectorAll('#summary-table th').forEach(th => th.addEventListener('click', () => {
        const key = th.dataset.sort;
        sortState.direction = sortState.key === key ? sortState.direction * -1 : 1;
        sortState.key = key;
        render();
      }));
    }
    function renderPrimerPanels(data) {
      const panels = groupByPrimer(data).map(([primer, group]) => {
        const stats = calculatePrimerStats(group);
        const risk = calculateRisk(stats);
        const hitPercent = group.length ? Math.round((stats.hits / group.length) * 100) : 0;
        return '<article class="primer-panel risk-' + risk.toLowerCase() + '-panel" id="' + primerPanelId(primer) + '"><h3>' + escapeHtml(primer) + ' ' + riskBadge(risk, riskExplanation(stats)) + '</h3>' +
          '<div class="primer-grid">' +
          '<div>' + t('rows') + '<strong><br>' + group.length + '</strong></div>' +
          '<div>' + t('hits') + '<strong><br>' + stats.hits + '</strong></div>' +
          '<div>' + t('no_hits') + '<strong><br>' + stats.noHits + '</strong></div>' +
          '<div>' + t('max_mismatches') + '<strong><br>' + stats.maxMismatches + '</strong></div>' +
          '<div>' + t('avg_identity') + '<strong><br>' + (stats.avgPercentIdentity === null ? t('not_available') : stats.avgPercentIdentity.toFixed(2)) + '</strong></div>' +
          '<div>' + t('two_plus_rate') + '<strong><br>' + formatPercent(stats.twoPlusMismatchRate) + '</strong></div>' +
          '<div>' + t('three_plus_rate') + '<strong><br>' + formatPercent(stats.threePlusMismatchRate) + '</strong></div>' +
          '<div>' + t(terminalLabelKey(stats.terminalRole)) + '<strong><br>' + formatPercent(stats.terminalMismatchRate) + '</strong></div>' +
          '</div><div class="section-note">' + t('hit_rate_label', { percent: hitPercent + '%' }) + '</div><div class="bar" aria-label="' + escapeHtml(t('hit_rate_aria', { percent: hitPercent + '%' })) + '"><span style="width:' + hitPercent + '%"></span></div>' +
          renderSequenceMap(group) +
          '<div class="plot-card"><h3>' + t('mismatch_distribution_title') + '</h3>' +
          '<div class="report-note report-note-compact"><div class="note-row"><strong>' + t('distribution_note_label') + '</strong><span>' + t('distribution_note') + '</span></div></div>' +
          distributionChart(stats.distribution) + '</div>' +
          renderDetailTable(primer, group) +
          '</article>';
      }).join('');
      document.getElementById('primer-panels').innerHTML = panels || '<div class="empty">' + t('no_rows_match') + '</div>';
      document.querySelectorAll('[data-detail-sort]').forEach(th => th.addEventListener('click', () => {
        const primer = th.dataset.primer;
        const state = primerState(primer);
        const key = th.dataset.detailSort;
        state.direction = state.sortKey === key ? state.direction * -1 : (key === 'Mismatches' ? -1 : 1);
        state.sortKey = key;
        render();
      }));
      document.querySelectorAll('[data-table-toggle]').forEach(button => button.addEventListener('click', () => {
        const state = primerState(button.dataset.tableToggle);
        state.expanded = !state.expanded;
        render();
      }));
      const rowByKey = new Map(data.map(row => [rowKey(row), row]));
      document.querySelectorAll('[data-alignment-key]').forEach(button => button.addEventListener('click', () => {
        const row = rowByKey.get(button.dataset.alignmentKey);
        if (row) showAlignmentModal(row);
      }));
    }
    function render() {
      const data = currentRows();
      renderCards(data);
      renderNgsPanelOverview(data);
      renderRiskLegend();
      renderSummary(data);
      renderPrimerPanels(data);
    }
    function populateFilterOptions() {
      fillSelect(filters.virus, uniqueValues('Virus_Type'), 'all_organisms');
      fillSelect(filters.primer, uniqueValues('Primer_Name'), 'all_primers');
      fillSelect(filters.fasta, uniqueValues('Fasta_File'), 'all_fasta');
      fillRiskSelect();
    }
    populateFilterOptions();
    Object.values(filters).forEach(control => control.addEventListener('input', render));
    document.querySelectorAll('[data-language]').forEach(button => button.addEventListener('click', () => {
      currentLanguage = button.dataset.language;
      const selected = { virus: filters.virus.value, primer: filters.primer.value, fasta: filters.fasta.value, risk: filters.risk.value };
      applyTranslations();
      populateFilterOptions();
      filters.virus.value = selected.virus;
      filters.primer.value = selected.primer;
      filters.fasta.value = selected.fasta;
      filters.risk.value = selected.risk;
      render();
    }));
    document.querySelectorAll('[data-tab]').forEach(button => button.addEventListener('click', () => {
      const selected = button.dataset.tab;
      document.querySelectorAll('[data-tab]').forEach(tab => tab.classList.toggle('active', tab.dataset.tab === selected));
      document.querySelectorAll('.tab-panel').forEach(panel => {
        panel.hidden = panel.id !== selected + '-tab';
      });
    }));
    document.getElementById('alignment-modal-close').addEventListener('click', closeAlignmentModal);
    document.getElementById('alignment-modal').addEventListener('click', event => {
      if (event.target.id === 'alignment-modal') closeAlignmentModal();
    });
    document.addEventListener('keydown', event => {
      if (event.key === 'Escape') closeAlignmentModal();
    });
    applyTranslations();
    render();
  </script>
</body>
</html>
"""
    return (
        template
        .replace("__REPORT_TITLE__", escaped_title)
        .replace("__REPORT_DATA__", data_json)
        .replace("__TRANSLATIONS_JSON__", translations_json)
        .replace("__FIELDS_JSON__", fields_json)
    )


def write_html_report(results: list[dict], output_file: str, previous_reports: list[dict] | None = None):
    report_html = build_html_report(results, previous_reports=previous_reports)
    with open(output_file, "w", encoding="utf-8") as htmlfile:
        htmlfile.write(report_html)
    print(f"HTML report successfully written to {output_file}")
