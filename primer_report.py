#!/usr/bin/env python3
"""HTML report generation for primer checker result rows."""

import csv
import html
import json
import os
import sys

from primer_analysis import CSV_FIELDNAMES

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
                report["rows"].append({field: row.get(field, "") for field in CSV_FIELDNAMES})
    except Exception as e:
        report["warnings"].append(f"Could not read CSV: {e}")
    return report


def load_previous_report_csvs(paths: list[str] | None) -> list[dict]:
    return [load_previous_report_csv(path) for path in paths or []]


def build_html_report(
    results: list[dict],
    title: str = "Primer Checker Report",
    previous_reports: list[dict] | None = None,
) -> str:
    """Build a self-contained HTML report with embedded result data."""
    report_rows = []
    for row in results:
        report_rows.append({field: row.get(field, "") for field in CSV_FIELDNAMES})

    previous_reports = previous_reports or []
    data_json = (
        json.dumps(report_rows, ensure_ascii=True)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )
    previous_json = (
        json.dumps(previous_reports, ensure_ascii=True)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )
    fields_json = json.dumps(CSV_FIELDNAMES)
    escaped_title = html.escape(title)
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
    .card, .plot-card, .previous-report-panel {
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
    .primer-panel {
      margin-top: 12px;
      padding: 14px;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
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
    .primer-sequence-map {
      width: 100%;
      overflow-x: auto;
      padding-bottom: 4px;
    }
    .timeline-plot {
      width: 100%;
      overflow-x: auto;
      padding-bottom: 4px;
    }
    .timeline-svg {
      max-width: 100%;
      min-width: 560px;
      height: auto;
      display: block;
    }
    .timeline-mismatch { stroke: var(--bad); fill: none; stroke-width: 2.5; }
    .timeline-ct { stroke: var(--accent); fill: none; stroke-width: 2.5; stroke-dasharray: 5 3; }
    .timeline-point-mismatch { fill: var(--bad); }
    .timeline-point-ct { fill: var(--accent); }
    .timeline-axis { stroke: #9aa6b2; }
    .timeline-label { fill: var(--muted); font-size: 11px; }
    .timeline-legend {
      display: flex;
      gap: 14px;
      flex-wrap: wrap;
      color: var(--muted);
      font-size: 13px;
      margin-top: 8px;
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
    .previous-reports {
      display: grid;
      gap: 12px;
    }
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
    .empty { padding: 18px; color: var(--muted); }
    @media (max-width: 700px) {
      header, main { padding-left: 16px; padding-right: 16px; }
      .distribution-row { grid-template-columns: 52px 1fr; }
      .distribution-row span:last-child { grid-column: 2; }
    }
  </style>
</head>
<body>
  <header>
    <h1>__REPORT_TITLE__</h1>
    <p class="meta">Static local report generated from primer checker result rows.</p>
  </header>
  <main>
    <section class="filters" aria-label="Report filters">
      <label>Organism or virus<select id="filter-virus"></select></label>
      <label>Primer<select id="filter-primer"></select></label>
      <label>FASTA file<select id="filter-fasta"></select></label>
      <label>Hit status<select id="filter-status"></select></label>
      <label>Maximum mismatches<input id="filter-mismatches" type="number" min="0" step="1" placeholder="Any"></label>
    </section>
    <section class="cards" id="overview-cards"></section>
    <section id="ngs-panel-overview"></section>
    <h2>Primer Overview</h2>
    <p class="section-note">Risk combines no-hit rate, average mismatches, mismatch burden, and terminal mismatch concentration near the 5' and 3' ends.</p>
    <div class="legend" id="risk-legend"></div>
    <div class="table-wrap"><table id="summary-table"></table></div>
    <h2>Attached Previous Reports</h2>
    <p class="section-note">Optional CSV reports attached at generation time for local comparison.</p>
    <section class="previous-reports" id="previous-reports"></section>
    <h2>Primer Investigation</h2>
    <section id="primer-panels"></section>
  </main>
  <div class="modal-backdrop" id="alignment-modal" role="dialog" aria-modal="true" aria-labelledby="alignment-modal-title">
    <div class="modal">
      <div class="modal-header">
        <div>
          <h2 id="alignment-modal-title">Primer alignment</h2>
          <p class="section-note" id="alignment-modal-meta"></p>
        </div>
        <button class="modal-close" type="button" id="alignment-modal-close">Close</button>
      </div>
      <div id="alignment-modal-body"></div>
    </div>
  </div>
  <script id="report-data" type="application/json">__REPORT_DATA__</script>
  <script id="previous-report-data" type="application/json">__PREVIOUS_REPORT_DATA__</script>
  <script>
    const rows = JSON.parse(document.getElementById('report-data').textContent);
    const previousReports = JSON.parse(document.getElementById('previous-report-data').textContent);
    const fields = __FIELDS_JSON__;
    const RISK_THRESHOLDS = {
      watchAvgMismatch: 1,
      highAvgMismatch: 2,
      criticalAvgMismatch: 3,
      highNoHitRate: 0.05,
      criticalNoHitRate: 0.20,
      highTwoPlusMismatchRate: 0.30,
      criticalThreePlusMismatchRate: 0.25,
      highTerminalMismatchRate: 0.25,
      criticalTerminalMismatchRate: 0.50,
      terminalBases: 5
    };
    const riskRank = { Low: 0, Watch: 1, High: 2, Critical: 3 };
    const distCategories = [
      { key: '0', label: '0 mismatches' },
      { key: '1', label: '1 mismatch' },
      { key: '2', label: '2 mismatches' },
      { key: '3', label: '3 mismatches' },
      { key: '4plus', label: '4+ mismatches' },
      { key: 'no_hit', label: 'No hit' }
    ];
    const filters = {
      virus: document.getElementById('filter-virus'),
      primer: document.getElementById('filter-primer'),
      fasta: document.getElementById('filter-fasta'),
      status: document.getElementById('filter-status'),
      mismatches: document.getElementById('filter-mismatches')
    };
    let sortState = { key: 'Risk', direction: -1 };
    const primerTableState = {};
    const detailColumns = [
      { key: 'Fasta_File', label: 'FASTA file' },
      { key: 'Virus_Type', label: 'Virus' },
      { key: 'Subject_Sequence_ID', label: 'Sample' },
      { key: 'Subject_Segment', label: 'Segment' },
      { key: 'Hit_Status', label: 'Status' },
      { key: 'Percent_Identity', label: 'Identity' },
      { key: 'Mismatches', label: 'Mismatches' },
      { key: 'Mismatch_Positions', label: 'Positions' },
      { key: 'Mismatch_Details', label: 'Mismatch bases' },
      { key: 'Metadata_Sample_ID', label: 'Metadata ID' },
      { key: 'Sample_Date', label: 'Sample date' },
      { key: 'Ct_Value', label: 'Ct' },
      { key: 'Ct_Source', label: 'Ct source' }
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
    function fillSelect(select, values, label) {
      select.innerHTML = '<option value="">All ' + label + '</option>' + values.map(value => '<option>' + escapeHtml(value) + '</option>').join('');
    }
    function escapeHtml(value) {
      return String(value).replace(/[&<>"']/g, character => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[character]));
    }
    function formatPercent(value) {
      return Number.isFinite(value) ? (value * 100).toFixed(1) + '%' : 'n/a';
    }
    function parseDateValue(value) {
      if (!value) return null;
      const text = String(value).trim();
      const dmy = text.match(/^(\\d{1,2})[./-](\\d{1,2})[./-](\\d{4})$/);
      const parsed = dmy ? new Date(Number(dmy[3]), Number(dmy[2]) - 1, Number(dmy[1])) : new Date(text);
      return Number.isFinite(parsed.getTime()) ? parsed : null;
    }
    function formatDateLabel(date) {
      return date.toISOString().slice(0, 10);
    }
    function currentRows() {
      const maxMismatches = filters.mismatches.value === '' ? null : Number(filters.mismatches.value);
      return rows.filter(row => {
        if (filters.virus.value && row.Virus_Type !== filters.virus.value) return false;
        if (filters.primer.value && row.Primer_Name !== filters.primer.value) return false;
        if (filters.fasta.value && row.Fasta_File !== filters.fasta.value) return false;
        if (filters.status.value && row.Hit_Status !== filters.status.value) return false;
        const mismatches = numeric(row.Mismatches);
        if (maxMismatches !== null && (mismatches === null || mismatches > maxMismatches)) return false;
        return true;
      });
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
    function calculateMismatchDistribution(group) {
      const distribution = { '0': 0, '1': 0, '2': 0, '3': 0, '4plus': 0, no_hit: 0 };
      for (const row of group) {
        if (row.Hit_Status === 'no_hit' || row.Percent_Identity === 'No hit') {
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
      const hits = group.filter(row => row.Hit_Status !== 'no_hit');
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
      const hitRows = group.filter(row => row.Hit_Status !== 'no_hit' && row.Percent_Identity !== 'No hit');
      const noHits = totalRows - hitRows.length;
      const mismatchValues = hitRows.map(row => parseNumber(row.Mismatches)).filter(value => value !== null);
      const identityValues = hitRows.map(row => parseNumber(row.Percent_Identity)).filter(value => value !== null);
      const avgMismatches = mismatchValues.length ? mismatchValues.reduce((a, b) => a + b, 0) / mismatchValues.length : 0;
      const maxMismatches = mismatchValues.length ? Math.max(...mismatchValues) : 0;
      const avgPercentIdentity = identityValues.length ? identityValues.reduce((a, b) => a + b, 0) / identityValues.length : null;
      const twoPlus = mismatchValues.filter(value => value >= 2).length;
      const threePlus = mismatchValues.filter(value => value >= 3).length;
      const positionData = calculateMismatchPositionCounts(group);
      const terminalBases = RISK_THRESHOLDS.terminalBases;
      const terminalMismatchEvents = positionData.counts.reduce((sum, count, index) => {
        const pos = index + 1;
        return sum + (pos <= terminalBases || pos > positionData.counts.length - terminalBases ? count : 0);
      }, 0);
      const totalMismatchEvents = positionData.counts.reduce((sum, count) => sum + count, 0);
      return {
        totalRows,
        hits: hitRows.length,
        noHits,
        noHitRate: totalRows ? noHits / totalRows : 0,
        avgMismatches,
        maxMismatches,
        avgPercentIdentity,
        twoPlusMismatchRate: hitRows.length ? twoPlus / hitRows.length : 0,
        threePlusMismatchRate: hitRows.length ? threePlus / hitRows.length : 0,
        terminalMismatchRate: totalMismatchEvents ? terminalMismatchEvents / totalMismatchEvents : 0,
        distribution: calculateMismatchDistribution(group)
      };
    }
    function calculateRisk(stats) {
      if (
        stats.noHitRate >= RISK_THRESHOLDS.criticalNoHitRate ||
        stats.avgMismatches >= RISK_THRESHOLDS.criticalAvgMismatch ||
        stats.threePlusMismatchRate >= RISK_THRESHOLDS.criticalThreePlusMismatchRate ||
        stats.terminalMismatchRate >= RISK_THRESHOLDS.criticalTerminalMismatchRate
      ) return 'Critical';
      if (
        stats.noHitRate >= RISK_THRESHOLDS.highNoHitRate ||
        stats.avgMismatches >= RISK_THRESHOLDS.highAvgMismatch ||
        stats.twoPlusMismatchRate >= RISK_THRESHOLDS.highTwoPlusMismatchRate ||
        stats.terminalMismatchRate >= RISK_THRESHOLDS.highTerminalMismatchRate
      ) return 'High';
      if (stats.avgMismatches >= RISK_THRESHOLDS.watchAvgMismatch || stats.maxMismatches >= 2 || stats.noHits > 0) return 'Watch';
      return 'Low';
    }
    function riskExplanation(stats) {
      return [
        'No-hit rate: ' + formatPercent(stats.noHitRate),
        'Average mismatches: ' + stats.avgMismatches.toFixed(2),
        'Max mismatches: ' + stats.maxMismatches,
        '2+ mismatch rate: ' + formatPercent(stats.twoPlusMismatchRate),
        '3+ mismatch rate: ' + formatPercent(stats.threePlusMismatchRate),
        "Terminal mismatch share (5' or 3' end): " + formatPercent(stats.terminalMismatchRate)
      ].join(' | ');
    }
    function riskBadge(risk, title) {
      return '<span class="risk-badge risk-' + risk.toLowerCase() + '" title="' + escapeHtml(title) + '">' + risk + '</span>';
    }
    function renderCards(data) {
      const hits = data.filter(row => row.Hit_Status === 'hit').length;
      const noHits = data.filter(row => row.Hit_Status === 'no_hit').length;
      const primers = new Set(data.map(row => row.Primer_Name)).size;
      const samples = new Set(data.map(row => row.Subject_Sequence_ID)).size;
      const metadataMatches = new Set(data.map(row => row.Metadata_Sample_ID).filter(Boolean)).size;
      const avgIdentity = average(data.map(row => row.Percent_Identity)) || 'n/a';
      document.getElementById('overview-cards').innerHTML = [
        ['Rows', data.length],
        ['Primers', primers],
        ['Samples', samples],
        ['Metadata matches', metadataMatches],
        ['Hits', hits],
        ['No hits', noHits],
        ['Avg identity', avgIdentity]
      ].map(([label, value]) => '<div class="card"><h3>' + label + '</h3><div class="metric">' + value + '</div></div>').join('');
    }
    function renderRiskLegend() {
      document.getElementById('risk-legend').innerHTML =
        '<span>Risk:</span>' +
        riskBadge('Low', 'Low risk under current thresholds') +
        riskBadge('Watch', 'Some no-hit or mismatch signal') +
        riskBadge('High', 'Elevated no-hit or mismatch signal') +
        riskBadge('Critical', 'Strong no-hit, mismatch, or terminal mismatch signal');
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
        if (!grouped.has(key)) grouped.set(key, { amplicon, assayId: row.Assay_ID || '', assayName: row.Assay_Name || row.Assay_ID || 'Unnamed NGS panel', rows: [], primers: new Set(), pools: new Set() });
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
          : 'Critical';
        const stats = calculatePrimerStats(group.rows);
        const directionStatus = 'Forward/LEFT: ' + (leftViable ? 'viable (' + bestLeft + ')' : 'NO viable primer') +
          ' | Reverse/RIGHT: ' + (rightViable ? 'viable (' + bestRight + ')' : 'NO viable primer');
        const title = group.assayName + ' | ' + group.amplicon + ' | ' + directionStatus + ' | Primers: ' + [...group.primers].join(', ');
        return '<div class="amplicon-tile risk-' + risk.toLowerCase() + '" title="' + escapeHtml(title) + '">' +
          '<strong>' + escapeHtml(group.amplicon) + '</strong>' +
          riskBadge(risk, title) +
          '<div>' + group.primers.size + ' primer' + (group.primers.size === 1 ? '' : 's') +
          (group.pools.size ? ' · pool ' + escapeHtml([...group.pools].join(', ')) : '') + '</div>' +
          '<div>L: ' + (leftViable ? bestLeft : 'missing') + ' · R: ' + (rightViable ? bestRight : 'missing') + '</div></div>';
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
      target.innerHTML = '<h2>NGS Amplicon Panel Overview</h2>' +
        '<p class="section-note">Each panel is shown separately. Alternative primers are redundant: an amplicon is critical only when it has no viable LEFT/forward primer or no viable RIGHT/reverse primer. Low and Watch primers are considered viable. Hover for details; filters also apply.</p>' +
        '<div class="legend">' + riskBadge('Low', 'No current warning signal') + riskBadge('Watch', 'Review') + riskBadge('High', 'Elevated risk') + riskBadge('Critical', 'Strong risk signal') + '</div>' +
        panels;
    }
    function distributionChart(distribution) {
      const total = Object.values(distribution).reduce((sum, value) => sum + value, 0);
      if (!total) return '<div class="empty">No mismatch-count data available.</div>';
      return '<div class="mismatch-distribution">' + distCategories.map(category => {
        const count = distribution[category.key] || 0;
        const pct = total ? count / total : 0;
        return '<div class="distribution-row" title="' + escapeHtml(category.label + ': ' + count + ' rows, ' + formatPercent(pct)) + '">' +
          '<strong>' + escapeHtml(category.label) + '</strong>' +
          '<div class="distribution-track"><div class="distribution-fill dist-' + category.key + '" style="width:' + Math.max(0, pct * 100) + '%"></div></div>' +
          '<span>' + count + ' (' + formatPercent(pct) + ')</span>' +
          '</div>';
      }).join('') + '</div>';
    }
    function renderSequenceMap(group) {
      const { primerSequence, counts, detailsByPosition, hitRows } = calculateMismatchPositionCounts(group);
      if (!primerSequence || !counts.length || !hitRows || counts.every(count => count === 0)) {
        return '<div class="plot-card"><h3>Mismatch distribution along primer sequence</h3><div class="empty">No mismatch-position data available.</div></div>';
      }
      const width = Math.max(560, counts.length * 26 + 96);
      const height = 142;
      const plotHeight = 58;
      const left = 72;
      const right = 18;
      const cellWidth = (width - left - right) / counts.length;
      const maxCount = Math.max(...counts, 1);
      const terminalBases = RISK_THRESHOLDS.terminalBases;
      const yTicks = [0, Math.ceil(maxCount / 2), maxCount].filter((value, index, arr) => arr.indexOf(value) === index);
      const tickMarks = yTicks.map(value => {
        const y = 78 - (value / maxCount) * plotHeight;
        return '<g><line x1="' + left + '" y1="' + y.toFixed(1) + '" x2="' + (width - right) + '" y2="' + y.toFixed(1) + '" stroke="#e2e8f0"></line>' +
          '<text class="axis-label" x="' + (left - 8) + '" y="' + (y + 3).toFixed(1) + '" text-anchor="end">' + value + '</text></g>';
      }).join('');
      const cells = counts.map((count, index) => {
        const pos = index + 1;
        const base = primerSequence[index] || '';
        const changes = topMismatchChanges(detailsByPosition.get(pos));
        const pct = hitRows ? count / hitRows : 0;
        const barHeight = count ? Math.max(4, (count / maxCount) * plotHeight) : 1;
        const x = left + index * cellWidth;
        const y = 78 - barHeight;
        const terminal = pos <= terminalBases || pos > counts.length - terminalBases;
        const opacity = count ? Math.min(1, 0.22 + (count / maxCount) * 0.78) : 0.12;
        return '<g><title>Position ' + pos + ', base ' + escapeHtml(base) + ': ' + count + ' mismatches (' + formatPercent(pct) + ' of hit rows)' + (changes ? '; ' + escapeHtml(changes) : '') + '</title>' +
          '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + Math.max(3, cellWidth - 4).toFixed(1) + '" height="' + barHeight.toFixed(1) + '" fill="' + (terminal ? '#b42318' : '#1f7a8c') + '" opacity="' + opacity.toFixed(2) + '"></rect>' +
          '<text class="base-label" x="' + (x + cellWidth / 2).toFixed(1) + '" y="102" text-anchor="middle">' + escapeHtml(base) + '</text>' +
          (pos === 1 || pos === counts.length || pos % 5 === 0 ? '<text class="axis-label" x="' + (x + cellWidth / 2).toFixed(1) + '" y="124" text-anchor="middle">' + pos + '</text>' : '') +
          '</g>';
      }).join('');
      return '<div class="plot-card"><h3>Mismatch distribution along primer sequence</h3>' +
        '<p class="section-note">Bars show mismatch frequency by primer position across hit rows. Hover over bars to see nucleotide-change details. The first and last ' + terminalBases + " bases mark 5' and 3' terminal regions.</p>" +
        '<div class="primer-sequence-map"><svg class="sequence-svg" viewBox="0 0 ' + width + ' ' + height + '" role="img" aria-label="Mismatch distribution along primer sequence">' +
        '<text class="axis-title" transform="translate(14 50) rotate(-90)" text-anchor="middle">count</text>' +
        tickMarks +
        '<text class="five-prime" x="' + (left - 25) + '" y="102">5&apos;</text><text class="three-prime" x="' + (width - 18) + '" y="102">3&apos;</text>' +
        '<line x1="' + left + '" y1="82" x2="' + (width - right) + '" y2="82" stroke="#9aa6b2"></line>' +
        cells +
        '</svg></div></div>';
    }
    function calculatePrimerTimeline(group) {
      const byDate = new Map();
      for (const row of group) {
        const date = parseDateValue(row.Sample_Date);
        if (!date) continue;
        const ctSource = row.Ct_Source || 'Ct_Value';
        const key = formatDateLabel(date) + '|' + ctSource;
        if (!byDate.has(key)) byDate.set(key, { date, ctSource, rows: 0, mismatches: [], cts: [], noHits: 0 });
        const bucket = byDate.get(key);
        bucket.rows += 1;
        if (row.Hit_Status === 'no_hit' || row.Percent_Identity === 'No hit') bucket.noHits += 1;
        const mismatches = parseNumber(row.Mismatches);
        if (mismatches !== null) bucket.mismatches.push(mismatches);
        const ct = parseNumber(row.Ct_Value);
        if (ct !== null) bucket.cts.push(ct);
      }
      return [...byDate.values()].sort((a, b) => a.date - b.date).map(bucket => ({
        date: bucket.date,
        label: formatDateLabel(bucket.date),
        ctSource: bucket.ctSource,
        rows: bucket.rows,
        noHits: bucket.noHits,
        avgMismatches: bucket.mismatches.length ? bucket.mismatches.reduce((a, b) => a + b, 0) / bucket.mismatches.length : null,
        avgCt: bucket.cts.length ? bucket.cts.reduce((a, b) => a + b, 0) / bucket.cts.length : null
      }));
    }
    function renderTimeline(group) {
      const points = calculatePrimerTimeline(group);
      if (!points.length) {
        return '<div class="plot-card"><h3>Sample timeline</h3><div class="empty">No sample-date metadata available for this primer.</div></div>';
      }
      const mismatchValues = points.map(point => point.avgMismatches).filter(value => value !== null);
      const ctValues = points.map(point => point.avgCt).filter(value => value !== null);
      if (!mismatchValues.length && !ctValues.length) {
        return '<div class="plot-card"><h3>Sample timeline</h3><div class="empty">Sample dates are available, but no numeric mismatch or Ct values can be plotted.</div></div>';
      }
      const width = Math.max(560, points.length * 72 + 90);
      const height = 230;
      const left = 48;
      const right = 54;
      const top = 24;
      const bottom = 50;
      const plotWidth = width - left - right;
      const plotHeight = height - top - bottom;
      const minTime = Math.min(...points.map(point => point.date.getTime()));
      const maxTime = Math.max(...points.map(point => point.date.getTime()));
      const timeRange = Math.max(1, maxTime - minTime);
      const maxMismatch = Math.max(1, ...mismatchValues, 1);
      const minCt = ctValues.length ? Math.min(...ctValues) : 0;
      const maxCt = ctValues.length ? Math.max(...ctValues) : 1;
      const ctRange = Math.max(1, maxCt - minCt);
      const x = point => left + ((point.date.getTime() - minTime) / timeRange) * plotWidth;
      const yMismatch = value => top + plotHeight - (value / maxMismatch) * plotHeight;
      const yCt = value => top + plotHeight - ((value - minCt) / ctRange) * plotHeight;
      const mismatchTicks = [0, maxMismatch / 2, maxMismatch];
      const ctTicks = ctValues.length ? [minCt, minCt + ctRange / 2, maxCt] : [];
      const yAxisTicks = mismatchTicks.map(value => {
        const y = yMismatch(value);
        return '<g><line x1="' + left + '" y1="' + y.toFixed(1) + '" x2="' + (width - right) + '" y2="' + y.toFixed(1) + '" stroke="#e2e8f0"></line>' +
          '<text class="timeline-label" x="' + (left - 8) + '" y="' + (y + 4).toFixed(1) + '" text-anchor="end">' + value.toFixed(maxMismatch < 2 ? 1 : 0) + '</text></g>';
      }).join('') + ctTicks.map(value => {
        const y = yCt(value);
        return '<text class="timeline-label" x="' + (width - right + 8) + '" y="' + (y + 4).toFixed(1) + '" text-anchor="start">' + value.toFixed(1) + '</text>';
      }).join('');
      const mismatchPolyline = points
        .filter(point => point.avgMismatches !== null)
        .map(point => x(point).toFixed(1) + ',' + yMismatch(point.avgMismatches).toFixed(1))
        .join(' ');
      const ctPolyline = points
        .filter(point => point.avgCt !== null)
        .map(point => x(point).toFixed(1) + ',' + yCt(point.avgCt).toFixed(1))
        .join(' ');
      const labels = points.map((point, index) => {
        if (index !== 0 && index !== points.length - 1 && index % Math.ceil(points.length / 6) !== 0) return '';
        return '<text class="timeline-label" x="' + x(point).toFixed(1) + '" y="' + (height - 16) + '" text-anchor="middle">' + escapeHtml(point.label) + '</text>';
      }).join('');
      const mismatchDots = points.filter(point => point.avgMismatches !== null).map(point =>
        '<circle class="timeline-point-mismatch" cx="' + x(point).toFixed(1) + '" cy="' + yMismatch(point.avgMismatches).toFixed(1) + '" r="4"><title>' +
        escapeHtml(point.label + ': avg mismatches ' + point.avgMismatches.toFixed(2) + ', rows ' + point.rows + ', no-hit rows ' + point.noHits) +
        '</title></circle>'
      ).join('');
      const ctDots = points.filter(point => point.avgCt !== null).map(point =>
        '<circle class="timeline-point-ct" cx="' + x(point).toFixed(1) + '" cy="' + yCt(point.avgCt).toFixed(1) + '" r="4"><title>' +
        escapeHtml(point.label + ': avg Ct ' + point.avgCt.toFixed(2) + ' from ' + point.ctSource + ', rows ' + point.rows) +
        '</title></circle>'
      ).join('');
      const ctSources = [...new Set(points.map(point => point.ctSource).filter(Boolean))].join(', ');
      return '<div class="plot-card"><h3>Sample timeline</h3>' +
        '<p class="section-note">Average mismatches and selected assay-specific Ct values by sample date for the current filter. Use this to check whether primer mismatches accumulate over time and whether Ct values remain stable.' + (ctSources ? ' Ct source: ' + escapeHtml(ctSources) + '.' : '') + '</p>' +
        '<div class="timeline-plot"><svg class="timeline-svg" viewBox="0 0 ' + width + ' ' + height + '" role="img" aria-label="Primer mismatch and Ct timeline">' +
        '<line class="timeline-axis" x1="' + left + '" y1="' + (top + plotHeight) + '" x2="' + (width - right) + '" y2="' + (top + plotHeight) + '"></line>' +
        '<line class="timeline-axis" x1="' + left + '" y1="' + top + '" x2="' + left + '" y2="' + (top + plotHeight) + '"></line>' +
        '<line class="timeline-axis" x1="' + (width - right) + '" y1="' + top + '" x2="' + (width - right) + '" y2="' + (top + plotHeight) + '"></line>' +
        '<text class="timeline-label" x="8" y="' + top + '">mismatch</text><text class="timeline-label" x="' + (width - 44) + '" y="' + top + '">Ct</text>' +
        yAxisTicks +
        (mismatchPolyline ? '<polyline class="timeline-mismatch" points="' + mismatchPolyline + '"></polyline>' : '') +
        (ctPolyline ? '<polyline class="timeline-ct" points="' + ctPolyline + '"></polyline>' : '') +
        mismatchDots + ctDots + labels +
        '</svg></div><div class="timeline-legend"><span><strong style="color:var(--bad)">solid</strong> avg mismatches</span><span><strong style="color:var(--accent)">dashed</strong> avg Ct</span></div></div>';
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
          Avg_Mismatches: stats.hits ? stats.avgMismatches.toFixed(2) : 'n/a',
          Max_Mismatches: stats.maxMismatches,
          Avg_Percent_Identity: stats.avgPercentIdentity === null ? 'n/a' : stats.avgPercentIdentity.toFixed(2),
          TwoPlus_Rate: formatPercent(stats.twoPlusMismatchRate),
          ThreePlus_Rate: formatPercent(stats.threePlusMismatchRate),
          Terminal_Mismatch_Share: formatPercent(stats.terminalMismatchRate)
        };
      });
    }
    function primerTableId(primer) {
      return 'primer-table-' + String(primer).replace(/[^A-Za-z0-9_-]/g, '_');
    }
    function primerState(primer) {
      if (!primerTableState[primer]) {
        primerTableState[primer] = { expanded: false, sortKey: 'Mismatches', direction: -1 };
      }
      return primerTableState[primer];
    }
    function valueForSort(row, key) {
      if (key === 'Mismatches') {
        if (row.Hit_Status === 'no_hit' || row.Percent_Identity === 'No hit') return Number.NEGATIVE_INFINITY;
        return parseNumber(row.Mismatches) ?? Number.NEGATIVE_INFINITY;
      }
      if (key === 'Percent_Identity') return parseNumber(row.Percent_Identity) ?? Number.NEGATIVE_INFINITY;
      return row[key] || '';
    }
    function sortedDetailRows(group, state) {
      return [...group].sort((a, b) => {
        if (state.sortKey !== 'Hit_Status') {
          const aNoHit = a.Hit_Status === 'no_hit' || a.Percent_Identity === 'No hit';
          const bNoHit = b.Hit_Status === 'no_hit' || b.Percent_Identity === 'No hit';
          if (aNoHit !== bNoHit) return aNoHit ? 1 : -1;
        }
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
      title.textContent = row.Subject_Sequence_ID || 'Sample alignment';
      meta.textContent = [row.Primer_Name, row.Fasta_File, row.Hit_Status, row.Mismatches ? row.Mismatches + ' mismatches' : '', row.Mismatch_Details || ''].filter(Boolean).join(' | ');
      if (row.Hit_Status === 'no_hit' || !subject) {
        body.innerHTML = '<div class="empty">No alignment is available for this row.</div>';
      } else {
        const matchLine = alignmentMatchLine(query, subject);
        body.innerHTML = '<div class="alignment-view">' +
          '<div class="alignment-row"><strong>Primer</strong><span>' + escapeHtml(query) + '</span></div>' +
          '<div class="alignment-row alignment-matchline"><strong></strong><span>' + escapeHtml(matchLine) + '</span></div>' +
          '<div class="alignment-row"><strong>Sample</strong><span>' + escapeHtml(subject) + '</span></div>' +
          '</div><p class="section-note">Mismatch details: ' + escapeHtml(row.Mismatch_Details || 'none') + '</p>';
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
          if (col.key === 'Hit_Status') return '<td><span class="badge ' + escapeHtml(row[col.key]) + '">' + escapeHtml(row[col.key]) + '</span></td>';
          if (col.key === 'Subject_Sequence_ID') return '<td><button class="alignment-button" type="button" data-alignment-key="' + escapeHtml(rowKey(row)) + '">' + escapeHtml(row[col.key] ?? '') + '</button></td>';
          return '<td>' + escapeHtml(row[col.key] ?? '') + '</td>';
        }).join('') +
        '</tr>').join('');
      return '<div class="table-actions">' +
          '<span>Sample rows sorted by highest mismatch count by default. Showing ' + visibleRows.length + ' of ' + sortedRows.length + ' rows.</span>' +
          (sortedRows.length > 10 ? '<button class="table-toggle" type="button" data-table-toggle="' + escapeHtml(primer) + '">' + (state.expanded ? 'Show top 10' : 'Show all ' + sortedRows.length) + '</button>' : '') +
        '</div>' +
        '<div class="table-wrap" style="margin-top:12px"><table class="primer-detail-table" id="' + tableId + '"><thead><tr>' +
          detailColumns.map(col => '<th data-primer="' + escapeHtml(primer) + '" data-detail-sort="' + col.key + '">' + escapeHtml(col.label) + (state.sortKey === col.key ? (state.direction > 0 ? ' ▲' : ' ▼') : '') + '</th>').join('') +
          '</tr></thead><tbody>' + rowsHtml + '</tbody></table></div>' +
        (hiddenCount ? '<div class="review-reason">' + hiddenCount + ' lower-priority rows are hidden. Expand to inspect all rows.</div>' : '');
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
      const columns = ['Risk', 'Primer_Name', 'Virus_Type', 'Primer_Segment', 'Rows', 'Hits', 'No_Hits', 'No_Hit_Rate', 'Avg_Mismatches', 'Max_Mismatches', 'Avg_Percent_Identity', 'TwoPlus_Rate', 'ThreePlus_Rate', 'Terminal_Mismatch_Share'];
      document.getElementById('summary-table').innerHTML =
        '<thead><tr>' + columns.map(col => '<th data-sort="' + col + '">' + col.replaceAll('_', ' ') + '</th>').join('') + '</tr></thead>' +
        '<tbody>' + groups.map(row => '<tr class="risk-' + row.Risk.toLowerCase() + '-row">' + columns.map(col => '<td>' + (col === 'Risk' ? riskBadge(row.Risk, row.Risk_Title) : escapeHtml(row[col])) + '</td>').join('') + '</tr>').join('') + '</tbody>';
      document.querySelectorAll('#summary-table th').forEach(th => th.addEventListener('click', () => {
        const key = th.dataset.sort;
        sortState.direction = sortState.key === key ? sortState.direction * -1 : 1;
        sortState.key = key;
        render();
      }));
    }
    function renderPreviousReports(data) {
      const target = document.getElementById('previous-reports');
      if (!previousReports.length) {
        target.innerHTML = '<div class="empty">No previous CSV reports were attached when this HTML was generated.</div>';
        return;
      }
      const currentByPrimer = new Map(summaryRows(data).map(row => [row.Primer_Name, row]));
      target.innerHTML = previousReports.map(report => {
        const reportRows = report.rows || [];
        const summary = summaryRows(reportRows);
        const totalRows = reportRows.length;
        const primers = new Set(reportRows.map(row => row.Primer_Name)).size;
        const noHits = reportRows.filter(row => row.Hit_Status === 'no_hit' || row.Percent_Identity === 'No hit').length;
        const avgMismatch = average(reportRows.map(row => row.Mismatches)) || 'n/a';
        const risky = summary.filter(row => riskRank[row.Risk] >= riskRank.Watch).slice(0, 5);
        const warnings = (report.warnings || []).map(warning => '<div class="empty">' + escapeHtml(warning) + '</div>').join('');
        const comparisonRows = summary.filter(row => currentByPrimer.has(row.Primer_Name)).slice(0, 8).map(previous => {
          const current = currentByPrimer.get(previous.Primer_Name);
          const currentMismatch = parseNumber(current.Avg_Mismatches);
          const previousMismatch = parseNumber(previous.Avg_Mismatches);
          const mismatchDelta = currentMismatch !== null && previousMismatch !== null ? currentMismatch - previousMismatch : null;
          return '<tr><td>' + escapeHtml(previous.Primer_Name) + '</td><td>' + escapeHtml(previous.Avg_Mismatches) + '</td><td>' + escapeHtml(current.Avg_Mismatches) + '</td><td class="' + (mismatchDelta > 0 ? 'delta-up' : mismatchDelta < 0 ? 'delta-down' : '') + '">' + (mismatchDelta === null ? 'n/a' : mismatchDelta.toFixed(2)) + '</td><td>' + escapeHtml(previous.No_Hit_Rate) + '</td><td>' + escapeHtml(current.No_Hit_Rate) + '</td></tr>';
        }).join('');
        return '<article class="previous-report-panel"><h3>' + escapeHtml(report.name || report.path || 'Previous report') + '</h3>' +
          warnings +
          '<div class="primer-grid"><div>Rows<strong><br>' + totalRows + '</strong></div><div>Primers<strong><br>' + primers + '</strong></div><div>No-hit rate<strong><br>' + formatPercent(totalRows ? noHits / totalRows : 0) + '</strong></div><div>Avg mismatches<strong><br>' + avgMismatch + '</strong></div></div>' +
          '<p class="section-note">Top risky primers: ' + (risky.length ? risky.map(row => escapeHtml(row.Primer_Name) + ' (' + row.Risk + ')').join(', ') : 'none') + '</p>' +
          (comparisonRows ? '<div class="table-wrap"><table class="comparison-table"><thead><tr><th>Primer</th><th>Previous avg mismatches</th><th>Current avg mismatches</th><th>Delta</th><th>Previous no-hit rate</th><th>Current no-hit rate</th></tr></thead><tbody>' + comparisonRows + '</tbody></table></div>' : '<div class="empty">No overlapping primer names found for comparison with current filtered rows.</div>') +
          '</article>';
      }).join('');
    }
    function renderPrimerPanels(data) {
      const panels = groupByPrimer(data).map(([primer, group]) => {
        const stats = calculatePrimerStats(group);
        const risk = calculateRisk(stats);
        const hitPercent = group.length ? Math.round((stats.hits / group.length) * 100) : 0;
        return '<article class="primer-panel risk-' + risk.toLowerCase() + '-panel"><h3>' + escapeHtml(primer) + ' ' + riskBadge(risk, riskExplanation(stats)) + '</h3>' +
          '<div class="primer-grid">' +
          '<div>Rows<strong><br>' + group.length + '</strong></div>' +
          '<div>Hits<strong><br>' + stats.hits + '</strong></div>' +
          '<div>No hits<strong><br>' + stats.noHits + '</strong></div>' +
          '<div>Avg mismatches<strong><br>' + (stats.hits ? stats.avgMismatches.toFixed(2) : 'n/a') + '</strong></div>' +
          '<div>Max mismatches<strong><br>' + stats.maxMismatches + '</strong></div>' +
          '<div>Avg identity<strong><br>' + (stats.avgPercentIdentity === null ? 'n/a' : stats.avgPercentIdentity.toFixed(2)) + '</strong></div>' +
          '<div>Terminal mismatch share<strong><br>' + formatPercent(stats.terminalMismatchRate) + '</strong></div>' +
          '</div><div class="bar" aria-label="Hit percentage"><span style="width:' + hitPercent + '%"></span></div>' +
          renderSequenceMap(group) +
          renderTimeline(group) +
          '<div class="plot-card"><h3>Mismatch count distribution for this primer</h3>' + distributionChart(stats.distribution) + '</div>' +
          renderDetailTable(primer, group) +
          '</article>';
      }).join('');
      document.getElementById('primer-panels').innerHTML = panels || '<div class="empty">No rows match the current filters.</div>';
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
      renderPreviousReports(data);
      renderPrimerPanels(data);
    }
    fillSelect(filters.virus, uniqueValues('Virus_Type'), 'organisms');
    fillSelect(filters.primer, uniqueValues('Primer_Name'), 'primers');
    fillSelect(filters.fasta, uniqueValues('Fasta_File'), 'FASTA files');
    fillSelect(filters.status, uniqueValues('Hit_Status'), 'statuses');
    Object.values(filters).forEach(control => control.addEventListener('input', render));
    document.getElementById('alignment-modal-close').addEventListener('click', closeAlignmentModal);
    document.getElementById('alignment-modal').addEventListener('click', event => {
      if (event.target.id === 'alignment-modal') closeAlignmentModal();
    });
    document.addEventListener('keydown', event => {
      if (event.key === 'Escape') closeAlignmentModal();
    });
    render();
  </script>
</body>
</html>
"""
    return (
        template
        .replace("__REPORT_TITLE__", escaped_title)
        .replace("__REPORT_DATA__", data_json)
        .replace("__PREVIOUS_REPORT_DATA__", previous_json)
        .replace("__FIELDS_JSON__", fields_json)
    )


def write_html_report(results: list[dict], output_file: str, previous_reports: list[dict] | None = None):
    if not results:
        print("No results to write to HTML report.", file=sys.stderr)
        return
    try:
        with open(output_file, "w", encoding="utf-8") as htmlfile:
            htmlfile.write(build_html_report(results, previous_reports=previous_reports))
        print(f"HTML report successfully written to {output_file}")
    except Exception as e:
        print(f"Error writing HTML report: {e}", file=sys.stderr)
