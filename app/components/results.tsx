"use client";
import { useLanguage } from "./language";
import { useEffect, useMemo, useRef, useState } from "react";
import {
  ArrowDownToLine,
  ArrowUpDown,
  Check,
  ChevronLeft,
  ChevronRight,
  SlidersHorizontal,
  X,
} from "lucide-react";
import {
  aggregatePrimers,
  alignmentColumns,
  download,
  emptyFilters,
  filterRows,
} from "../lib/results";
import type { Analysis, Filters, ResultRow } from "../lib/results";

export default function Results({ analysis }: { analysis: Analysis }) {
  const { t, locale } = useLanguage();
  const [filters, setFilters] = useState<Filters>(emptyFilters);
  const [tab, setTab] = useState("primers");
  const [sort, setSort] = useState({ column: "primer", direction: 1 });
  const [page, setPage] = useState(0);
  const [selected, setSelected] = useState<ResultRow | null>(null);
  const dialog = useRef<HTMLDialogElement>(null);
  const filtered = useMemo(
    () => filterRows(analysis.rows, filters),
    [analysis.rows, filters],
  );
  const primers = useMemo(() => aggregatePrimers(filtered), [filtered]);
  const displayRows = tab === "primers" ? primers : filtered;
  const ordered = [...displayRows].sort((a, b) => {
    const av = (a as Record<string, string | number>)[sort.column] ?? "";
    const bv = (b as Record<string, string | number>)[sort.column] ?? "";
    return (
      (typeof av === "number" && typeof bv === "number"
        ? av - bv
        : String(av).localeCompare(String(bv), undefined, { numeric: true })) *
      sort.direction
    );
  });
  const pages = Math.max(1, Math.ceil(ordered.length / 25));
  const visiblePage = Math.min(page, pages - 1);
  const visible = ordered.slice(visiblePage * 25, (visiblePage + 1) * 25);
  const s = analysis.summary;
  const m = analysis.manifest;
  const segments = [
    ...new Set(
      analysis.rows
        .map((r) => r.Primer_Segment || r.Subject_Segment)
        .filter(Boolean),
    ),
  ];
  const assays = [
    ...new Map(
      analysis.rows.map((r) => [r.Assay_ID, r.Assay_Name || r.Assay_ID]),
    ).entries(),
  ];
  const mismatchCounts = [
    ...new Set(
      analysis.rows
        .filter((r) => r.Hit_Status === "hit")
        .map((r) => Number(r.Mismatches)),
    ),
  ].sort((a, b) => a - b);
  function update(key: keyof Filters, value: string) {
    setFilters({ ...filters, [key]: value });
    setPage(0);
  }
  function changeTab(value: string) {
    setTab(value);
    setPage(0);
    setSort({
      column: value === "primers" ? "primer" : "Subject_Sequence_ID",
      direction: 1,
    });
  }
  useEffect(() => {
    if (selected) dialog.current?.showModal();
  }, [selected]);
  function heading(label: string, column: string) {
    return (
      <th
        key={column}
        aria-sort={
          sort.column === column
            ? sort.direction === 1
              ? "ascending"
              : "descending"
            : "none"
        }
      >
        <button
          onClick={() =>
            setSort({
              column,
              direction: sort.column === column ? -sort.direction : 1,
            })
          }
        >
          {t(label)}
          <ArrowUpDown size={11} />
        </button>
      </th>
    );
  }
  return (
    <section className="results-section" aria-labelledby="results-title">
      <div className="results-heading">
        <div>
          <div className="eyebrow">
            <Check size={12} />
            {t("ANALYSIS COMPLETE")}
          </div>
          <h2 id="results-title">{t("Compatibility results")}</h2>
          <p>
            {m.selection.virus}
            {m.selection.flu_type ? ` · ${m.selection.flu_type}` : ""}{" "}
            <span>·</span> {m.selection.assay_type.toUpperCase()} <span>·</span>{" "}
            {new Date(m.analysis_utc).toLocaleString(locale)}
          </p>
        </div>
        <div className="download-group">
          <button
            className="button secondary"
            onClick={() =>
              download(
                analysis.downloads.csv,
                "primer-results.csv",
                "text/csv;charset=utf-8",
              )
            }
          >
            <ArrowDownToLine size={14} /> CSV
          </button>
          <button
            className="button secondary"
            onClick={() =>
              download(
                analysis.downloads.html,
                "primer-report.html",
                "text/html;charset=utf-8",
              )
            }
          >
            <ArrowDownToLine size={14} />
            {t("HTML report")}
          </button>
        </div>
      </div>
      {analysis.warnings.map((w) => (
        <div className="warning-banner" key={w}>
          {t(w)}
        </div>
      ))}
      <div className="stats-grid">
        <Stat
          label={t("Files / sequence records")}
          value={`${s.files} / ${s.sequence_records}`}
          note={t("{count} analyzed sample identifiers", { count: s.samples })}
        />
        <Stat
          label={t("Primers evaluated")}
          value={s.primers}
          note={t("{count} comparisons", {
            count: s.comparisons.toLocaleString(locale),
          })}
        />
        <Stat
          label={t("Successful hits")}
          value={s.hits}
          note={t("{count} without a hit", { count: s.comparisons - s.hits })}
        />
        <Stat
          label={t("Hits with mismatches")}
          value={s.mismatch_comparisons}
          note={t("{count} sample identifiers affected", {
            count: s.samples_affected,
          })}
          accented
        />
      </div>
      <div className="results-card card">
        <div className="table-toolbar">
          <div className="tabs" role="tablist" aria-label={t("Result views")}>
            <button
              role="tab"
              id="primers-tab"
              aria-controls="result-table"
              aria-selected={tab === "primers"}
              onClick={() => changeTab("primers")}
            >
              {t("By primer")} <span>{primers.length}</span>
            </button>
            <button
              role="tab"
              id="samples-tab"
              aria-controls="result-table"
              aria-selected={tab === "samples"}
              onClick={() => changeTab("samples")}
            >
              {t("By sample")} <span>{filtered.length}</span>
            </button>
          </div>
          <span className="filter-label">
            <SlidersHorizontal size={14} />
            {t("Filter results")}
          </span>
        </div>
        <div className="filters">
          <label>
            {t("Primer")}
            <input
              value={filters.primer}
              placeholder={t("Search primer…")}
              onChange={(e) => update("primer", e.target.value)}
            />
          </label>
          <label>
            {t("Sample")}
            <input
              value={filters.sample}
              placeholder={t("Search sample…")}
              onChange={(e) => update("sample", e.target.value)}
            />
          </label>
          <label>
            {t("Segment")}
            <select
              value={filters.segment}
              onChange={(e) => update("segment", e.target.value)}
            >
              <option value="">{t("All segments")}</option>
              {segments.map((x) => (
                <option key={x}>{x}</option>
              ))}
            </select>
          </label>
          <label>
            {t("Assay")}
            <select
              value={filters.assay}
              onChange={(e) => update("assay", e.target.value)}
            >
              <option value="">{t("All assays")}</option>
              {assays.map(([id, name]) => (
                <option key={id} value={id}>
                  {name}
                </option>
              ))}
            </select>
          </label>
          <label>
            {t("Mismatches")}
            <select
              value={filters.mismatches}
              onChange={(e) => update("mismatches", e.target.value)}
            >
              <option value="">{t("Any count")}</option>
              <option value="any">{t("At least one")}</option>
              {mismatchCounts.map((x) => (
                <option key={x} value={x}>
                  {x}
                </option>
              ))}
            </select>
          </label>
          <label>
            {t("Hit status")}
            <select
              value={filters.status}
              onChange={(e) => update("status", e.target.value)}
            >
              <option value="">{t("All results")}</option>
              <option value="hit">{t("Hit")}</option>
              <option value="no_hit">{t("No hit")}</option>
            </select>
          </label>
        </div>
        <div
          role="tabpanel"
          id="result-table"
          aria-labelledby={`${tab}-tab`}
          className="table-scroll"
        >
          <table>
            <caption className="sr-only">
              {t(
                "Primer compatibility results. Select a sample result to inspect its alignment.",
              )}
            </caption>
            <thead>
              <tr>
                {tab === "primers"
                  ? [
                      heading("Primer", "primer"),
                      heading("Assay", "assay"),
                      heading("Segment / pool", "segment"),
                      heading("Tested", "tested"),
                      heading("Perfect", "perfect"),
                      heading("Mismatches", "affected"),
                      heading("No hit", "noHit"),
                      heading("Max. mismatches", "maximum"),
                      heading("Affected", "percent"),
                    ]
                  : [
                      heading("Sample / file", "Subject_Sequence_ID"),
                      heading("Primer", "Primer_Name"),
                      heading("Segment", "Subject_Segment"),
                      heading("Identity", "Percent_Identity"),
                      heading("Mismatches", "Mismatches"),
                      heading("Positions", "Mismatch_Positions"),
                      <th key="inspect">{t("Alignment")}</th>,
                    ]}
              </tr>
            </thead>
            <tbody>
              {tab === "primers"
                ? (visible as ReturnType<typeof aggregatePrimers>).map((p) => (
                    <tr key={p.key}>
                      <td>
                        <button
                          className="table-link"
                          onClick={() => {
                            setFilters({ ...filters, primer: p.primer });
                            changeTab("samples");
                          }}
                        >
                          {p.primer}
                        </button>
                      </td>
                      <td className="muted">{p.assay || "—"}</td>
                      <td>
                        {p.segment || "—"}
                        <span className="cell-secondary">{p.pool}</span>
                      </td>
                      <td>{p.tested}</td>
                      <td>
                        <span className="match-count">{p.perfect}</span>
                      </td>
                      <td>
                        <span
                          className={p.affected ? "mismatch-count" : "muted"}
                        >
                          {p.affected}
                        </span>
                      </td>
                      <td>{p.noHit}</td>
                      <td>{p.noHit === p.tested ? "—" : p.maximum}</td>
                      <td>{p.percent.toFixed(1)}%</td>
                    </tr>
                  ))
                : (visible as ResultRow[]).map((r, i) => (
                    <tr
                      key={`${r.Fasta_File}-${r.Subject_Sequence_ID}-${r.Assay_ID}-${r.Primer_Name}-${i}`}
                    >
                      <td className="mono">
                        {r.Subject_Sequence_ID}
                        <span className="cell-secondary">{r.Fasta_File}</span>
                      </td>
                      <td>
                        {r.Primer_Name}
                        <span className="cell-secondary">{r.Assay_Name}</span>
                      </td>
                      <td>{r.Subject_Segment || r.Primer_Segment || "—"}</td>
                      <td>
                        {r.Hit_Status === "hit" ? (
                          `${Number(r.Percent_Identity).toFixed(1)}%`
                        ) : (
                          <span className="no-hit">{t("No hit")}</span>
                        )}
                      </td>
                      <td>
                        <span
                          className={
                            Number(r.Mismatches) > 0
                              ? "mismatch-count"
                              : "match-count"
                          }
                        >
                          {r.Hit_Status === "hit" ? r.Mismatches : "—"}
                        </span>
                      </td>
                      <td className="mono">{r.Mismatch_Positions || "—"}</td>
                      <td>
                        <button
                          className="table-link"
                          onClick={() => setSelected(r)}
                        >
                          {t("Inspect")}
                          <ChevronRight size={13} />
                        </button>
                      </td>
                    </tr>
                  ))}
            </tbody>
          </table>
          {!visible.length && (
            <div className="empty-table">
              <h3>{t("No results match these filters")}</h3>
              <button
                className="text-button"
                onClick={() => {
                  setFilters(emptyFilters);
                  setPage(0);
                }}
              >
                {t("Clear filters")}
              </button>
            </div>
          )}
        </div>
        <div className="table-pagination">
          <span>
            {t("{start}–{end} of {count} {kind}", {
              start: ordered.length ? visiblePage * 25 + 1 : 0,
              end: Math.min((visiblePage + 1) * 25, ordered.length),
              count: ordered.length,
              kind: t(tab === "primers" ? "primers" : "comparisons"),
            })}
          </span>
          <div>
            <button
              className="icon-button"
              aria-label={t("Previous page")}
              disabled={!visiblePage}
              onClick={() => setPage(visiblePage - 1)}
            >
              <ChevronLeft size={17} />
            </button>
            <span>
              {t("Page {page} of {pages}", { page: visiblePage + 1, pages })}
            </span>
            <button
              className="icon-button"
              aria-label={t("Next page")}
              disabled={visiblePage >= pages - 1}
              onClick={() => setPage(visiblePage + 1)}
            >
              <ChevronRight size={17} />
            </button>
          </div>
        </div>
      </div>
      <p className="interpretation-note">
        {t(
          "Tables summarize the filtered primer/sequence comparisons. “Affected” means at least one mismatch; no-hit results are shown separately. These are descriptive findings, not predictions of assay performance. Each FASTA record is counted separately, using its filename and sequence ID.",
        )}
      </p>
      <details className="provenance">
        <summary>
          {t("Analysis provenance")}{" "}
          <span>
            {t("Database {version} · App {app}", {
              version: m.database.version,
              app: m.application_version,
            })}{" "}
            · {m.blast.version}
          </span>
        </summary>
        <dl>
          <dt>{t("Selection")}</dt>
          <dd>
            {m.selection.virus} / {m.selection.flu_type || "—"} /{" "}
            {m.selection.assay_type} /{" "}
            {m.selection.assay_id || t("all matching assays")}
          </dd>
          <dt>{t("Database source")}</dt>
          <dd>
            {m.database.filename ||
              t(m.database.is_dummy ? "Dummy database" : "Installed database")}
          </dd>
          <dt>{t("Database SHA-256")}</dt>
          <dd className="mono">{m.database.sha256}</dd>
          <dt>{t("Git commit")}</dt>
          <dd className="mono">{m.git_commit}</dd>
          <dt>{t("Input files")}</dt>
          <dd>
            {m.files
              .map((f) =>
                t("{filename} ({count} records)", {
                  filename: f.filename,
                  count: f.records,
                }),
              )
              .join(", ")}
          </dd>
        </dl>
        <button
          className="button secondary"
          onClick={() =>
            download(
              analysis.downloads.manifest,
              "analysis-provenance.json",
              "application/json",
            )
          }
        >
          <ArrowDownToLine size={14} />
          {t("Download provenance")}
        </button>
      </details>
      <dialog
        ref={dialog}
        className="alignment-dialog"
        aria-labelledby="alignment-title"
        onClose={() => setSelected(null)}
        onClick={(e) => {
          if (e.target === dialog.current) dialog.current?.close();
        }}
      >
        {selected && (
          <>
            <div className="dialog-heading">
              <div>
                <div className="eyebrow">{t("ALIGNMENT INSPECTOR")}</div>
                <h2 id="alignment-title">{selected.Primer_Name}</h2>
                <p>{selected.Subject_Sequence_ID}</p>
              </div>
              <button
                className="icon-button"
                aria-label={t("Close alignment")}
                onClick={() => dialog.current?.close()}
              >
                <X size={19} />
              </button>
            </div>
            <div className="alignment-meta">
              <span>{selected.Virus_Type}</span>
              <span>{selected.Assay_Name}</span>
              <span>
                {selected.Hit_Status === "hit"
                  ? t(
                      "Local BLAST hit: {strand} strand · {start}–{end} (1-based, inclusive)",
                      {
                        strand: t(
                          Number(selected.Subject_End) >=
                            Number(selected.Subject_Start)
                            ? "Forward"
                            : "Reverse",
                        ),
                        start: selected.Subject_Start,
                        end: selected.Subject_End,
                      },
                    )
                  : t("No BLAST hit")}
              </span>
            </div>
            <p className="sequence-label">
              {t("Original primer sequence (5′ → 3′)")}
            </p>
            <div className="original-sequence mono">
              {selected.Primer_Sequence}
            </div>
            {selected.Hit_Status === "hit" ? (
              <>
                <Alignment row={selected} />
                <p className="section-note">
                  {t(
                    "Both rows follow the primer’s 5′ → 3′ direction. For a reverse-strand hit, the subject is reverse-complemented. Coordinates above refer to the local BLAST hit before extension to the full primer.",
                  )}
                </p>
                <p className="alignment-legend">
                  <span />
                  {t(
                    "Highlighted columns show base differences or gaps. Primer positions start at 1 at the 5′ end; the 3′ end is on the right. Insertions do not advance primer numbering.",
                  )}
                </p>
                <dl className="alignment-details">
                  <dt>{t("Mismatch positions")}</dt>
                  <dd>{selected.Mismatch_Positions || t("None")}</dd>
                  <dt>{t("Base differences")}</dt>
                  <dd className="mono">
                    {selected.Mismatch_Details || t("None")}
                  </dd>
                  <dt>{t("Identity")}</dt>
                  <dd>{Number(selected.Percent_Identity).toFixed(2)}%</dd>
                </dl>
              </>
            ) : (
              <p className="warning-banner">
                {t(
                  "BLAST did not return a hit for this primer/sequence pair. There is no alignment to inspect; this is separate from a measured mismatch count.",
                )}
              </p>
            )}
          </>
        )}
      </dialog>
    </section>
  );
}
function Stat({
  label,
  value,
  note,
  accented = false,
}: {
  label: string;
  value: string | number;
  note: string;
  accented?: boolean;
}) {
  return (
    <div className={`stat ${accented ? "accented" : ""}`}>
      <span>{label}</span>
      <strong>{value}</strong>
      <small>{note}</small>
    </div>
  );
}
function Alignment({ row }: { row: ResultRow }) {
  const { t } = useLanguage();
  const columns = alignmentColumns(row);
  return (
    <div
      className="alignment-scroll"
      aria-label={t("Primer and subject alignment")}
    >
      <div className="alignment-grid">
        <div className="alignment-labels">
          <span>{t("Position")}</span>
          <span>{t("Primer")}</span>
          <span>{t("Subject")}</span>
        </div>
        {columns.map((c, i) => (
          <div
            className={`alignment-column ${c.mismatch ? "mismatch" : ""}`}
            key={i}
          >
            <span>{c.base === "-" ? "·" : c.position}</span>
            <b>{c.base}</b>
            <b>{c.subject}</b>
          </div>
        ))}
      </div>
    </div>
  );
}
