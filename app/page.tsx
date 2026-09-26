"use client";
import { useLanguage, LanguageSwitch } from "./components/language";
import Documentation from "./components/documentation";
import { FastaFormatHelp } from "./components/input-format-help";

import {
  useEffect,
  useMemo,
  useRef,
  useState,
  useSyncExternalStore,
} from "react";
import Link from "next/link";
import {
  ArrowDownToLine,
  ArrowRight,
  Check,
  ChevronRight,
  Dna,
  FileText,
  FlaskConical,
  Info,
  LoaderCircle,
  Plus,
  ShieldCheck,
  Upload,
  X,
} from "lucide-react";
import type { Analysis, Catalog } from "./lib/results";
import { validateFiles } from "./lib/results";
import Results from "./components/results";
import DatabasePicker from "./components/database-picker";

const fallbackLimits = {
  upload_bytes: 3_000_000,
  files: 10,
  records: 200,
  comparisons: 2000,
  blast_calls: 300,
  base_comparisons: 50_000_000,
};
function subscribeView(callback: () => void) {
  window.addEventListener("hashchange", callback);
  return () => window.removeEventListener("hashchange", callback);
}
export default function Home() {
  const { t, locale } = useLanguage();
  const documentation = useSyncExternalStore(
    subscribeView,
    () => window.location.hash.startsWith("#documentation"),
    () => false,
  );
  useEffect(() => {
    if (documentation)
      document.getElementById(window.location.hash.slice(1))?.scrollIntoView();
  }, [documentation]);
  const customActive = useRef(false);
  const [preflightAttempt, setPreflightAttempt] = useState(0);
  const [bundled, setBundled] = useState<Catalog | null>(null);
  const [customCatalog, setCustomCatalog] = useState<
    Catalog | null | undefined
  >(undefined);
  const catalog = customCatalog === undefined ? bundled : customCatalog;
  const [database, setDatabase] = useState<File | null>(null);
  const [databaseReset, setDatabaseReset] = useState(0);
  const [preflight, setPreflight] = useState<{
    request: FormData;
    error?: string;
    workload?: {
      records: number;
      primers: number;
      comparisons: number;
      blast_calls: number;
      base_comparisons: number;
    };
    warnings?: string[];
  } | null>(null);
  const [catalogError, setCatalogError] = useState("");
  const [files, setFiles] = useState<File[]>([]);
  const [metadata, setMetadata] = useState<File | null>(null);
  const [virus, setVirus] = useState("SARS-CoV-2");
  const [subtype, setSubtype] = useState("H3");
  const [assayType, setAssayType] = useState("pcr");
  const [assay, setAssay] = useState("");
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [dragging, setDragging] = useState(false);
  const [result, setResult] = useState<Analysis | null>(null);
  const input = useRef<HTMLInputElement>(null);
  const resultHeading = useRef<HTMLDivElement>(null);

  const [connectionAttempt, setConnectionAttempt] = useState(0);
  useEffect(() => {
    const controller = new AbortController();
    fetch("/api/catalog", { signal: controller.signal })
      .then(async (response) => {
        if (!response.ok)
          throw new Error(
            "The primer database is unavailable. Check that the API is running and try again.",
          );
        return (await response.json()) as Catalog;
      })
      .then((data) => {
        setBundled(data);
        setCatalogError("");
        if (!customActive.current)
          setVirus((current) =>
            data.viruses.some((v) => v.id === current)
              ? current
              : data.viruses[0]?.id || "",
          );
      })
      .catch((error) => {
        if (!controller.signal.aborted)
          setCatalogError(
            error instanceof Error
              ? error.message
              : "Could not connect to the API.",
          );
      });
    return () => controller.abort();
  }, [connectionAttempt]);
  const selectedVirus = catalog?.viruses.find((v) => v.id === virus);
  const organism =
    virus === "influenza"
      ? subtype === "B"
        ? "Influenza-B"
        : "Influenza-A"
      : virus;
  const assays =
    catalog?.assays.filter(
      (a) =>
        a.organism === organism &&
        (assayType === "all" || a.type === assayType),
    ) || [];
  const limits = catalog?.limits || fallbackLimits;
  const bytes = files.reduce(
    (n, f) => n + f.size,
    (metadata?.size || 0) + (database?.size || 0),
  );
  const uploadProblem = validateFiles(files, metadata, limits, database);
  const uploadData = useMemo(() => {
    const data = new FormData();
    files.forEach((file) => data.append("files", file));
    if (metadata) data.append("metadata", metadata);
    if (database) data.append("database", database);
    data.append("virus", virus);
    data.append("assay_type", assayType);
    if (virus === "influenza") data.append("flu_type", subtype);
    if (assay) data.append("assay_id", assay);
    return data;
  }, [files, metadata, database, virus, assayType, subtype, assay]);
  const checked = preflight?.request === uploadData ? preflight : null;
  const readyToAnalyze =
    !!catalog && !!files.length && !uploadProblem && !!checked?.workload;
  useEffect(() => {
    if (!catalog || !files.length || uploadProblem) return;
    const controller = new AbortController();
    const timer = setTimeout(async () => {
      try {
        const response = await fetch("/api/preflight", {
          method: "POST",
          body: uploadData,
          signal: AbortSignal.any([
            controller.signal,
            AbortSignal.timeout(30_000),
          ]),
        });
        const payload = await response.json().catch(() => null);
        if (!response.ok || !payload?.workload)
          throw new Error(
            payload?.error?.message ||
              "Could not check this batch. Try again or use fewer files.",
          );
        if (!controller.signal.aborted)
          setPreflight({ request: uploadData, ...payload });
      } catch (e) {
        if (!controller.signal.aborted)
          setPreflight({
            request: uploadData,
            error:
              e instanceof Error ? e.message : "Could not check this batch.",
          });
      }
    }, 400);
    return () => {
      clearTimeout(timer);
      controller.abort();
    };
  }, [catalog, files.length, uploadData, uploadProblem, preflightAttempt]);

  function changeDatabase(
    file: File | null,
    data: Catalog | null,
    useReference = false,
  ) {
    customActive.current = !useReference;
    setDatabase(file);
    setCustomCatalog(useReference ? undefined : data);
    setPreflight(null);
    setAssay("");
    setError("");
    if (data) {
      const selected =
        data.viruses.find((v) => v.id === virus) || data.viruses[0];
      setVirus(selected?.id || "");
      setSubtype(selected?.subtypes[0] || "A");
      setAssayType(data.assays[0]?.type || "pcr");
    }
  }

  function addFiles(incoming: File[]) {
    const next = [...files, ...incoming];
    const problem = validateFiles(next, metadata, limits, database);
    if (problem) {
      setError(problem);
      return;
    }
    setFiles(next);
    setError("");
  }
  async function useExample() {
    try {
      const response = await fetch("/example.fasta");
      if (!response.ok) throw new Error();
      setFiles([
        new File([await response.text()], "example.fasta", {
          type: "text/plain",
        }),
      ]);
      customActive.current = false;
      setCustomCatalog(undefined);
      setDatabase(null);
      setDatabaseReset((n) => n + 1);
      setVirus("SARS-CoV-2");
      setAssayType("pcr");
      setAssay("");
      setMetadata(null);
      setError("");
    } catch {
      setError("The example could not be loaded. Try uploading a FASTA file.");
    }
  }
  async function analyze(event: React.FormEvent) {
    event.preventDefault();
    const problem = validateFiles(files, metadata, limits, database);
    if (problem || !files.length) {
      setError(problem || "Choose at least one FASTA file.");
      return;
    }
    if (!readyToAnalyze) return;
    setBusy(true);
    setError("");
    try {
      const response = await fetch("/api/analyze", {
        method: "POST",
        body: uploadData,
        signal: AbortSignal.timeout(285_000),
      });
      const payload = await response.json().catch(() => null);
      if (!response.ok)
        throw new Error(
          payload?.error?.message ||
            (response.status === 413
              ? "The upload or result is too large. Try fewer files."
              : response.status === 504
                ? "Analysis timed out. Try fewer files or a single assay."
                : "The analysis service could not complete the request. Please try again."),
        );
      if (!payload?.rows)
        throw new Error(
          "The service returned an incomplete result. Please try again.",
        );
      setResult(payload);
      requestAnimationFrame(() => {
        if (!resultHeading.current?.closest("[hidden]")) {
          resultHeading.current?.scrollIntoView({ behavior: "smooth" });
          resultHeading.current?.focus();
        }
      });
    } catch (e) {
      setError(
        e instanceof DOMException && e.name === "TimeoutError"
          ? "Analysis timed out. Try fewer files or one assay."
          : e instanceof Error
            ? e.message
            : "Could not connect to the analysis service.",
      );
    } finally {
      setBusy(false);
    }
  }
  return (
    <>
      <a className="skip-link" href="#main">
        {t("Skip to content")}
      </a>
      <header className="topbar">
        <div className="topbar-inner">
          <Link
            className="brand"
            href="/"
            aria-label={t("Primer Checker home")}
          >
            <span className="brand-mark">
              <Dna size={22} />
            </span>
            <span>
              Primer Checker<span className="brand-divider">/</span>
              <span className="brand-subtitle">
                {t("Sequence compatibility")}
              </span>
            </span>
          </Link>
          <div className="site-controls">
            <nav aria-label={t("Main navigation")}>
              <a
                href="#analysis"
                aria-current={!documentation ? "page" : undefined}
              >
                {t("Analysis")}
              </a>
              <a
                href="#documentation"
                aria-current={documentation ? "page" : undefined}
              >
                {t("Documentation")}
              </a>
              <a
                className="github-link"
                href="https://github.com/RasmusKoRiis/primer-checker"
                target="_blank"
                rel="noreferrer"
              >
                GitHub <ArrowRight size={13} />
              </a>
            </nav>
            <LanguageSwitch />
          </div>
        </div>
      </header>
      <main id="main" className="workspace">
        <div id="analysis" hidden={documentation}>
          <div className="page-heading">
            <div>
              <div className="eyebrow">
                <span className="status-dot" />
                {t("CONSENSUS SEQUENCE ANALYSIS")}
              </div>
              <h1>{t(result ? "Analysis workspace" : "New analysis")}</h1>
              <p>
                {t(
                  "Evaluate diagnostic and sequencing primer compatibility against viral consensus sequences.",
                )}
              </p>
            </div>
            <div className="db-stamp">
              <span>{t("PRIMER DATABASE")}</span>
              <strong>
                {catalog?.database.version ||
                  t(customCatalog === null ? "Not selected" : "Connecting…")}
              </strong>
              <small>
                {catalog
                  ? database?.name || t("Reference library")
                  : t("Choose a valid database")}
              </small>
            </div>
          </div>
          <div className="workflow" aria-label={t("Analysis workflow")}>
            <span className={files.length ? "complete" : "current"}>
              <i>{files.length ? <Check size={12} /> : "1"}</i>
              {t("Upload sequences")}
            </span>
            <ChevronRight size={14} />
            <span className={files.length ? "current" : ""}>
              <i>2</i>
              {t("Configure analysis")}
            </span>
            <ChevronRight size={14} />
            <span className={result ? "complete" : ""}>
              <i>{result ? <Check size={12} /> : "3"}</i>
              {t("Explore results")}
            </span>
          </div>
          {catalogError && (
            <div className="error-banner" role="alert">
              {t(catalogError)}{" "}
              <button onClick={() => setConnectionAttempt((n) => n + 1)}>
                {t("Retry connection")}
              </button>
            </div>
          )}
          <form onSubmit={analyze} aria-busy={busy}>
            <fieldset disabled={busy} className="form-reset">
              <DatabasePicker
                key={databaseReset}
                bundled={bundled}
                disabled={busy}
                onChange={changeDatabase}
              />
              <div className="input-grid">
                <section
                  className="card upload-card"
                  aria-labelledby="sequences-title"
                >
                  <div className="card-heading">
                    <div className="section-icon">
                      <FileText size={18} />
                    </div>
                    <div>
                      <h2 id="sequences-title">{t("Sequence files")}</h2>
                      <p>
                        {t("One or more consensus sequences in FASTA format")}
                      </p>
                    </div>
                    <span className="label-tag">{t("REQUIRED")}</span>
                  </div>
                  <input
                    ref={input}
                    id="fasta-upload"
                    className="sr-only"
                    type="file"
                    accept=".fasta,.fa,.fas,.fna"
                    multiple
                    aria-label={t("Upload FASTA files")}
                    onChange={(e) => {
                      addFiles(Array.from(e.target.files || []));
                      e.target.value = "";
                    }}
                  />
                  <div
                    className={`dropzone ${dragging ? "dragging" : ""}`}
                    onDragOver={(e) => {
                      e.preventDefault();
                      if (!busy) setDragging(true);
                    }}
                    onDragLeave={() => setDragging(false)}
                    onDrop={(e) => {
                      e.preventDefault();
                      setDragging(false);
                      if (!busy) addFiles(Array.from(e.dataTransfer.files));
                    }}
                  >
                    <div className="upload-icon">
                      <Upload size={25} strokeWidth={1.5} />
                    </div>
                    <h3>{t("Drop your FASTA files here")}</h3>
                    <p>{t("or browse files from your computer")}</p>
                    <button
                      type="button"
                      className="button secondary"
                      onClick={() => input.current?.click()}
                    >
                      <Plus size={15} />
                      {t("Choose files")}
                    </button>
                    <small>
                      {t(".fasta, .fa, .fna, .fas")}
                      <span>·</span>
                      {t("Up to {count} files, 3 MB combined", {
                        count: limits.files,
                      })}
                    </small>
                  </div>
                  {files.length > 0 && (
                    <ul className="file-list">
                      {files.map((file, i) => (
                        <li key={`${file.name}-${i}`}>
                          <FileText size={16} />
                          <span>{file.name}</span>
                          <small>{(file.size / 1000).toFixed(1)} kB</small>
                          <button
                            type="button"
                            className="icon-button"
                            aria-label={t("Remove {name}", { name: file.name })}
                            onClick={() =>
                              setFiles(files.filter((_, n) => n !== i))
                            }
                          >
                            <X size={15} />
                          </button>
                        </li>
                      ))}
                    </ul>
                  )}
                  <div className="example-row">
                    <span>{t("Just exploring?")}</span>
                    <button
                      type="button"
                      className="text-button"
                      onClick={useExample}
                    >
                      {t("Use a synthetic example")}
                      <ArrowRight size={13} />
                    </button>
                    <a
                      href="/example.fasta"
                      download
                      aria-label={t("Download synthetic FASTA example")}
                    >
                      <ArrowDownToLine size={15} />
                    </a>
                  </div>
                  <FastaFormatHelp />
                  <div className="metadata-section">
                    <div>
                      <h3>
                        {t("Sample metadata")}{" "}
                        <span className="optional">{t("Optional")}</span>
                      </h3>
                      <p>
                        {t(
                          "Attach a CSV with SampleID, Sample_Date, and Ct values.",
                        )}
                      </p>
                    </div>
                    <label
                      className="button secondary metadata-label"
                      htmlFor="metadata-upload"
                    >
                      <Plus size={14} />
                      {t("Add CSV")}
                      <input
                        id="metadata-upload"
                        type="file"
                        accept=".csv"
                        className="sr-only"
                        onChange={(e) => {
                          const next = e.target.files?.[0] || null;
                          const problem = validateFiles(
                            files,
                            next,
                            limits,
                            database,
                          );
                          if (problem) setError(problem);
                          else {
                            setMetadata(next);
                            setError("");
                          }
                          e.target.value = "";
                        }}
                      />
                    </label>
                  </div>
                  {metadata && (
                    <div className="metadata-file">
                      <FileText size={15} />
                      <span>{metadata.name}</span>
                      <button
                        type="button"
                        className="icon-button"
                        aria-label={t("Remove metadata")}
                        onClick={() => setMetadata(null)}
                      >
                        <X size={14} />
                      </button>
                    </div>
                  )}
                </section>
                <section
                  className="card settings-card"
                  aria-labelledby="settings-title"
                >
                  <div className="card-heading">
                    <div className="section-icon">
                      <FlaskConical size={18} />
                    </div>
                    <div>
                      <h2 id="settings-title">{t("Analysis settings")}</h2>
                      <p>{t("Select the primers to evaluate")}</p>
                    </div>
                  </div>
                  <label className="field">
                    {t("Virus")}
                    <select
                      value={virus}
                      onChange={(e) => {
                        setVirus(e.target.value);
                        setAssay("");
                        const options = catalog?.viruses.find(
                          (v) => v.id === e.target.value,
                        )?.subtypes;
                        if (options?.length)
                          setSubtype(
                            options.includes("H3") ? "H3" : options[0],
                          );
                      }}
                      disabled={!catalog || busy}
                    >
                      {catalog?.viruses.map((v) => (
                        <option key={v.id} value={v.id}>
                          {v.name}
                        </option>
                      ))}
                    </select>
                  </label>
                  {selectedVirus?.subtypes.length ? (
                    <label className="field">
                      {t("Influenza subtype")}
                      <select
                        value={subtype}
                        onChange={(e) => {
                          setSubtype(e.target.value);
                          setAssay("");
                        }}
                      >
                        {selectedVirus.subtypes.map((s) => (
                          <option key={s}>{s}</option>
                        ))}
                      </select>
                    </label>
                  ) : null}
                  <div className="field">
                    <span id="assay-type-label">{t("Assay type")}</span>
                    <div
                      className="segmented"
                      role="group"
                      aria-labelledby="assay-type-label"
                    >
                      {["pcr", "ngs", "all"].map((type) => (
                        <button
                          type="button"
                          key={type}
                          aria-pressed={assayType === type}
                          className={assayType === type ? "active" : ""}
                          onClick={() => {
                            setAssayType(type);
                            setAssay("");
                          }}
                        >
                          {type === "all"
                            ? t("All assays")
                            : type.toUpperCase()}
                        </button>
                      ))}
                    </div>
                  </div>
                  <label className="field">
                    {t("Primer scheme / panel")}
                    <select
                      value={assay}
                      onChange={(e) => setAssay(e.target.value)}
                      disabled={!catalog || busy}
                    >
                      <option value="">
                        {t(
                          assayType === "ngs"
                            ? "All matching panels"
                            : "All matching schemes",
                        )}
                      </option>
                      {assays.map((a) => (
                        <option key={`${a.type}-${a.id}`} value={a.id}>
                          {a.name}
                          {virus === "influenza" &&
                          ["H1", "H3"].includes(subtype)
                            ? ""
                            : t(" ({count} primers)", { count: a.primers })}
                        </option>
                      ))}
                    </select>
                  </label>
                  <div className="settings-note">
                    <Info size={15} />
                    <p>
                      {t(
                        virus === "influenza"
                          ? "Use segment labels in FASTA headers, such as 01-HA|sample. H1/H3 include the matching subtype and untagged influenza A primers."
                          : "Sequences are compared with the selected primer database using BLASTn and IUPAC-aware mismatch matching.",
                      )}
                    </p>
                  </div>
                  <div className="settings-bottom">
                    <span className="tiny-dot" />
                    {t("Shared with the command-line analysis engine")}
                  </div>
                </section>
              </div>
              <section
                className="batch-check card"
                aria-label={t("Batch limits")}
              >
                <div>
                  <strong>{t("Check before analysis")}</strong>
                  <p>
                    {t(
                      "Maximum {mb} MB total, including the database · {files} files · {records} records · {comparisons} comparisons · {searches} BLAST searches · {bases} million bases × primers.",
                      {
                        mb: (limits.upload_bytes / 1_000_000).toFixed(0),
                        files: limits.files,
                        records: limits.records,
                        comparisons: limits.comparisons.toLocaleString(locale),
                        searches: limits.blast_calls || 300,
                        bases: (
                          (limits.base_comparisons || 50_000_000) / 1_000_000
                        ).toFixed(0),
                      },
                    )}
                  </p>
                </div>
                <div aria-live="polite">
                  {uploadProblem ? (
                    <p className="batch-error">{t(uploadProblem)}</p>
                  ) : !catalog ? (
                    <p>
                      {t(
                        "Create or select a valid primer database to check your batch.",
                      )}
                    </p>
                  ) : !files.length ? (
                    <p>
                      {t(
                        "Add sequences to check the workload. Large batches can use the CLI.",
                      )}
                    </p>
                  ) : checked?.error ? (
                    <p className="batch-error">
                      {t(checked.error)}{" "}
                      <button
                        type="button"
                        className="text-button"
                        onClick={() => {
                          setPreflight(null);
                          setPreflightAttempt((n) => n + 1);
                        }}
                      >
                        {t("Retry check")}
                      </button>
                    </p>
                  ) : checked?.workload ? (
                    <p className="batch-ready">
                      <Check size={16} />
                      {t(
                        "Ready: {records} records · {primers} primers · {comparisons} comparisons · {searches} BLAST searches",
                        {
                          records: checked.workload.records,
                          primers: checked.workload.primers,
                          comparisons:
                            checked.workload.comparisons.toLocaleString(locale),
                          searches: checked.workload.blast_calls,
                        },
                      )}
                    </p>
                  ) : (
                    <p>
                      <LoaderCircle size={15} className="spin" />
                      {t("Checking files and selected primers…")}
                    </p>
                  )}
                  {!!catalog &&
                    !uploadProblem &&
                    checked?.warnings?.map((warning, i) => (
                      <p key={i}>{t(warning)}</p>
                    ))}
                </div>
              </section>
              <div className="run-bar">
                <div className="privacy-note">
                  <ShieldCheck size={18} />
                  <p>
                    {t(
                      "Do not upload confidential, identifiable, or otherwise restricted data to this public deployment.",
                    )}
                    <span>
                      {t(
                        "Uploads are processed temporarily and removed after analysis.",
                      )}
                    </span>
                  </p>
                </div>
                <div className="run-action">
                  <span>
                    {files.length} {t(files.length === 1 ? "file" : "files")}{" "}
                    {t("selected ·")} {(bytes / 1000).toFixed(0)} kB
                  </span>
                  <button
                    className="button primary"
                    type="submit"
                    disabled={!readyToAnalyze || busy}
                  >
                    {busy ? (
                      <>
                        <LoaderCircle className="spin" size={17} />
                        {t("Analyzing…")}
                      </>
                    ) : (
                      <>
                        {t("Analyze sequences")}
                        <ArrowRight size={17} />
                      </>
                    )}
                  </button>
                </div>
              </div>
            </fieldset>
            {error && (
              <div className="error-banner" role="alert">
                {t(error)}
              </div>
            )}
            {busy && (
              <div className="analysis-status" role="status">
                <LoaderCircle className="spin" size={19} />
                <div>
                  <strong>{t("Analysis request in progress")}</strong>
                  <p>
                    {t(
                      "The server validates your files, runs BLASTn, and builds the reports. Larger panels may take a few minutes.",
                    )}
                  </p>
                </div>
              </div>
            )}
          </form>
          {result ? (
            <div ref={resultHeading} tabIndex={-1} className="results-anchor">
              <Results key={result.manifest.analysis_utc} analysis={result} />
            </div>
          ) : (
            <section
              className="before-results"
              aria-label={t("About the results")}
            >
              <div>
                <span className="mini-number">01</span>
                <h3>{t("Compare every primer")}</h3>
                <p>
                  {t(
                    "Review matches and mismatch patterns across your selected sequences.",
                  )}
                </p>
              </div>
              <div>
                <span className="mini-number">02</span>
                <h3>{t("Inspect the alignment")}</h3>
                <p>
                  {t(
                    "Explore individual bases, substitutions, and primer mismatch positions.",
                  )}
                </p>
              </div>
              <div>
                <span className="mini-number">03</span>
                <h3>{t("Keep a reproducible record")}</h3>
                <p>
                  {t(
                    "Download CSV results, a standalone HTML report, and analysis provenance.",
                  )}
                </p>
              </div>
            </section>
          )}
        </div>
        <div id="documentation" hidden={!documentation}>
          <Documentation />
        </div>
        <footer>
          <span>
            <Dna size={15} /> Primer Checker{" "}
            <span className="footer-version">
              {t("v")}
              {catalog?.application_version || "0.1.0"}
            </span>
          </span>
          <span>{t("Consensus sequences · PCR & NGS compatibility")}</span>
          <a href="https://github.com/RasmusKoRiis/primer-checker/issues">
            {t("Report an issue")}
            <ArrowRight size={12} />
          </a>
        </footer>
      </main>
    </>
  );
}
