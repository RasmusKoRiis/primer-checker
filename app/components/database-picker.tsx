"use client";
import { useLanguage } from "./language";
import { DatabaseFormatHelp } from "./input-format-help";

import { useEffect, useRef, useState } from "react";
import {
  Database,
  Download,
  PixelLoader,
  Plus,
  Trash2,
  Upload,
} from "./pixel-icons";
import type { Catalog } from "../lib/results";
import { download } from "../lib/results";
import { makeDatabase, type PrimerDraft } from "../lib/database";
import dummyDatabase from "../../primer_db/dummy_primers.json";

const emptyPrimer = (key: number): PrimerDraft => ({
  key,
  name: "",
  sequence: "",
  role: "primer",
  segment: "",
  pool: "",
  subtype: "",
});

export default function DatabasePicker({
  bundled,
  onChange,
  disabled,
}: {
  bundled: Catalog | null;
  onChange: (
    file: File | null,
    catalog: Catalog | null,
    useReference?: boolean,
  ) => void;
  disabled: boolean;
}) {
  const { t } = useLanguage();
  const [mode, setMode] = useState("bundled");
  const [name, setName] = useState("");
  const [organism, setOrganism] = useState("SARS-CoV-2");
  const [version, setVersion] = useState("1.0");
  const [assayType, setAssayType] = useState("pcr");
  const [primers, setPrimers] = useState<PrimerDraft[]>([emptyPrimer(0)]);
  const [ready, setReady] = useState<{ file: File; catalog: Catalog } | null>(
    null,
  );
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const request = useRef(0);
  const controller = useRef<AbortController | null>(null);
  useEffect(
    () => () => {
      request.current++;
      controller.current?.abort();
    },
    [],
  );
  const nextKey = useRef(1);
  const input = useRef<HTMLInputElement>(null);
  const limits = bundled?.limits;

  function invalidate() {
    request.current++;
    controller.current?.abort();
    setReady(null);
    setError("");
    setLoading(false);
    onChange(null, null);
  }
  function changeMode(next: string) {
    invalidate();
    setMode(next);
    if (next === "bundled") onChange(null, bundled, true);
  }
  function updatePrimer(key: number, field: keyof PrimerDraft, value: string) {
    invalidate();
    setPrimers((current) =>
      current.map((p) => (p.key === key ? { ...p, [field]: value } : p)),
    );
  }
  async function validate(file: File) {
    invalidate();
    if (
      !/\.json$/i.test(file.name) ||
      !file.size ||
      file.size > (limits?.database_bytes || 250_000)
    ) {
      setError("Choose a nonempty JSON database of at most 250 kB.");
      return;
    }
    const id = request.current;
    const currentController = new AbortController();
    controller.current = currentController;
    setLoading(true);
    try {
      const body = new FormData();
      body.append("database", file);
      const response = await fetch("/api/database", {
        method: "POST",
        body,
        signal: AbortSignal.any([
          currentController.signal,
          AbortSignal.timeout(30_000),
        ]),
      });
      const data = await response.json().catch(() => null);
      if (!response.ok || !data?.assays)
        throw new Error(
          data?.error?.message ||
            "Could not validate the database. Check the API connection and try again.",
        );
      if (id !== request.current) return;
      setReady({ file, catalog: data });
      onChange(file, data);
    } catch (e) {
      if (id === request.current)
        setError(
          e instanceof Error ? e.message : "Database validation failed.",
        );
    } finally {
      if (id === request.current) setLoading(false);
    }
  }
  function create() {
    try {
      const text = makeDatabase({
        name,
        organism,
        version,
        assayType,
        primers,
      });
      void validate(
        new File([text], "custom-primers.json", { type: "application/json" }),
      );
    } catch (e) {
      setError(e instanceof Error ? e.message : "Check the primer fields.");
    }
  }
  return (
    <section className="card database-card" aria-labelledby="database-title">
      <div className="card-heading">
        <div className="section-icon">
          <Database size={20} />
        </div>
        <div>
          <h2 id="database-title">{t("Primer database")}</h2>
          <p>{t("Try dummy data or bring your own primers")}</p>
        </div>
      </div>
      <div
        className="database-modes"
        role="group"
        aria-label={t("Database source")}
      >
        {[
          [
            "bundled",
            bundled?.database.is_dummy === false
              ? "Installed database"
              : "Dummy database",
          ],
          ["upload", "Upload JSON"],
          ["build", "Build a database"],
        ].map(([value, label]) => (
          <button
            key={value}
            type="button"
            className={`button ${mode === value ? "selected" : "secondary"}`}
            aria-pressed={mode === value}
            disabled={disabled}
            onClick={() => changeMode(value)}
          >
            {t(label)}
          </button>
        ))}
      </div>
      <DatabaseFormatHelp />
      {mode === "bundled" && (
        <div>
          {bundled?.database.is_dummy && (
            <div className="warning-banner" role="note">
              <strong>{t("Dummy data — testing only")}</strong>
              <p>
                {t(
                  "These invented sequences demonstrate the software. They are not a validated primer set. Upload or build your own database to check your sequences.",
                )}
              </p>
              <button
                type="button"
                className="button secondary"
                onClick={() =>
                  download(
                    JSON.stringify(dummyDatabase, null, 2) + "\n",
                    "dummy-primers.json",
                    "application/json",
                  )
                }
              >
                <Download size={16} />
                {t("Download dummy database")}
              </button>
            </div>
          )}
          <p className="database-hint">
            {t("Database version")}:{" "}
            {bundled?.database.version || t("loading…")}
          </p>
        </div>
      )}
      {mode === "upload" && (
        <div className="database-upload">
          <p className="database-hint">
            {t(
              "Upload a self-contained JSON database: normalized schemes with primer sequences, or a legacy organism → primer → sequence dictionary. Up to 250 kB and {count} primers. BED/FASTA asset references require the CLI.",
              { count: limits?.database_primers || 500 },
            )}
          </p>
          <input
            ref={input}
            className="sr-only"
            type="file"
            accept=".json,application/json"
            aria-label={t("Upload primer database")}
            onChange={(e) => {
              const file = e.target.files?.[0];
              e.target.value = "";
              if (file) void validate(file);
            }}
          />
          <button
            type="button"
            className="button secondary"
            onClick={() => input.current?.click()}
          >
            <Upload size={16} />
            {t("Choose database")}
          </button>
        </div>
      )}
      {mode === "build" && (
        <div className="database-builder">
          <p className="database-hint">
            {t(
              "Enter primers in 5′ → 3′ orientation, including reverse primers; do not reverse-complement them. DNA IUPAC ambiguity codes are accepted. Create the database to use it in this analysis and download a reusable JSON file.",
            )}
          </p>
          <div className="database-fields">
            <label className="field">
              {t("Database name")}
              <input
                value={name}
                maxLength={200}
                placeholder={t("My primer scheme")}
                onChange={(e) => {
                  invalidate();
                  setName(e.target.value);
                }}
              />
            </label>
            <label className="field">
              {t("Organism")}
              <input
                list="database-organisms"
                value={organism}
                maxLength={200}
                onChange={(e) => {
                  invalidate();
                  setOrganism(e.target.value);
                }}
              />
            </label>
            <datalist id="database-organisms">
              {[
                "SARS-CoV-2",
                "Influenza-A",
                "Influenza-B",
                "Influenza-C",
                "Influenza-D",
                "RSV-A",
                "RSV-B",
              ].map((v) => (
                <option key={v} value={v} />
              ))}
            </datalist>
            <label className="field">
              {t("Database version")}
              <input
                value={version}
                maxLength={200}
                onChange={(e) => {
                  invalidate();
                  setVersion(e.target.value);
                }}
              />
            </label>
            <label className="field">
              {t("Database assay type")}
              <select
                value={assayType}
                onChange={(e) => {
                  invalidate();
                  setAssayType(e.target.value);
                }}
              >
                <option value="pcr">PCR</option>
                <option value="ngs">NGS</option>
              </select>
            </label>
          </div>
          <div className="primer-drafts">
            {primers.map((primer, i) => (
              <div className="primer-draft" key={primer.key}>
                <div className="primer-draft-title">
                  <strong>
                    {t("Primer")} {i + 1}
                  </strong>
                  <button
                    type="button"
                    className="icon-button"
                    aria-label={t("Remove primer {number}", { number: i + 1 })}
                    disabled={primers.length === 1}
                    onClick={() => {
                      invalidate();
                      setPrimers(primers.filter((p) => p.key !== primer.key));
                    }}
                  >
                    <Trash2 size={16} />
                  </button>
                </div>
                <div className="primer-sequence-fields">
                  <label className="field">
                    {t("Name")}
                    <input
                      aria-label={t("Primer {number} name", { number: i + 1 })}
                      value={primer.name}
                      maxLength={200}
                      placeholder="Target_F"
                      onChange={(e) =>
                        updatePrimer(primer.key, "name", e.target.value)
                      }
                    />
                  </label>
                  <label className="field">
                    {t("Sequence · 5′ → 3′")}
                    <textarea
                      aria-label={t("Primer {number} sequence", {
                        number: i + 1,
                      })}
                      value={primer.sequence}
                      rows={2}
                      maxLength={1000}
                      spellCheck={false}
                      placeholder="ACGT…"
                      onChange={(e) =>
                        updatePrimer(primer.key, "sequence", e.target.value)
                      }
                    />
                  </label>
                </div>
                <div className="primer-metadata-fields">
                  <label className="field">
                    {t("Role")}
                    <select
                      aria-label={t("Primer {number} role", { number: i + 1 })}
                      value={primer.role}
                      onChange={(e) =>
                        updatePrimer(primer.key, "role", e.target.value)
                      }
                    >
                      <option value="primer">{t("Primer")}</option>
                      <option value="forward">{t("Forward primer")}</option>
                      <option value="reverse">{t("Reverse primer")}</option>
                      <option value="probe">{t("Probe")}</option>
                    </select>
                  </label>
                  <label className="field">
                    {t("Segment")}
                    <input
                      aria-label={t("Primer {number} segment", {
                        number: i + 1,
                      })}
                      value={primer.segment}
                      list="database-segments"
                      maxLength={32}
                      placeholder={t("For example PB2, HA, NA")}
                      onChange={(e) =>
                        updatePrimer(primer.key, "segment", e.target.value)
                      }
                    />
                  </label>
                  <label className="field">
                    {t("Pool (optional)")}
                    <input
                      aria-label={t("Primer {number} pool", { number: i + 1 })}
                      value={primer.pool}
                      maxLength={200}
                      onChange={(e) =>
                        updatePrimer(primer.key, "pool", e.target.value)
                      }
                    />
                  </label>
                  {/^influenza-/i.test(organism.trim()) && (
                    <label className="field">
                      {t("Subtype tags (optional)")}
                      <input
                        aria-label={t("Primer {number} subtype", {
                          number: i + 1,
                        })}
                        value={primer.subtype}
                        maxLength={1000}
                        placeholder={t("For example H5N1, H7N9")}
                        onChange={(e) =>
                          updatePrimer(primer.key, "subtype", e.target.value)
                        }
                      />
                    </label>
                  )}
                </div>
              </div>
            ))}
          </div>
          <datalist id="database-segments">
            {["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS", "HEF", "P3"].map(
              (s) => (
                <option key={s} value={s} />
              ),
            )}
          </datalist>
          <p className="database-hint">
            {t(
              "Up to {primers} primers, {bases} bases each. For influenza, enter a segment matching your FASTA headers. Any segment label with 1–32 letters or digits is accepted. Separate subtype tags with commas; leave blank for primers shared by all subtypes of that influenza type.",
              {
                primers: limits?.database_primers || 500,
                bases: limits?.primer_length || 200,
              },
            )}
          </p>
          <div className="database-actions">
            <button
              type="button"
              className="button secondary"
              disabled={primers.length >= (limits?.database_primers || 500)}
              onClick={() => {
                invalidate();
                setPrimers([...primers, emptyPrimer(nextKey.current++)]);
              }}
            >
              <Plus size={16} />
              {t("Add primer")}
            </button>
            <button
              type="button"
              className="button primary"
              disabled={loading}
              onClick={create}
            >
              {t("Create database")}
            </button>
          </div>
        </div>
      )}
      {loading && (
        <p className="database-hint" role="status">
          <PixelLoader size={16} />
          {t("Validating database…")}
        </p>
      )}
      {error && (
        <div className="error-banner" role="alert">
          {t(error)}
        </div>
      )}
      {ready && (
        <div className="database-ready" role="status">
          <div>
            <strong>
              {t("{filename} is ready for analysis", {
                filename: ready.file.name,
              })}
            </strong>
            <p>
              {t("{count} primers · version {version}", {
                count: ready.catalog.assays.reduce((n, a) => n + a.primers, 0),
                version: ready.catalog.database.version,
              })}
            </p>
          </div>
          <button
            type="button"
            className="button secondary"
            onClick={async () =>
              download(
                await ready.file.text(),
                ready.file.name,
                "application/json",
              )
            }
          >
            <Download size={16} />
            {t("Download database JSON")}
          </button>
        </div>
      )}
      {ready?.catalog.warnings?.map((warning, i) => (
        <p className="database-hint" key={i}>
          {t(warning)}
        </p>
      ))}
    </section>
  );
}
