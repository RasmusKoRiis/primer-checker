"use client";

import { useState } from "react";
import { Download } from "./pixel-icons";
import examples from "../lib/input-examples.json";
import { download } from "../lib/results";
import { useLanguage } from "./language";

export function FastaFormatHelp({ expanded = false }: { expanded?: boolean }) {
  const { t } = useLanguage();
  const [kind, setKind] = useState<keyof typeof examples.fasta>("general");
  const content = (
    <div className="format-guide">
      <p>
        {t(
          "Start each record with > and a unique sequence ID. Put the nucleotide sequence on the following line or lines.",
        )}
      </p>
      <div
        className="format-options"
        role="group"
        aria-label={t("FASTA example type")}
      >
        {(["general", "influenza"] as const).map((value) => (
          <button
            key={value}
            type="button"
            className="button secondary"
            aria-pressed={kind === value}
            onClick={() => setKind(value)}
          >
            {t(value === "general" ? "General FASTA" : "Influenza FASTA")}
          </button>
        ))}
      </div>
      <pre tabIndex={0} aria-label={t("FASTA example")}>
        <code>{examples.fasta[kind]}</code>
      </pre>
      {kind === "influenza" && (
        <div className="format-callout">
          <strong>{t("Influenza needs segment labels")}</strong>
          <p>
            {t(
              "Use a pipe-separated tag such as 01-PB2, 06-NA, or 08-NS before any spaces. The parser reads the segment after one or two digits and a hyphen; the digits do not determine the segment. Any label with 1–32 letters or digits is supported when it matches the primer database. Primers are compared only with that segment.",
            )}
          </p>
          <p>
            <code>{">sample_001 HA"}</code> / <code>{">sample_001|HA"}</code> —{" "}
            {t(
              "these do not identify the segment. Use >01-HA|sample_001 instead.",
            )}
          </p>
        </div>
      )}
      <ul>
        <li>
          {t(
            "Only the first word after > is the ID: >sample_001 optional description becomes sample_001. Descriptions after a space are ignored.",
          )}
        </li>
        <li>
          {t(
            "IDs must be unique within each file and at most 200 characters. Keep the full header within 1,000 characters. Avoid reserved prefixes gi|, lcl|, ref|, and gb|.",
          )}
        </li>
        <li>
          {t(
            "Use nucleotide IUPAC letters (for example A, C, G, T, N, R, Y), with no spaces or gaps inside sequence lines. Every header needs a sequence. FASTQ is not supported.",
          )}
        </li>
      </ul>
      <button
        type="button"
        className="button secondary"
        onClick={() =>
          download(examples.fasta[kind], `${kind}-template.fasta`, "text/plain")
        }
      >
        <Download size={16} />
        {t("Download FASTA template")}
      </button>
      <p className="format-footnote">
        {t(
          "Synthetic formatting examples only. Replace the IDs and sequences with your own data.",
        )}
      </p>
    </div>
  );
  if (expanded) return content;
  return (
    <div className="input-format-help">
      <p>
        {t("Header example:")} <code>{">sample_001"}</code> · {t("Influenza:")}{" "}
        <code>{">01-HA|sample_001"}</code>
      </p>
      <details>
        <summary>{t("FASTA header rules and examples")}</summary>
        {content}
      </details>
    </div>
  );
}

export function DatabaseFormatHelp({
  expanded = false,
}: {
  expanded?: boolean;
}) {
  const { t } = useLanguage();
  const [kind, setKind] =
    useState<keyof typeof examples.databases>("normalized");
  const text = JSON.stringify(examples.databases[kind], null, 2) + "\n";
  const content = (
    <div className="format-guide">
      <p>
        {t(
          "Start from a template below, or use Build a database to create the JSON without editing it by hand. Download the file, replace the example primers, then choose Upload a database.",
        )}
      </p>
      <p>
        {t(
          "Enter every primer as the oligo sequence in 5′ → 3′ direction, including reverse primers. Do not reverse-complement it before entry. The forward/reverse role is metadata; the analysis searches both strands of the uploaded sequence.",
        )}
      </p>
      <div
        className="format-options"
        role="group"
        aria-label={t("Database example type")}
      >
        {(
          [
            ["normalized", "PCR / NGS template"],
            ["influenza", "Influenza template"],
            ["legacy", "Simple dictionary"],
          ] as const
        ).map(([value, label]) => (
          <button
            key={value}
            type="button"
            className="button secondary"
            aria-pressed={kind === value}
            onClick={() => setKind(value)}
          >
            {t(label)}
          </button>
        ))}
      </div>
      <pre tabIndex={0} aria-label={t("Primer database JSON example")}>
        <code>{text}</code>
      </pre>
      <button
        type="button"
        className="button secondary"
        onClick={() =>
          download(text, `${kind}-primers.json`, "application/json")
        }
      >
        <Download size={16} />
        {t("Download JSON template")}
      </button>
      <p className="format-footnote">
        {t(
          "Synthetic primers for showing the format, not a validated assay. Replace the organism, names, and sequences before using your own data.",
        )}
      </p>
      {kind === "legacy" ? (
        <p>
          {t(
            "The simple format maps organism → primer name → sequence. It is treated as PCR and has no explicit scheme, pool, or subtype fields. Use the full template for NGS or influenza metadata.",
          )}
        </p>
      ) : (
        <dl className="format-fields">
          <dt>
            <code>schema_version</code> · <code>database_version</code>
          </dt>
          <dd>
            {t(
              'Keep schema_version as "1.0". Set database_version to your own version label.',
            )}
          </dd>
          <dt>
            <code>schemes[]</code>
          </dt>
          <dd>
            {t(
              "Each scheme needs scheme_id, display_name, organism, version, and primers. Give every scheme a unique scheme_id. Use assay_type: pcr or ngs (pcr is the default). Add another scheme object to include another assay or organism.",
            )}
          </dd>
          <dt>
            <code>primers[]</code>
          </dt>
          <dd>
            {t(
              "Each primer needs id, name, sequence, role, and segment. IDs and names must be distinct within a scheme. Write both forward and reverse primer sequences as ordered, 5′ → 3′, using IUPAC nucleotide letters. Use role values such as forward, reverse, or probe; use an empty segment for non-influenza primers.",
            )}
          </dd>
          <dt>
            <code>segment</code> · <code>subtype_tags</code>
          </dt>
          <dd>
            {t(
              'For influenza, specify the type in organism (for example "Influenza-A", "Influenza-B", "Influenza-C", or "Influenza-D"). Segment labels come from your database, for example "PB2", "PB1", "PA", "HA", "NP", "NA", "M", or "NS", and must match the FASTA headers. subtype_tags accepts labels such as ["H5N1", "H7N9"]. A selected tag matches exactly: "H5" does not automatically include "H5N1". List both tags to include a primer in both selections. An explicit empty list [] makes a primer shared by all subtypes of its influenza type; omitted tags can be inferred from legacy primer names. The subtype menu shows only types and tags present in the database; a type such as A selects all its primers, while other types use labels such as B/VICTORIA.',
            )}
          </dd>
          <dt>
            <code>pool</code> · <code>strand</code>
          </dt>
          <dd>
            {t(
              'Optional: pool labels multiplex groups; strand may be "+", "-", or "". Upload inline primer sequences; files referring to external BED or FASTA assets require the CLI.',
            )}
          </dd>
        </dl>
      )}
      <p>
        {t(
          "Save as UTF-8 JSON with double quotes and no comments or trailing commas. Maximum 250 kB, 500 primers, and 200 bases per primer. The upload check validates your file before it can be used.",
        )}
      </p>
    </div>
  );
  if (expanded) return content;
  return (
    <details className="input-format-help">
      <summary>
        {t("View primer database format and download templates")}
      </summary>
      {content}
    </details>
  );
}
