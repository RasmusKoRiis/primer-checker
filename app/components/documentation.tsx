"use client";
import english from "../../report_text/english.json";
import norwegian from "../../report_text/norwegian.json";
import { useLanguage } from "./language";
import { DatabaseFormatHelp, FastaFormatHelp } from "./input-format-help";

type ReportKey = keyof typeof english;
export default function Documentation() {
  const { language, t } = useLanguage();
  const report = language === "no" ? norwegian : english;
  const list = (keys: ReportKey[], ordered = false) => {
    const items = keys.map((key) => <li key={key}>{report[key]}</li>);
    return ordered ? <ol>{items}</ol> : <ul>{items}</ul>;
  };
  return (
    <article className="documentation-view">
      <div className="page-heading">
        <div>
          <div className="eyebrow">{t("METHODS & REPORT GUIDE")}</div>
          <h1>{report.doc_title}</h1>
          <p>{report.doc_intro}</p>
        </div>
      </div>
      <nav
        className="documentation-contents"
        aria-label={t("Documentation sections")}
      >
        <a href="#documentation-fasta">{t("FASTA headers")}</a>
        <a href="#documentation-database">{t("Primer database format")}</a>
        <a href="#documentation-method">{report.doc_workflow_title}</a>
        <a href="#documentation-web">{t("Website results")}</a>
        <a href="#documentation-html">{t("Downloaded HTML report")}</a>
        <a href="#documentation-levels">{report.doc_risk_title}</a>
        <a href="#documentation-limits">{report.doc_limits_title}</a>
      </nav>
      <section className="card documentation-section" id="documentation-fasta">
        <h2>{t("FASTA headers")}</h2>
        <FastaFormatHelp expanded />
      </section>
      <section
        className="card documentation-section"
        id="documentation-database"
      >
        <h2>{t("Primer database format")}</h2>
        <DatabaseFormatHelp expanded />
      </section>
      <section className="card documentation-section" id="documentation-method">
        <h2>{report.doc_workflow_title}</h2>
        <p>
          {t(
            "The website and CLI use the same Python analysis engine. Your selected database, virus, subtype, assay type, and scheme determine which primers are tested. Influenza primers with a segment label are reported only against matching segment labels in FASTA headers, such as 01-HA|sample.",
          )}
        </p>
        {list(
          [
            "doc_workflow_1",
            "doc_workflow_2",
            "doc_workflow_3",
            "doc_workflow_4",
            "doc_workflow_5",
          ],
          true,
        )}
        <div className="documentation-formula">
          <strong>{t("Identity calculation")}</strong>
          <p>{report.doc_faq_4}</p>
          <code>{t("20 primer bases, 2 mismatches → 90% identity")}</code>
        </div>
        <p>
          {t(
            "BLASTn uses reward 2, penalty −3, word size 4, and DUST filtering. A no-hit result means BLAST did not return a usable alignment under these search settings; the web interface does not add a percentage-identity cutoff. The provenance download records the actual BLAST version and parameters.",
          )}
        </p>
        <h3>{report.doc_glossary_title}</h3>
        {list([
          "doc_term_hit",
          "doc_term_no_hit",
          "doc_term_mismatch",
          "doc_term_identity",
          "doc_term_terminal",
        ])}
      </section>
      <section className="card documentation-section" id="documentation-web">
        <h2>{t("Website results")}</h2>
        <p>
          {t(
            "Start with the summary cards, then use By primer and By sample to examine individual comparisons. Click a primer to narrow the sample table and choose Inspect to see its alignment.",
          )}
        </p>
        <dl className="documentation-definitions">
          <dt>{t("Summary cards")}</dt>
          <dd>
            {t(
              "The top cards describe the complete analysis and stay unchanged when table filters are applied. Files and sequence records count the uploaded inputs. Each combination of filename and sequence ID is counted separately. Successful hits exclude no-hit comparisons. Hits with mismatches counts hit comparisons with at least one mismatch.",
            )}
          </dd>
          <dt>{t("By primer")}</dt>
          <dd>
            {t(
              "Tested is the number of visible primer/sequence comparisons. Perfect means a hit with zero mismatches. Mismatches counts hits with at least one mismatch; No hit is separate. Max. mismatches is the highest mismatch count among hits. Affected is 100 × hits with mismatches / all tested comparisons, including no-hit comparisons in the denominator. Primers from different assays are kept separate.",
            )}
          </dd>
          <dt>{t("By sample")}</dt>
          <dd>
            {t(
              "Each row is one primer checked against one relevant sequence. Identity and mismatch counts are shown only for hits. Positions are numbered from the primer’s 5′ end. A change such as 9:T>A means primer base T corresponds to sequence base A at position 9.",
            )}
          </dd>
          <dt>{t("Filters and downloads")}</dt>
          <dd>
            {t(
              "Website filters change the tables and their row counts. CSV, HTML, and provenance downloads contain the complete analysis, not only the currently filtered table. The HTML report has its own filters, charts, and automatic review levels.",
            )}
          </dd>
          <dt>{t("Alignment")}</dt>
          <dd>
            {t(
              "The inspector shows the primer above the observed sequence, with mismatches highlighted and positions relative to the primer. Forward or reverse describes the subject alignment orientation. No-hit comparisons have no alignment and must not be read as zero mismatches.",
            )}
          </dd>
        </dl>
        <p className="documentation-note">
          {t(
            "The website tables show descriptive counts. The Low, Watch, High, and Critical review levels described below appear in the downloaded HTML report.",
          )}
        </p>
      </section>
      <section className="card documentation-section" id="documentation-html">
        <h2>{t("Downloaded HTML report")}</h2>
        <p>
          {t(
            "Download HTML report to open a self-contained report in your browser. It includes its own English/Norsk switch and Help and documentation tab. The guidance below uses the same text files as that report.",
          )}
        </p>
        <h3>{report.doc_use_title}</h3>
        {list(
          [
            "doc_use_1",
            "doc_use_2",
            "doc_use_3",
            "doc_use_4",
            "doc_use_5",
            "doc_use_6",
          ],
          true,
        )}
        <h3>{report.doc_cards_title}</h3>
        <p>{report.doc_filters_text}</p>
        {list([
          "doc_cards_1",
          "doc_cards_2",
          "doc_cards_3",
          "doc_cards_4",
          "doc_cards_5",
        ])}
        <h3>{report.doc_overview_title}</h3>
        <p>{report.doc_overview_text}</p>
        {list([
          "doc_overview_1",
          "doc_overview_2",
          "doc_overview_3",
          "doc_overview_4",
          "doc_overview_5",
          "doc_overview_6",
          "doc_overview_7",
        ])}
        <h3>{report.doc_panels_title}</h3>
        <p>{report.doc_panels_text}</p>
        <h3>{report.doc_chart_title}</h3>
        <p>{report.doc_chart_text}</p>
        <h3>{report.doc_distribution_title}</h3>
        <p>{report.doc_distribution_text}</p>
        <h3>{report.doc_detail_title}</h3>
        <p>{report.doc_detail_text}</p>
        <h3>{report.doc_alignment_title}</h3>
        <p>{report.doc_alignment_text}</p>
        <h3>{report.doc_ngs_title}</h3>
        <p>{report.doc_ngs_text}</p>
      </section>
      <section className="card documentation-section" id="documentation-levels">
        <h2>{report.doc_risk_title}</h2>
        <p>{report.doc_risk_text}</p>
        {list([
          "doc_risk_critical",
          "doc_risk_high",
          "doc_risk_watch",
          "doc_risk_low",
        ])}
        <p>{report.doc_faq_5}</p>
        <p>{report.doc_faq_6}</p>
        <h3>{report.doc_faq_title}</h3>
        {([1, 2, 3, 4, 5, 6, 7] as const).map((n) => (
          <details key={n}>
            <summary>{report[`doc_faq_${n}_q`]}</summary>
            <p>{report[`doc_faq_${n}`]}</p>
          </details>
        ))}
      </section>
      <section className="card documentation-section" id="documentation-limits">
        <h2>{report.doc_limits_title}</h2>
        {list([
          "doc_limits_1",
          "doc_limits_2",
          "doc_limits_3",
          "doc_limits_4",
          "doc_limits_5",
        ])}
        <h3>{t("Uploads and reproducibility")}</h3>
        <p>
          {t(
            "The bundled dummy database contains invented sequences for software testing only. Use the synthetic FASTA example to try it. For your own analysis, upload a self-contained JSON library or build and download your own database. The preflight check validates file sizes and workload before Analyze is enabled; the server repeats those checks during analysis. Large workloads should use the CLI.",
          )}
        </p>
        <p>
          {t(
            "Uploads are processed temporarily. Download your custom database and the CSV, HTML, and provenance JSON to retain a reproducible record. Provenance includes database and input fingerprints, selections, time, application version, and BLAST settings. Language selection changes presentation only; primer names, sequences, calculations, and CSV column names remain unchanged.",
          )}
        </p>
      </section>
      <a className="button secondary" href="#analysis">
        {t("Back to analysis")}
      </a>
    </article>
  );
}
