"use client";

import { ArrowDown, ArrowUpRight } from "./pixel-icons";
import { useLanguage } from "./language";

export default function AnalysisIntro() {
  const { t } = useLanguage();
  return (
    <section className="analysis-intro" aria-labelledby="intro-title">
      <div className="intro-copy">
        <div className="intro-kicker">
          <span />
          {t("CONSENSUS SEQUENCE ANALYSIS")}
        </div>
        <h1 id="intro-title">
          <span>{t("Every base.")}</span>
          <span>{t("In focus.")}</span>
        </h1>
        <div className="intro-summary">
          <p>
            {t(
              "Evaluate PCR and sequencing primer compatibility against viral consensus sequences.",
            )}
          </p>
          <div className="intro-links">
            <a className="intro-start" href="#analysis-workspace">
              {t("Start an analysis")}
              <ArrowDown size={16} />
            </a>
            <a className="intro-guide" href="#documentation-method">
              {t("Explore the method")}
              <ArrowUpRight size={16} />
            </a>
          </div>
        </div>
      </div>
    </section>
  );
}
