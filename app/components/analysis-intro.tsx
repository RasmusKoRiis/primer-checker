"use client";

import { useLanguage } from "./language";
import FocusHeadline from "./focus-headline";

export default function AnalysisIntro() {
  const { t } = useLanguage();
  return (
    <section className="analysis-intro" aria-labelledby="intro-title">
      <div className="intro-copy">
        <div className="intro-kicker">
          <span />
          {t("CONSENSUS SEQUENCE ANALYSIS")}
        </div>
        <FocusHeadline
          firstLine={t("Every base.")}
          secondLine={t("In focus.")}
        />
        <div className="intro-summary">
          <p>
            {t(
              "Evaluate PCR and sequencing primer compatibility against viral consensus sequences.",
            )}
          </p>
        </div>
      </div>
    </section>
  );
}
