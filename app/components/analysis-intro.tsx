"use client";

import { ArrowDown, ArrowUpRight } from "lucide-react";
import { useLanguage } from "./language";

// A decorative study of nucleotide typography in perspective, not analysis data.
function SequenceStudy() {
  const center = { x: 290, y: 155 };
  const depths = [0, 110, 255, 440, 700, 1050];
  return (
    <svg viewBox="0 0 580 310" fill="none" aria-hidden="true" focusable="false">
      <g className="sequence-guides" stroke="currentColor" strokeWidth="0.7">
        <path d="M8 10 290 155 572 10M8 300 290 155 572 300" />
        {depths.map((depth) => {
          const scale = 1 / (1 + depth / 400);
          return (
            <rect
              key={depth}
              x={290 - 282 * scale}
              y={155 - 145 * scale}
              width={564 * scale}
              height={290 * scale}
            />
          );
        })}
      </g>
      <g
        className="sequence-type"
        fill="currentColor"
        textAnchor="middle"
        dominantBaseline="central"
      >
        {depths.map((depth, layer) => {
          const scale = 1 / (1 + depth / 400);
          const project = (x: number, y: number) =>
            `${center.x + (x - center.x) * scale} ${center.y + (y - center.y) * scale}`;
          return (
            <g key={depth} opacity={1 - layer * 0.1}>
              {Array.from("ACGTACGT").map((base, i) => (
                <g key={i}>
                  <text
                    transform={`translate(${project(47 + i * 69, 40)}) scale(${scale} ${scale * 0.65}) rotate(180)`}
                  >
                    {base}
                  </text>
                  <text
                    className={
                      layer === 1 && i === 5 ? "sequence-accent" : undefined
                    }
                    transform={`translate(${project(47 + i * 69, 270)}) scale(${scale} ${scale * 0.65})`}
                  >
                    {base}
                  </text>
                </g>
              ))}
              {Array.from("AGT").map((base, i) => (
                <g key={i}>
                  <text
                    transform={`translate(${project(38, 86 + i * 69)}) scale(${scale * 0.65} ${scale}) rotate(90)`}
                  >
                    {base}
                  </text>
                  <text
                    transform={`translate(${project(542, 86 + i * 69)}) scale(${scale * 0.65} ${scale}) rotate(-90)`}
                  >
                    {base}
                  </text>
                </g>
              ))}
            </g>
          );
        })}
      </g>
      <path
        className="sequence-center"
        d="M282 155h16m-8-8v16"
        stroke="currentColor"
        strokeWidth="1"
      />
    </svg>
  );
}

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
            <ArrowUpRight size={15} />
          </a>
        </div>
      </div>
      <div className="intro-study">
        <div className="study-topline" aria-hidden="true">
          <span>5′ — A C G T — 3′</span>
          <span>01 / PC</span>
        </div>
        <div className="study-art">
          <SequenceStudy />
        </div>
        <div className="study-caption">
          <span>{t("A study in sequence")}</span>
          <span>{t("Illustration")}</span>
        </div>
      </div>
    </section>
  );
}
