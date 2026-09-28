"use client";

import { useEffect, useRef } from "react";

function CyclingBase({ letter }: { letter: string }) {
  const slot = useRef<HTMLSpanElement>(null);
  const canvas = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const target = slot.current;
    const overlay = canvas.current;
    if (!target || !overlay) return;
    const context = overlay.getContext("2d");
    const source = document.createElement("canvas");
    const sourceContext = source.getContext("2d");
    const pixels = document.createElement("canvas");
    const pixelContext = pixels.getContext("2d");
    if (!context || !sourceContext || !pixelContext) return;

    const motion = window.matchMedia("(prefers-reduced-motion: reduce)");
    const letters = [letter, "c", "g", "t"];
    let timer = 0;
    let inView = false;
    const stop = () => {
      window.clearTimeout(timer);
      delete target.dataset.cycling;
    };
    const start = () => {
      stop();
      if (!inView || motion.matches || document.hidden) return;
      const box = target.getBoundingClientRect();
      const overlayBox = overlay.getBoundingClientRect();
      if (!box.width || !box.height) return;
      overlay.width = source.width = Math.ceil(overlayBox.width);
      overlay.height = source.height = Math.ceil(overlayBox.height);
      const style = getComputedStyle(target);
      const font = `${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
      const mismatchColor = style.getPropertyValue("--accent").trim();
      sourceContext.font = font;
      const metrics = sourceContext.measureText(letter);
      const size = parseFloat(style.fontSize);
      const ascent = metrics.fontBoundingBoxAscent ?? size * 0.8;
      const descent = metrics.fontBoundingBoxDescent ?? size * 0.2;
      const baseline =
        box.top - overlayBox.top + (box.height - ascent - descent) / 2 + ascent;
      const started = performance.now();
      let lastFrame = "";
      const draw = () => {
        const elapsed = performance.now() - started;
        const phase = (elapsed % 1600) / 1600;
        const index = Math.floor(elapsed / 1600);
        const glyph =
          letters[(index + (phase >= 0.65 ? 1 : 0)) % letters.length];
        const block =
          phase < 0.35
            ? 1
            : phase < 0.5
              ? 3
              : phase < 0.65
                ? 7
                : phase < 0.8
                  ? 4
                  : phase < 0.9
                    ? 2
                    : 1;
        const frame = `${glyph}:${block}`;
        if (frame !== lastFrame) {
          sourceContext.clearRect(0, 0, source.width, source.height);
          sourceContext.fillStyle =
            glyph === letter ? style.color : mismatchColor;
          const glyphMetrics = sourceContext.measureText(glyph);
          const inkWidth =
            glyphMetrics.actualBoundingBoxLeft +
            glyphMetrics.actualBoundingBoxRight;
          const x =
            (source.width - inkWidth) / 2 + glyphMetrics.actualBoundingBoxLeft;
          sourceContext.fillText(glyph, x, baseline);
          pixels.width = Math.max(1, Math.ceil(source.width / block));
          pixels.height = Math.max(1, Math.ceil(source.height / block));
          pixelContext.drawImage(source, 0, 0, pixels.width, pixels.height);
          context.clearRect(0, 0, overlay.width, overlay.height);
          context.imageSmoothingEnabled = false;
          context.drawImage(pixels, 0, 0, overlay.width, overlay.height);
          lastFrame = frame;
        }
        // The stepped effect only needs a few updates per second.
        timer = window.setTimeout(draw, 120);
      };
      draw();
      target.dataset.cycling = "true";
    };
    const visibility = new IntersectionObserver(([entry]) => {
      inView = entry.isIntersecting;
      start();
    });
    const size = new ResizeObserver(start);
    visibility.observe(target);
    size.observe(target);
    motion.addEventListener("change", start);
    document.addEventListener("visibilitychange", start);
    return () => {
      stop();
      visibility.disconnect();
      size.disconnect();
      motion.removeEventListener("change", start);
      document.removeEventListener("visibilitychange", start);
    };
  }, [letter]);

  return (
    <span className="focus-mismatch" ref={slot}>
      <span className="focus-letter">{letter}</span>
      <canvas className="focus-canvas" ref={canvas} aria-hidden="true" />
    </span>
  );
}

export default function FocusHeadline({
  firstLine,
  secondLine,
}: {
  firstLine: string;
  secondLine: string;
}) {
  // Both translations use “base”; keep the full sentence in the translation file.
  const base = firstLine.lastIndexOf("base.");
  const mismatch = base + 1;
  return (
    <h1
      id="intro-title"
      className="focus-headline"
      aria-label={`${firstLine} ${secondLine}`}
    >
      <span className="focus-line">
        {base < 0 ? (
          firstLine
        ) : (
          <>
            {firstLine.slice(0, mismatch)}
            <CyclingBase letter={firstLine[mismatch]} />
            {firstLine.slice(mismatch + 1)}
          </>
        )}
      </span>
      <span className="focus-line">{secondLine}</span>
    </h1>
  );
}
