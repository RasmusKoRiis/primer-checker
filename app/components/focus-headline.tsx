"use client";

import { useEffect, useRef } from "react";

export default function FocusHeadline({
  firstLine,
  secondLine,
}: {
  firstLine: string;
  secondLine: string;
}) {
  const heading = useRef<HTMLHeadingElement>(null);
  const canvas = useRef<HTMLCanvasElement>(null);
  // Both translations use “base”; find the accent without splitting a word for translation.
  const mismatch = firstLine.lastIndexOf("base.") + 1;

  useEffect(() => {
    const title = heading.current;
    const overlay = canvas.current;
    const motion = window.matchMedia("(prefers-reduced-motion: reduce)");
    if (!title || !overlay || motion.matches || document.hidden) return;
    const context = overlay.getContext("2d");
    const source = document.createElement("canvas");
    const sourceContext = source.getContext("2d");
    const pixels = document.createElement("canvas");
    const pixelContext = pixels.getContext("2d");
    if (!context || !sourceContext || !pixelContext) return;

    const bounds = title.getBoundingClientRect();
    if (!bounds.width || !bounds.height) return;
    overlay.width = source.width = Math.ceil(bounds.width);
    overlay.height = source.height = Math.ceil(bounds.height);
    const accent = getComputedStyle(title).getPropertyValue("--accent");
    const lines = Array.from(
      title.querySelectorAll<HTMLElement>(".focus-line"),
    );
    const glyphs = lines.flatMap((line, lineIndex) => {
      const style = getComputedStyle(line);
      const box = line.getBoundingClientRect();
      const text = line.textContent || "";
      const font = `${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
      sourceContext.font = font;
      const metrics = sourceContext.measureText(text);
      const size = parseFloat(style.fontSize);
      const ascent = metrics.fontBoundingBoxAscent ?? size * 0.8;
      const descent = metrics.fontBoundingBoxDescent ?? size * 0.2;
      const baseline =
        box.top - bounds.top + (box.height - ascent - descent) / 2 + ascent;
      const spacing = parseFloat(style.letterSpacing) || 0;
      return Array.from(text, (letter, index) => ({
        letter,
        font,
        color: lineIndex === 0 && index === mismatch ? accent : style.color,
        x:
          box.left -
          bounds.left +
          sourceContext.measureText(text.slice(0, index)).width +
          index * spacing,
        y: baseline,
        mismatch: lineIndex === 0 && index === mismatch,
        width: sourceContext.measureText(letter).width + spacing,
        height: parseFloat(style.fontSize),
      }));
    });

    let frame = 0;
    let lastStage = -1;
    const started = performance.now();
    const finish = () => {
      cancelAnimationFrame(frame);
      delete title.dataset.focusing;
    };
    const draw = (now: number) => {
      const elapsed = now - started;
      if (elapsed >= 1850) {
        finish();
        return;
      }
      const stage = Math.min(6, Math.floor(elapsed / 230));
      if (stage !== lastStage) {
        sourceContext.clearRect(0, 0, source.width, source.height);
        for (const glyph of glyphs) {
          sourceContext.font = glyph.font;
          sourceContext.fillStyle = glyph.color;
          sourceContext.fillText(glyph.letter, glyph.x, glyph.y);
        }
        const block = [10, 8, 6, 4, 3, 2, 1][stage];
        pixels.width = Math.max(1, Math.ceil(source.width / block));
        pixels.height = Math.max(1, Math.ceil(source.height / block));
        pixelContext.drawImage(source, 0, 0, pixels.width, pixels.height);
        context.clearRect(0, 0, overlay.width, overlay.height);
        context.imageSmoothingEnabled = false;
        context.drawImage(pixels, 0, 0, overlay.width, overlay.height);

        // One small region remains coarse until the last pass: the hidden mismatch.
        const hit = glyphs.find((glyph) => glyph.mismatch);
        if (hit && stage >= 3 && stage < 6) {
          const top = Math.max(0, Math.floor(hit.y - hit.height * 0.75));
          const width = Math.ceil(hit.width);
          const height = Math.ceil(hit.height * 0.8);
          pixels.width = Math.max(1, Math.ceil(width / 7));
          pixels.height = Math.max(1, Math.ceil(height / 7));
          pixelContext.drawImage(
            source,
            hit.x,
            top,
            width,
            height,
            0,
            0,
            pixels.width,
            pixels.height,
          );
          context.clearRect(hit.x, top, width, height);
          context.drawImage(pixels, hit.x, top, width, height);
          context.strokeStyle = accent;
          context.lineWidth = 1;
          context.strokeRect(
            Math.floor(hit.x) - 2.5,
            top - 2.5,
            width + 5,
            height + 5,
          );
        }
        lastStage = stage;
      }
      frame = requestAnimationFrame(draw);
    };
    draw(started);
    title.dataset.focusing = "true";
    motion.addEventListener("change", finish);
    window.addEventListener("resize", finish);
    document.addEventListener("visibilitychange", finish);
    return () => {
      finish();
      motion.removeEventListener("change", finish);
      window.removeEventListener("resize", finish);
      document.removeEventListener("visibilitychange", finish);
    };
  }, [firstLine, secondLine, mismatch]);

  return (
    <h1 id="intro-title" className="focus-headline" ref={heading}>
      <span className="focus-line">
        {firstLine.slice(0, mismatch)}
        <span className="focus-mismatch">{firstLine[mismatch]}</span>
        {firstLine.slice(mismatch + 1)}
      </span>
      <span className="focus-line">{secondLine}</span>
      <canvas className="focus-canvas" ref={canvas} aria-hidden="true" />
    </h1>
  );
}
