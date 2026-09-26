import { expect, it } from "vitest";
import norwegian from "../../app/lib/norwegian.json";
import { translate } from "../../app/lib/translate";

it("keeps translated placeholders intact so names and counts are never lost", () => {
  const placeholders = (text: string) =>
    [...text.matchAll(/\{(\w+)\}/g)].map((m) => m[1]).sort();
  for (const [english, translated] of Object.entries(norwegian)) {
    expect(translated.trim(), english).not.toBe("");
    expect(placeholders(translated), english).toEqual(placeholders(english));
  }
  expect(
    translate("no", "{filename} is ready for analysis", {
      filename: "My assay.json",
    }),
  ).toBe("My assay.json er klar til analyse");
});
it("localizes dynamic validation messages and preserves unrecognized technical text", () => {
  expect(translate("no", "Choose at most 10 FASTA files.")).toBe(
    "Velg maksimalt 10 FASTA-filer.",
  );
  expect(
    translate("no", "Use at most 200 sequence records per analysis."),
  ).toBe("Bruk maksimalt 200 sekvensoppføringer per analyse.");
  expect(translate("no", "unrecognized technical detail")).toBe(
    "unrecognized technical detail",
  );
  expect(translate("en", "{count} comparisons", { count: 6 })).toBe(
    "6 comparisons",
  );
});
