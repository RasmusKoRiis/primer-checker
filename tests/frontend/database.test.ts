import { expect, it } from "vitest";
import { makeDatabase, type DatabaseDraft } from "../../app/lib/database";
import { validateFiles } from "../../app/lib/results";

const draft: DatabaseDraft = {
  name: "My NGS primers",
  organism: "Influenza-A",
  version: "1",
  assayType: "ngs",
  primers: [
    {
      key: 1,
      name: "HA_F",
      sequence: "acgt ryn\n",
      role: "forward",
      segment: "HA",
      pool: "1",
      subtype: "H3",
    },
  ],
};
it("exports a portable scheme with normalized bases and selection metadata", () => {
  const db = JSON.parse(makeDatabase(draft));
  expect(db.schemes[0]).toMatchObject({
    assay_type: "ngs",
    organism: "Influenza-A",
  });
  expect(db.schemes[0].primers[0]).toMatchObject({
    sequence: "ACGTRYN",
    segment: "HA",
    pool: "1",
    subtype_tags: ["H3"],
  });
});
it("rejects incomplete, invalid, and ambiguous primer entries", () => {
  expect(() => makeDatabase({ ...draft, name: "" })).toThrow(/name/);
  expect(() =>
    makeDatabase({ ...draft, primers: [draft.primers[0], draft.primers[0]] }),
  ).toThrow(/distinct/);
  expect(() =>
    makeDatabase({
      ...draft,
      primers: [{ ...draft.primers[0], sequence: "ACGU-" }],
    }),
  ).toThrow(/IUPAC/);
});
it("includes the custom database in the combined browser upload limit", () => {
  expect(
    validateFiles(
      [{ name: "sample.fa", size: 2_900_000 }],
      { files: 10, upload_bytes: 3_000_000, records: 200, comparisons: 2000 },
      { name: "custom.json", size: 100_001 },
    ),
  ).toMatch(/3 MB/);
});
