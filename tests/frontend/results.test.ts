import { describe, expect, it } from "vitest";
import { aggregatePrimers, emptyFilters, filterRows, validateFiles } from "../../app/lib/results";
import type { ResultRow } from "../../app/lib/results";
const base = { Primer_Name: "primer", Primer_Sequence: "ACTG", Assay_ID: "scheme", Assay_Name: "Scheme", Primer_Segment: "HA", Subject_Segment: "HA", Primer_Pool: "1", Subject_Sequence_ID: "sample", Hit_Status: "hit", Mismatches: 0 } as ResultRow;
const rows = [base, { ...base, Mismatches: 2 }, { ...base, Hit_Status: "no_hit", Mismatches: "" }, { ...base, Assay_ID: "other", Mismatches: 1 }];
describe("result summaries", () => {
  it("keeps no-hits separate and distinct assays separate", () => {
    const [first, second] = aggregatePrimers(rows);
    expect(first).toMatchObject({ tested: 3, perfect: 1, affected: 1, noHit: 1, maximum: 2 });
    expect(first.percent).toBeCloseTo(100 / 3);
    expect(second.tested).toBe(1);
  });
  it("does not classify a no-hit as a zero-mismatch result", () => {
    expect(filterRows(rows, { ...emptyFilters, mismatches: "0" })).toEqual([base]);
    expect(filterRows(rows, { ...emptyFilters, mismatches: "any", assay: "scheme", segment: "HA" })).toHaveLength(1);
    expect(filterRows(rows, { ...emptyFilters, status: "no_hit" })).toHaveLength(1);
  });
  it("filters sample and primer case-insensitively", () => {
    expect(filterRows(rows, { ...emptyFilters, primer: "PRIM", sample: "SAMP" })).toHaveLength(4);
    expect(filterRows(rows, { ...emptyFilters, sample: "missing" })).toHaveLength(0);
  });
});
describe("upload validation", () => {
  const limits = { upload_bytes: 3_000_000, files: 10, records: 200, comparisons: 2000 };
  it("bounds the combined metadata and FASTA payload", () => {
    expect(validateFiles([{ name: "a.fasta", size: 2_000_000 }], { name: "m.csv", size: 1_000_001 }, limits)).toMatch(/3 MB/);
    expect(validateFiles([{ name: "a.fasta", size: 100 }], null, limits)).toBeNull();
  });
  it("rejects unsupported, empty, and duplicate files", () => {
    expect(validateFiles([{ name: "a.fastq", size: 100 }], null, limits)).toMatch(/FASTQ/);
    expect(validateFiles([{ name: "a.fa", size: 0 }], null, limits)).toMatch(/empty/);
    expect(validateFiles([{ name: "a.fa", size: 5 }, { name: "a.fa", size: 6 }], null, limits)).toMatch(/distinct/);
  });
});
