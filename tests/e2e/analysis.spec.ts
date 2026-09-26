import { expect, test } from "@playwright/test";
import path from "node:path";

test("upload, analyze with real BLAST, inspect mismatch, and download reports", async ({
  page,
}) => {
  await page.goto("/");
  await expect(
    page.getByRole("combobox", { name: "Virus", exact: true }),
  ).toBeEnabled();
  await page
    .getByLabel("Upload FASTA files")
    .setInputFiles(path.resolve("public/example.fasta"));
  await page
    .getByRole("button", { name: "Analyze sequences", exact: true })
    .click();
  await expect(
    page.getByRole("heading", { name: "Compatibility results" }),
  ).toBeVisible({ timeout: 45_000 });
  await expect(page.getByText("6 comparisons", { exact: true })).toBeVisible();
  await page.getByRole("tab", { name: "By sample 6", exact: true }).click();
  await page
    .getByRole("combobox", { name: "Mismatches", exact: true })
    .selectOption("any");
  await expect(
    page.getByRole("button", { name: "Inspect", exact: true }),
  ).toHaveCount(1);
  await page.getByRole("button", { name: "Inspect", exact: true }).click();
  await expect(page.getByRole("dialog")).toBeVisible();
  await expect(
    page.getByRole("dialog").getByText("9:T>A", { exact: true }),
  ).toBeVisible();
  await expect(page.locator(".alignment-column.mismatch")).toHaveCount(1);
  await page.getByRole("button", { name: "Close alignment" }).click();
  for (const [label, filename] of [
    ["CSV", "primer-results.csv"],
    ["HTML report", "primer-report.html"],
  ]) {
    const downloaded = page.waitForEvent("download");
    await page.getByRole("button", { name: label, exact: true }).click();
    expect((await downloaded).suggestedFilename()).toBe(filename);
  }
});

test("mobile upload and influenza configuration remain usable", async ({
  page,
}) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("/");
  await page
    .getByRole("combobox", { name: "Virus", exact: true })
    .selectOption("influenza");
  await page.getByLabel("Influenza subtype").selectOption("H1");
  await expect(page.getByLabel("Influenza subtype")).toHaveValue("H1");
  await expect(
    page.getByRole("button", { name: "Choose files", exact: true }),
  ).toBeVisible();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
});

test("build, download, re-upload, and analyze a custom primer database", async ({
  page,
}) => {
  await page.goto("/");
  await page
    .getByRole("button", { name: "Build a database", exact: true })
    .click();
  await page
    .getByRole("textbox", { name: "Database name", exact: true })
    .fill("My assay");
  await page
    .getByRole("combobox", { name: "Organism", exact: true })
    .fill("Custom-virus");
  await page.getByLabel("Primer 1 name", { exact: true }).fill("target_F");
  await page
    .getByLabel("Primer 1 sequence", { exact: true })
    .fill("ACGTTGCAAGCTTAGCGATCGATGCTAGCA");
  await page
    .getByRole("button", { name: "Create database", exact: true })
    .click();
  await expect(
    page.getByText("custom-primers.json is ready for analysis", {
      exact: true,
    }),
  ).toBeVisible();
  const downloadEvent = page.waitForEvent("download");
  await page
    .getByRole("button", { name: "Download database JSON", exact: true })
    .click();
  const download = await downloadEvent;
  expect(download.suggestedFilename()).toBe("custom-primers.json");
  const saved = await download.path();
  expect(saved).toBeTruthy();
  const { readFile } = await import("node:fs/promises");
  const database = await readFile(saved!, "utf8");
  expect(JSON.parse(database).schemes[0].primers[0].sequence).toBe(
    "ACGTTGCAAGCTTAGCGATCGATGCTAGCA",
  );
  await page.getByRole("button", { name: "Upload JSON", exact: true }).click();
  await page
    .getByLabel("Upload primer database", { exact: true })
    .setInputFiles({
      name: "my-primers.json",
      mimeType: "application/json",
      buffer: Buffer.from(database),
    });
  await expect(
    page.getByText("my-primers.json is ready for analysis", { exact: true }),
  ).toBeVisible();
  await page.getByLabel("Upload FASTA files").setInputFiles({
    name: "custom.fa",
    mimeType: "text/plain",
    buffer: Buffer.from(">sample\nACGTTGCAAGCTTAGCGATCGATGCTAGCA\n"),
  });
  await expect(
    page.getByText(
      "Ready: 1 records · 1 primers · 1 comparisons · 1 BLAST searches",
      { exact: true },
    ),
  ).toBeVisible();
  await page
    .getByRole("button", { name: "Analyze sequences", exact: true })
    .click();
  await expect(
    page.getByRole("heading", { name: "Compatibility results" }),
  ).toBeVisible({ timeout: 45_000 });
  await expect(page.getByText("1 comparisons", { exact: true })).toBeVisible();
  await page.getByRole("tab", { name: "By sample 1", exact: true }).click();
  await page
    .getByRole("combobox", { name: "Mismatches", exact: true })
    .selectOption("0");
  await expect(
    page.getByRole("button", { name: "Inspect", exact: true }),
  ).toHaveCount(1);
  await expect(
    page.getByRole("cell", { name: "target_F My assay", exact: true }),
  ).toBeVisible();
});

test("oversized uploads and excessive records cannot start analysis", async ({
  page,
}) => {
  await page.goto("/");
  await expect(
    page.getByRole("combobox", { name: "Virus", exact: true }),
  ).toBeEnabled();
  const requests: string[] = [];
  page.on("request", (request) => {
    if (/\/api\/(preflight|analyze)$/.test(request.url()))
      requests.push(request.url());
  });
  await page.getByLabel("Upload FASTA files").setInputFiles({
    name: "large.fa",
    mimeType: "text/plain",
    buffer: Buffer.alloc(3_000_001, "A"),
  });
  await expect(
    page.getByRole("alert").filter({ hasText: "Combined uploads" }),
  ).toContainText("Combined uploads exceed 3 MB");
  await expect(
    page.getByRole("button", { name: "Analyze sequences", exact: true }),
  ).toBeDisabled();
  expect(requests).toHaveLength(0);
  await page.getByLabel("Upload FASTA files").setInputFiles({
    name: "records.fa",
    mimeType: "text/plain",
    buffer: Buffer.from(
      Array.from({ length: 201 }, (_, i) => `>sample${i}\nACGT\n`).join(""),
    ),
  });
  await expect(page.getByLabel("Batch limits")).toContainText(
    "Use at most 200 sequence records",
  );
  await expect(
    page.getByRole("button", { name: "Analyze sequences", exact: true }),
  ).toBeDisabled();
  expect(requests.some((url) => url.endsWith("/api/analyze"))).toBe(false);
});
