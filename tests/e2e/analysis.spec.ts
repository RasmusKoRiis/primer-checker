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
  await page.getByRole("button", { name: "Norsk", exact: true }).click();
  await expect(
    page.getByRole("heading", {
      name: "Kompatibilitetsresultater",
      exact: true,
    }),
  ).toBeVisible();
  await expect(
    page.getByRole("combobox", { name: "Mismatch", exact: true }),
  ).toHaveValue("any");
  await page.getByRole("button", { name: "Undersøk", exact: true }).click();
  await expect(
    page.getByRole("dialog").getByText("9:T>A", { exact: true }),
  ).toBeVisible();
  await page
    .getByRole("button", { name: "Lukk sekvenssammenstillingen", exact: true })
    .click();
  await page.getByRole("button", { name: "English", exact: true }).click();

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

test("documentation and language switching preserve the analysis and database draft", async ({
  page,
}) => {
  await page.goto("/");
  await page
    .getByRole("button", { name: "Use a synthetic example", exact: true })
    .click();
  await page
    .getByRole("button", { name: "Build a database", exact: true })
    .click();
  await page
    .getByRole("textbox", { name: "Database name", exact: true })
    .fill("My preserved assay");
  await page.getByLabel("Primer 1 name", { exact: true }).fill("My_F");
  await page.getByRole("link", { name: "Documentation", exact: true }).click();
  await expect(
    page.getByRole("heading", { name: "Help and documentation", exact: true }),
  ).toBeVisible();
  await expect(
    page.getByText("20 primer bases, 2 mismatches → 90% identity", {
      exact: true,
    }),
  ).toBeVisible();
  await page
    .getByRole("link", { name: "Automatic review levels", exact: true })
    .click();
  await expect(
    page.getByRole("heading", { name: "Automatic review levels", exact: true }),
  ).toBeVisible();
  await expect(
    page.getByText(/Critical: no-hit rate is at least 20%/),
  ).toBeVisible();
  await page.getByRole("button", { name: "Norsk", exact: true }).click();
  await expect(page.locator("html")).toHaveAttribute("lang", "nb");
  await expect(
    page.getByRole("heading", { name: "Hjelp og dokumentasjon", exact: true }),
  ).toBeVisible();
  await expect(
    page.getByText(/Kritisk: andelen uten treff er minst 20/),
  ).toBeVisible();
  await page.getByRole("link", { name: "Analyse", exact: true }).click();
  await expect(
    page.getByRole("textbox", { name: "Databasenavn", exact: true }),
  ).toHaveValue("My preserved assay");
  await expect(page.getByLabel("Primer 1 navn", { exact: true })).toHaveValue(
    "My_F",
  );
  await expect(
    page.getByRole("button", { name: "Fjern example.fasta", exact: true }),
  ).toBeVisible();
  await page.reload();
  await expect(
    page.getByRole("heading", { name: "Ny analyse", exact: true }),
  ).toBeVisible();
  await expect(
    page.getByRole("button", { name: "Norsk", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await page.getByLabel("Last opp FASTA-filer").setInputFiles({
    name: "too-big.fa",
    mimeType: "text/plain",
    buffer: Buffer.alloc(3_000_001, "A"),
  });
  await expect(
    page.getByRole("alert").filter({ hasText: "overskrider totalt 3 MB" }),
  ).toBeVisible();
  await page.getByRole("button", { name: "English", exact: true }).click();
  await expect(page.locator("html")).toHaveAttribute("lang", "en");
  await expect(
    page.getByRole("heading", { name: "New analysis", exact: true }),
  ).toBeVisible();
});

test("Norwegian documentation opens directly and remains usable on mobile", async ({
  page,
}) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("/#documentation-method");
  await expect(
    page.getByRole("heading", {
      name: "How the values are calculated",
      exact: true,
    }),
  ).toBeVisible();
  await page.getByRole("button", { name: "Norsk", exact: true }).click();
  await expect(
    page.getByRole("heading", { name: "Slik beregnes verdiene", exact: true }),
  ).toBeVisible();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
  await page.getByRole("link", { name: "Analyse", exact: true }).click();
  await expect(
    page.getByRole("button", { name: "Velg filer", exact: true }),
  ).toBeVisible();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
});

test("format guides download usable templates and explain influenza headers in both languages", async ({
  page,
}) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("/");
  const sequenceCard = page.getByRole("region", {
    name: "Sequence files",
    exact: true,
  });
  const databaseCard = page.getByRole("region", {
    name: "Primer database",
    exact: true,
  });
  await expect(
    sequenceCard.getByText(">sample_001", { exact: true }),
  ).toBeVisible();
  await sequenceCard
    .getByText("FASTA header rules and examples", { exact: true })
    .click();
  await sequenceCard
    .getByRole("button", { name: "Influenza FASTA", exact: true })
    .click();
  await expect(
    sequenceCard.getByLabel("FASTA example", { exact: true }),
  ).toContainText(">01-HA|sample_001");
  const fastaEvent = page.waitForEvent("download");
  await sequenceCard
    .getByRole("button", { name: "Download FASTA template", exact: true })
    .click();
  const fastaDownload = await fastaEvent;
  expect(fastaDownload.suggestedFilename()).toBe("influenza-template.fasta");

  await databaseCard
    .getByText("View primer database format and download templates", {
      exact: true,
    })
    .click();
  await databaseCard
    .getByRole("button", { name: "Influenza template", exact: true })
    .click();
  const databaseEvent = page.waitForEvent("download");
  await databaseCard
    .getByRole("button", { name: "Download JSON template", exact: true })
    .click();
  const databaseDownload = await databaseEvent;
  expect(databaseDownload.suggestedFilename()).toBe("influenza-primers.json");
  const { readFile } = await import("node:fs/promises");
  const database = await readFile((await databaseDownload.path())!, "utf8");
  expect(JSON.parse(database).schemes[0].organism).toBe("Influenza-A");
  await page.getByRole("button", { name: "Upload JSON", exact: true }).click();
  await page
    .getByLabel("Upload primer database", { exact: true })
    .setInputFiles({
      name: databaseDownload.suggestedFilename(),
      mimeType: "application/json",
      buffer: Buffer.from(database),
    });
  await expect(
    page.getByText("influenza-primers.json is ready for analysis", {
      exact: true,
    }),
  ).toBeVisible();
  await page.getByLabel("Upload FASTA files").setInputFiles({
    name: fastaDownload.suggestedFilename(),
    mimeType: "text/plain",
    buffer: await readFile((await fastaDownload.path())!),
  });
  await page
    .getByRole("combobox", { name: "Influenza subtype", exact: true })
    .selectOption("H1");
  await expect(
    page.getByRole("button", { name: "Analyze sequences", exact: true }),
  ).toBeEnabled();
  await expect(page.getByLabel("Batch limits")).toContainText(
    "3 records · 2 primers · 2 comparisons",
  );
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
  await page.getByRole("button", { name: "Norsk", exact: true }).click();
  await expect(
    page.getByText("Influensa trenger segmentmerking", { exact: true }),
  ).toBeVisible();
  await expect(
    page.getByRole("button", { name: "Last ned JSON-mal", exact: true }),
  ).toBeVisible();
  await page.getByRole("link", { name: "Dokumentasjon", exact: true }).click();
  await page
    .getByRole("link", { name: "Format for primerdatabaser", exact: true })
    .click();
  await expect(
    page.getByRole("heading", {
      name: "Format for primerdatabaser",
      exact: true,
    }),
  ).toBeVisible();
  await page
    .getByRole("button", { name: "Enkelt oppslagsformat", exact: true })
    .click();
  await expect(
    page
      .locator("#documentation")
      .getByLabel("JSON-eksempel på primerdatabase", { exact: true }),
  ).toContainText('"Example-virus"');
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
});
