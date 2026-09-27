import { expect, test } from "@playwright/test";
import path from "node:path";
import { pathToFileURL } from "node:url";

test("reverse oligos, terminal bases, and gaps agree in the website and downloaded report", async ({
  page,
  context,
}, testInfo) => {
  await page.goto("/");
  await expect(
    page.getByRole("combobox", { name: "Virus", exact: true }),
  ).toBeEnabled();
  await page.getByRole("button", { name: "Upload JSON", exact: true }).click();
  await page
    .getByLabel("Upload primer database", { exact: true })
    .setInputFiles(path.resolve("fixtures/reverse_primer/primers.json"));
  await expect(
    page.getByText("primers.json is ready for analysis", { exact: true }),
  ).toBeVisible();
  await page
    .getByLabel("Upload FASTA files")
    .setInputFiles(path.resolve("fixtures/reverse_primer/sequences.fasta"));
  await expect(page.getByLabel("Batch limits")).toContainText(
    "14 records · 1 primers · 14 comparisons",
  );
  await page
    .getByRole("button", { name: "Analyze sequences", exact: true })
    .click();
  await expect(
    page.getByRole("heading", { name: "Compatibility results", exact: true }),
  ).toBeVisible();
  await page.getByRole("tab", { name: "By sample 14", exact: true }).click();

  for (const [sample, details, highlights, coordinates] of [
    ["perfect_minus", "None", 0, "80–41"],
    ["first_minus", "1:A>T", 1, "79–41"],
    ["last_minus", "40:A>C", 1, "80–42"],
    ["both_minus", "1:A>T,40:A>C", 2, "79–42"],
    ["deletion_minus", "19:T>-", 1, "79–41"],
    ["insertion_minus", "18:->G", 1, "81–41"],
  ] as const) {
    await page
      .getByRole("textbox", { name: "Sample", exact: true })
      .fill(sample);
    await page.getByRole("button", { name: "Inspect", exact: true }).click();
    const dialog = page.getByRole("dialog");
    await expect(dialog).toContainText(
      `Local BLAST hit: Reverse strand · ${coordinates} (1-based, inclusive)`,
    );
    await expect(dialog.locator(".alignment-details")).toContainText(details);
    await expect(dialog).toContainText(
      "Both rows follow the primer’s 5′ → 3′ direction.",
    );
    await expect(dialog.locator(".alignment-column.mismatch")).toHaveCount(
      highlights,
    );
    await expect(dialog.locator(".alignment-column").last()).toContainText(
      "40",
    );
    if (sample === "last_minus") {
      await expect(dialog.locator(".alignment-column").last()).toHaveClass(
        /mismatch/,
      );
    }
    if (sample === "insertion_minus") {
      await expect(dialog.locator(".alignment-column")).toHaveCount(41);
      await expect(dialog.locator(".alignment-column.mismatch")).toHaveText(
        "·-G",
      );
      await page.screenshot({
        path: testInfo.outputPath("reverse-insertion-website.png"),
      });
    }
    await page
      .getByRole("button", { name: "Close alignment", exact: true })
      .click();
  }
  await page.getByRole("button", { name: "Norsk", exact: true }).click();
  await page.getByRole("button", { name: "Undersøk", exact: true }).click();
  await expect(page.getByRole("dialog")).toContainText(
    "Lokalt BLAST-treff: Revers tråd · 81–41",
  );
  await expect(page.getByRole("dialog")).toContainText(
    "Begge rader følger primerens 5′ → 3′-retning.",
  );
  await page
    .getByRole("button", { name: "Lukk sekvenssammenstillingen", exact: true })
    .click();
  await page.getByRole("button", { name: "English", exact: true }).click();

  const downloaded = page.waitForEvent("download");
  await page.getByRole("button", { name: "HTML report", exact: true }).click();
  const reportPath = testInfo.outputPath("reverse-primer-report.html");
  await (await downloaded).saveAs(reportPath);
  const report = await context.newPage();
  await report.goto(pathToFileURL(reportPath).href);
  // Only the last/both variants affect the 3′ end: 4 of 14 hits, on either strand.
  await expect(report.locator("#summary-table tbody tr")).toContainText(
    "28.6%",
  );
  await report
    .getByRole("button", { name: "Show all 14", exact: true })
    .click();
  for (const [sample, details, coordinates] of [
    ["last_minus", "40:A>C", "80–42"],
    ["insertion_minus", "18:->G", "81–41"],
  ] as const) {
    await report.getByRole("button", { name: sample, exact: true }).click();
    await expect(report.locator("#alignment-modal-meta")).toContainText(
      `Reverse strand · ${coordinates}`,
    );
    const body = report.locator("#alignment-modal-body");
    await expect(body).toContainText(details);
    await expect(body).toContainText(
      "For a reverse-strand hit, the subject is reverse-complemented.",
    );
    const aligned = await body.locator(".alignment-row span").allTextContents();
    expect(aligned[0].length).toBe(aligned[2].length);
    expect(aligned[1].split(" ").length - 1).toBe(1);
    if (sample === "insertion_minus") {
      expect(aligned[0]).toBe("ACGTTGCAAGCTTAGCGA-TCGATGCTAGCAGTTCGACCTA");
      expect(aligned[2]).toBe("ACGTTGCAAGCTTAGCGAGTCGATGCTAGCAGTTCGACCTA");
      await report.screenshot({
        path: testInfo.outputPath("reverse-insertion-report.png"),
      });
    }
    await report.getByRole("button", { name: "Close", exact: true }).click();
  }
  await report.getByRole("button", { name: "Norsk", exact: true }).click();
  await report
    .getByRole("button", { name: "insertion_minus", exact: true })
    .click();
  await expect(report.locator("#alignment-modal-meta")).toContainText(
    "Revers tråd · 81–41",
  );
  await expect(report.locator("#alignment-modal-body")).toContainText(
    "Ved treff på revers tråd er prøvesekvensen reverskomplementert.",
  );
});
