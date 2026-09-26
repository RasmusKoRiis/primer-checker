import { expect, test } from "@playwright/test";
import path from "node:path";

test("upload, analyze with real BLAST, inspect mismatch, and download reports", async ({ page }) => {
  await page.goto("/");
  await expect(page.getByLabel("Virus", { exact: true })).toBeEnabled();
  await page.getByLabel("Upload FASTA files").setInputFiles(path.resolve("public/example.fasta"));
  await page.getByRole("button", { name: "Analyze sequences", exact: true }).click();
  await expect(page.getByRole("heading", { name: "Compatibility results" })).toBeVisible();
  await expect(page.getByText("6 comparisons", { exact: true })).toBeVisible();
  await page.getByRole("tab", { name: "By sample 6", exact: true }).click();
  await page.getByRole("combobox", { name: "Mismatches", exact: true }).selectOption("any");
  await expect(page.getByRole("button", { name: "Inspect", exact: true })).toHaveCount(1);
  await page.getByRole("button", { name: "Inspect", exact: true }).click();
  await expect(page.getByRole("dialog")).toBeVisible();
  await expect(page.getByRole("dialog").getByText("9:T>A", { exact: true })).toBeVisible();
  await expect(page.locator(".alignment-column.mismatch")).toHaveCount(1);
  await page.getByRole("button", { name: "Close alignment" }).click();
  for (const [label, filename] of [["CSV", "primer-results.csv"], ["HTML report", "primer-report.html"]]) {
    const downloaded = page.waitForEvent("download");
    await page.getByRole("button", { name: label, exact: true }).click();
    expect((await downloaded).suggestedFilename()).toBe(filename);
  }
});

test("mobile upload and influenza configuration remain usable", async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("/");
  await page.getByRole("combobox", { name: "Virus", exact: true }).selectOption("influenza");
  await page.getByLabel("Influenza subtype").selectOption("H1");
  await expect(page.getByLabel("Influenza subtype")).toHaveValue("H1");
  await expect(page.getByRole("button", { name: "Choose files", exact: true })).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
});
