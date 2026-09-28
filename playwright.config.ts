import { defineConfig } from "@playwright/test";
const port = Number(process.env.PLAYWRIGHT_PORT || 3001);
export default defineConfig({
  testDir: "./tests/e2e",
  timeout: 60_000,
  workers: 1,
  use: { baseURL: `http://127.0.0.1:${port}`, trace: "retain-on-failure" },
  webServer: [
    {
      command:
        "python3 -m uvicorn api.index:app --host 127.0.0.1 --port 8000 --no-access-log",
      url: "http://127.0.0.1:8000/api/catalog",
      reuseExistingServer: !process.env.CI,
      timeout: 30_000,
    },
    {
      command: `npm run start -- --port ${port}`,
      url: `http://127.0.0.1:${port}`,
      reuseExistingServer: !process.env.CI,
      timeout: 30_000,
    },
  ],
});
