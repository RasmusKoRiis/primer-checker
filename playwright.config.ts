import { defineConfig } from "@playwright/test";
export default defineConfig({
  testDir: "./tests/e2e",
  timeout: 60_000,
  workers: 1,
  use: { baseURL: "http://127.0.0.1:3001", trace: "retain-on-failure" },
  webServer: [
    { command: "python3 -m uvicorn api.index:app --host 127.0.0.1 --port 8000 --no-access-log", url: "http://127.0.0.1:8000/api/catalog", reuseExistingServer: !process.env.CI, timeout: 30_000 },
    { command: "npm run start -- --port 3001", url: "http://127.0.0.1:3001", reuseExistingServer: !process.env.CI, timeout: 30_000 },
  ],
});
