import { spawn } from "node:child_process";
import { existsSync } from "node:fs";
const python = process.env.PYTHON || (existsSync(".venv/bin/python") ? ".venv/bin/python" : "python3");
const processes = [
  spawn(python, ["-m", "uvicorn", "api.index:app", "--host", "127.0.0.1", "--port", "8000", "--reload", "--reload-dir", "api", "--reload-dir", "web_service", "--no-access-log"], { stdio: "inherit" }),
  spawn(process.execPath, ["node_modules/next/dist/bin/next", "dev"], { stdio: "inherit" }),
];
let stopping = false;
function stop(code = 0) {
  if (stopping) return;
  stopping = true;
  for (const child of processes) child.kill("SIGTERM");
  process.exitCode = code;
}
for (const child of processes) {
  child.on("error", (error) => { console.error(error.message); stop(1); });
  child.on("exit", (code) => stop(code || 0));
}
process.on("SIGINT", () => stop());
process.on("SIGTERM", () => stop());
