# Web application implementation note

Discovery baseline: `main` at `33af65b`, 2026-09-26. No repository-specific
AGENTS.md instructions were present. The workspace was empty and was cloned
from the existing GitHub repository before inspection.

## Existing behavior

- `primer_checker.py` is the CLI and compatibility import surface.
- `primer_analysis.py` owns database loading/validation, PCR/NGS selection,
  influenza segment/subtype filtering, metadata matching, alignment
  reconstruction, IUPAC mismatch calculation, and CSV output.
- `primer_report.py` renders the existing standalone HTML report, including
  its existing risk interpretation and bilingual report text.
- The batch script infers virus/subtype from filenames and calls the same engine.
- Runtime Python code currently uses only the standard library. BLAST+ is an
  external executable; pytest is the only existing test dependency.
- BLAST runs once per primer per FASTA, using an argument array, a temporary
  query file, `-subject`, reward 2, penalty -3, word size 4, dust yes, and
  max_target_seqs 5000. Alignment selection and mismatch definitions must stay
  unchanged. A failed process currently prints stderr and returns no hits.

## Baseline verification

All 32 existing tests pass on Python 3.10.5. The three-file synthetic batch
fixture, using its bundled legacy database and local BLAST 2.15.0+, took 2.033s
including process startup and CSV/HTML generation: 3 rows, all hits, 910-byte
CSV and 102,493-byte HTML. Baseline output is outside the repository in a
temporary directory. No algorithmic optimization is justified by this tiny
measurement; large NGS panels will need a separate measurement.

## Implementation approach

- Keep the CLI and engine in place. Add optional execution controls to the
  engine for strict BLAST errors, a shared deadline, and quiet web operation;
  CLI defaults retain existing behavior.
- Add a Python `web_service/` adapter and FastAPI entrypoint in `api/index.py`.
  Validate bounded multipart uploads, call existing engine functions, and
  return structured rows plus CSV/HTML downloads in one bounded JSON response.
  Each request uses its own temporary directory, with no persistent jobs.
- Add Next.js/TypeScript in `app/` with native configuration, tables, filters,
  alignment inspection, and browser-side downloads. No analytical calculations
  move into JavaScript; table aggregation is descriptive only.
- Use one repository-root Vercel project: Next.js frontend plus a file-based
  Python ASGI function. Explicitly select Next.js to avoid Python framework
  auto-detection taking over the frontend, and rewrite API paths to that function.
- Package only the pinned official Linux x86_64 `blastn` and any required
  libraries, verify its checksum, executable permissions, and Linux dependencies
  during build. Resolve BLASTN_PATH, then the bundle, then PATH.

## Deployment risks and limits

Vercel's current documentation lists a 4.5 MB request **and response** limit,
300s maximum duration on Hobby with Fluid compute, 2 GB/1 vCPU, and a 500 MB
uncompressed Python bundle limit. Enforce lower application limits, including
result size and comparison/work limits, leaving headroom for multipart/JSON.
Do not simulate backend progress. Long-running surveillance analysis stays CLI.

Python packaging/build ordering and native BLAST compatibility require Linux
build verification; macOS success alone is insufficient. If native execution
cannot be verified on Vercel, use the same FastAPI app on a small container
service and retain the frontend on Vercel. Do not introduce paid infrastructure
or DNS changes. Vercel credentials have not been found locally; GitHub access
is available. A preview may require the account owner to connect Vercel.

Sources checked 2026-09-26:
- https://vercel.com/docs/functions/limitations
- https://vercel.com/docs/functions/runtimes/python
- https://vercel.com/docs/functions/runtimes/python/api-directory
- https://nextjs.org/docs/app/getting-started/installation
