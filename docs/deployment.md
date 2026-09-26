# Deployment

## Verified Hobby preview (2026-09-26)

Vercel project **primer-checker** is linked to this workspace under
**rasmus-projects1** (Rasmus' projects), using the existing **Hobby** plan.
Project settings are Next.js, Node **22.x**, Fluid compute, and the default
Hobby function resources. Python **3.12** comes from `.python-version`.
No paid service or plan upgrade was added.

- [Project dashboard](https://vercel.com/rasmus-projects1/primer-checker)
- [Verified preview](https://primer-checker-hx0caunf0-rasmus-projects1.vercel.app)
- [Build details](https://vercel.com/rasmus-projects1/primer-checker/D52Z1yNrxzsgcTS2rrhebjk4HE7c)
- Deployed application commit: `2889368fe3ff0612fd5b70a5fb8d170edef3227f`.

The preview retains Vercel Authentication, so open it in a browser signed in
to the owning Vercel account. `npx vercel curl` also supports authenticated
checks without publishing the preview or exposing credentials.

Verified on the actual hosted function:

- `/api/health`: BLAST **2.15.0+** executes successfully; the build verified its
  **33.5 MB** Linux bundle and four shared libraries.
- `/api/catalog`: all seven bundled PCR/NGS schemes load with the expected
  database fingerprint and upload/workload limits.
- PCR on `public/example.fasta`: six comparisons, six hits, one mismatch
  (`9:T>A`), with valid CSV, HTML, and provenance downloads.
- VMIDT 2.2 NGS panel: 68 primers and 136 comparisons complete successfully.
- Custom influenza JSON upload and H1 preflight: three records and two eligible
  comparisons, confirming segment and subtype filtering.
- A 201-record batch returns **413**; malformed FASTA returns **422**.
- The frontend HTML is served successfully, with the current input guides and
  without the removed metadata upload controls. Local browser tests cover the
  interactive workflow; the in-app browser requires a separate Vercel sign-in
  to interact with the protected hosted preview.

Only the preview is assigned. There is no active production deployment or
custom domain. Vercel classified the first deployment as production despite
`--target=preview`; it briefly assigned the default `.vercel.app` domains. That
specific deployment was removed after creating the correctly classified preview.
The production URL now returns **404**. This first-deployment behavior is also
tracked in [Vercel issue #17069](https://github.com/vercel/vercel/issues/17069).
Always inspect the returned target; the CLI flag alone is insufficient for a
brand-new project. No PR merge or DNS changes were made.

Git-based automatic deployments are not connected yet. The Vercel-linked
GitHub account lacks write/admin access to `RasmusKoRiis/primer-checker`, so
Vercel rejected that connection. The feature branch `feat/webapp-vercel` is
published to `RasmusKRiis/primer-checker`, with a draft PR into the original
repository's `main`. Choose that fork for automatic previews, or give the
Vercel-connected GitHub account write/admin access to the original repository
and allow the Vercel GitHub App to access it. Once
connected, keep the production branch set to **main**.

## One project, two runtimes

Use the **repository root** as the Vercel project root, Next.js as its framework,
Node **22.x**, and Python **3.12** (`.python-version`). The Next.js page is
prerendered; analysis runs in the Python ASGI function `api/index.py`.
`vercel.json` selects Next.js explicitly so FastAPI auto-detection cannot take
over the frontend, rewrites `/api/*` to the Python entrypoint, sets its duration,
and restricts its bundle. There are no Next.js API proxy functions in production.

The install command is `python3 scripts/package_blast.py && npm ci`; the build
command is `npm run build`. BLAST is prepared during **installation**, before
Python function file collection, rather than downloaded on incoming requests.
The function explicitly includes `bin/`, both shared engine/report modules,
`web_service/`, `primer_db/` (including NGS assets), and `report_text/`.
Virtual environments, Node dependencies, frontend build files, tests, legacy
archives, and local output are excluded from the Python bundle.

Vercel still supports file-based Python ASGI functions alongside a frontend.
Its newer Services architecture is an alternative, but is unnecessary for this
small application. See [Python functions in /api](https://vercel.com/docs/functions/runtimes/python/api-directory)
and the [Python runtime](https://vercel.com/docs/functions/runtimes/python).

## Connect GitHub and update the preview

The project already exists. Do not import a duplicate project or start a paid
trial. To enable automatic deployments:

1. Grant the Vercel-connected GitHub account write/admin access to
   `RasmusKoRiis/primer-checker` and allow the Vercel GitHub App to access it.
   Alternatively, connect `RasmusKRiis/primer-checker` after choosing the fork.
2. In the existing project's Git settings, connect the chosen repository and
   verify **main** is the production branch. Keep the repository root and
   committed Next.js build settings.
3. Use `feat/webapp-vercel` for previews. A fork-based PR into the original
   repository can require the repository/Vercel owner to authorize deployment.
4. Verify `/api/health`, the catalog, an analysis, and downloads after updates.
   Inspect the deployment target and assigned aliases before sharing its URL.
5. Production remains a separate approval/merge decision. Do not promote the
   development branch or merge the PR without approval.

For direct updates from this already linked workspace:

```bash
npx vercel whoami
npx vercel deploy --target=preview --scope rasmus-projects1 --yes
# Inspect the returned URL; its target must be preview (null in the REST API).
npx vercel inspect <preview-url> --scope rasmus-projects1
npx vercel curl /api/health --deployment <preview-url> --scope rasmus-projects1
```

For a new checkout, sign in and link to the existing project:

```bash
npx vercel login
npx vercel link --yes --scope rasmus-projects1 --project primer-checker
```

Local `.vercel/` connection data and `.env.local` are ignored and must not be
committed. Vercel may create a local OIDC token when linking; the application
itself needs no API key. The Hobby plan is for
[personal, non-commercial projects](https://vercel.com/docs/plans/hobby), with
usage caps. Existing per-analysis limits do not enforce the account's monthly
allowance. Review team usage before opening the service widely.

Vercel Hobby checks commit authorship/account ownership, and fork pull requests
require deployment authorization. See [Vercel Git integration](https://vercel.com/docs/git)
and [GitHub integration](https://vercel.com/docs/git/vercel-for-github).

## BLAST bundle

The pinned release is **NCBI BLAST+ 2.15.0**, matching the discovery baseline.
It uses the official `ncbi-blast-2.15.0+-x64-linux.tar.gz` archive. The archive is
about 246 MiB to download; only the 32,743,960-byte `blastn` executable and its
non-glibc shared dependencies are deployed. No reference BLAST database or
`makeblastdb` executable is needed: the engine uses `-subject`.

The verified Amazon Linux bundle totals **33,483,200 bytes (33.5 MB)**. Its four
libraries are `libbz2.so.1`, `libgcc_s.so.1`, `libgomp.so.1`, and `libz.so.1`.
The CI artifact records each hash. Vercel's current build-image libraries may
have different patch versions; the script verifies their executability together.

`scripts/package_blast.py` verifies the archive and extracted binary with pinned
SHA-256 digests. The original archive's published NCBI MD5 was also checked
during implementation. It reads only the expected regular-file tar member,
sets mode 0755, resolves dependencies using `ldd`, copies non-glibc libraries to
`bin/lib/`, runs `blastn -version`, and performs a small synthetic alignment.
It fails on missing dependencies, checksum drift, an unexpected version, or a
bundle over 60 MB. `bin/manifest.json` records the final library hashes and sizes.

The application's subprocess environment adds the adjacent `lib/` directory to
`LD_LIBRARY_PATH` per invocation. glibc and its loader are deliberately supplied
by the execution environment. Vercel's current build image is Amazon Linux 2023;
CI tests the same distribution. If its library inventory changes, fix the
build-image/library mismatch before deployment rather than suppressing `ldd`
errors. [Build image documentation](https://vercel.com/docs/builds/build-image).

Executable resolution:

1. `BLASTN_PATH`, when set (an invalid override is an error, not a fallback).
2. Repository-relative `bin/blastn`, if present.
3. `blastn` from PATH.

On macOS use your native BLAST installation; do not copy the Linux binary into
the local `bin/` directory. To verify the deployment bundle with Docker:

```bash
docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work amazonlinux:2023 \
  /bin/sh -c 'dnf install -y python3 libgomp bzip2-libs && python3 scripts/package_blast.py'
```

This command intentionally creates a **Linux-only** local bundle. Move it aside
or set `BLASTN_PATH` to your native executable before running the CLI on macOS.
Do not commit binaries or generated reports. For an offline repeat build,
provide `--archive /path/to/the/pinned/archive.tar.gz`; checksums still apply.

## Environment variables

| Variable | Default / use |
| --- | --- |
| `BLASTN_PATH` | Optional explicit executable; normally leave unset on Vercel. |
| `PRIMER_DATABASE_PATH` | Optional operator-controlled database path; defaults to the included unified database. Its assets must also be bundled. Never accepted from request input. |
| `APP_GIT_COMMIT` | Optional local reproducibility field. |
| `VERCEL_GIT_COMMIT_SHA` | Automatically supplied by Vercel and preferred over `APP_GIT_COMMIT`. |
| `PYTHON` | Optional local interpreter override for `npm run dev`. |

No API keys, persistent database, blob store, or other paid service is required.
FastAPI reads process environment variables; it does not automatically load
Next.js `.env.local` files. Export Python settings in the shell or set them in
the deployment dashboard. `.env.example` documents the available overrides.

## Limits and privacy

As checked on 2026-09-26, Vercel Functions allow **4.5 MB request and response
payloads**. Hobby with Fluid compute allows **300 seconds**, **2 GB / 1 vCPU**,
and a standard **500 MB uncompressed Python bundle**. These are platform limits,
not an assurance that every input within them finishes. See
[Vercel function limits](https://vercel.com/docs/functions/limitations).

This application uses lower limits: 3,000,000 combined upload bytes, 4,000,000
wire/request and serialized-response bytes, 10 FASTA files, 200 sequence records,
2,000 primer/record comparisons, 300 BLAST subprocesses, 50 million total
sequence bases × selected primers, and a 240-second analysis deadline.
The combined upload budget includes custom database JSON. A custom database is
separately capped at 250,000 bytes, 500 primers, and 200 bases per primer. Metadata is limited to 2,000 rows, 100 columns, and
500
characters per cell. Large responses are rejected with guidance to narrow the
analysis. The browser checks file sizes before sending data. A debounced
`/api/preflight` request then uses the same validation and primer selection as
analysis, counting the batch without running BLAST. Changing inputs invalidates
the check and disables Analyze until a new check passes. `/api/analyze` repeats
all limits so preflight cannot be bypassed by calling the API directly.

These are per-request bounds, not global rate limiting or a monthly compute
budget. Repeated small analyses can still use the Hobby allowance. Check Vercel
project usage before sharing widely; no distributed counter or paid service is
introduced here. Preflight validates and discards its uploads; the browser sends
them again when the user starts analysis.

Uploads are untrusted UTF-8 nucleotide FASTA, optional CSV, and optional
self-contained JSON primer databases. JSON uploads never resolve asset paths or
modify the installed database; supported formats are documented in the README.
ZIP, FASTQ,
gapped alignments, duplicate record identifiers within a file, and reserved
BLAST ID prefixes are rejected. Influenza needs supported segment tokens such
as `01-HA|sample`, `03-M|sample`, or `08-NS|sample`. Duplicate sample IDs across
different files remain separate records in the website. Direct API callers may
link records with sample metadata.

Input and query files use temporary directories/files, with cleanup on success
and failures. Requests are not stored in a database or application logs. HTTP
errors contain safe messages rather than tracebacks or filesystem paths.
Results stay in the browser's memory until navigation/reload; downloaded files
are deliberately retained on the user's computer. Infrastructure providers
still process requests and may keep operational metadata under their own
policies; this is not a confidential-data deployment.

Use public, synthetic, or anonymized sequence data. Do not upload confidential,
identifiable, restricted, or identifiable surveillance/NIPH data. The interface
contains the same notice. No suitability for confidential workloads is claimed.

## Performance and smallest fallback

The baseline three-file synthetic CLI batch took 2.033s. A later run of the
68-primer VMIDT 2.2 panel against two short synthetic records took 56.243s on the
implementation Mac (136 rows, roughly 413 KB JSON before provenance columns
were added). These are fixture measurements, not cloud throughput benchmarks.
The canonical engine still runs one BLAST per primer/file. In particular,
all-panel NGS analysis can be slow. No search parameters or mismatch definitions
were changed to improve apparent speed.

If a Vercel preview cannot package or execute BLAST cleanly, keep the frontend
on Vercel and run **this same FastAPI application** on a small Linux x86_64
container service such as Cloud Run. Install the pinned bundle at image build
time, launch `uvicorn api.index:app --host 0.0.0.0 --port 8080`, and proxy `/api/*`
to that service. Do not deploy or provision a paid service without a separate
decision. For large surveillance batches, the existing CLI is already usable.

For a future larger web workload, direct-to-object-storage uploads with short
expiry and asynchronous jobs would remove the request-body bottleneck. That
would require explicit storage retention, authentication, and job lifecycle
design; none is introduced in this release.

## Future custom domain

After `rasmuskriis.no` is registered and production is approved, add
`primers.rasmuskriis.no` in Vercel's project Domains settings. Follow the exact
DNS verification/CNAME records Vercel provides for that project, then confirm
TLS issuance. Do not assume a fixed CNAME target. No DNS records are changed
by this implementation.
