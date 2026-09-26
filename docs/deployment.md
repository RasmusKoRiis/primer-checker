# Deployment

## Current verification and remaining account steps

The feature branch is `feat/webapp-vercel`. The available GitHub account is
`RasmusKRiis`, which has read-only access to `RasmusKoRiis/primer-checker`.
The branch is therefore published to the fork `RasmusKRiis/primer-checker` for
a pull request into the original repository's `main` branch. Nothing is merged.

Vercel CLI is logged out in the implementation environment. No Vercel project,
preview deployment, production deployment, or DNS change has been made.
The production Next.js build, local production routes, and the API are tested.
GitHub Actions additionally packages and executes BLAST in Amazon Linux 2023
and runs browser tests against the production build. A **real Vercel preview
is still required** to confirm function collection, routing, bundle contents,
and execution in Vercel itself; a Linux CI pass is not that verification.

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

## Connect GitHub and create a preview

1. Sign in to the intended Vercel account and grant its GitHub integration access
   to `RasmusKoRiis/primer-checker`.
2. Import the repository as project **primer-checker**, root **.**, framework
   **Next.js**, Node **22.x**. Keep the committed install/build settings.
3. Set the production branch to **main**. Enable Fluid compute and keep the
   standard Hobby memory allocation. No project-specific secrets are required.
4. Preview the feature branch or its pull request. For this fork-based PR,
   the repository/Vercel owner must authorize the preview. Alternatively, an
   account with upstream write access can fetch the feature branch and push it
   under the same name to the original repository. Do not promote it to production.
5. Confirm the deployment logs show the BLAST bundle verification. In the preview,
   open `/api/health` and `/api/catalog`, upload `public/example.fasta`, run PCR,
   inspect the single mismatch, and download CSV/HTML/provenance. Check an NGS
   panel and a malformed upload too. Confirm no `node_modules` or `.next` directory
   is included in the Python function bundle.
6. Production is a separate review/merge decision: after approval, merge into
   `main` and let that branch deploy. This task does not authorize the merge.

CLI alternative after account connection:

```bash
npx vercel login
npx vercel link
npx vercel pull --environment=preview
npx vercel                 # preview only; do not add --prod
```

Keep production set to `main`; verify it in project settings before deploying.
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
2,000 primer/record comparisons, 1,200 BLAST subprocesses, and a 240-second
analysis deadline. Metadata is limited to 2,000 rows, 100 columns, and 500
characters per cell. Large responses are rejected with guidance to narrow the
analysis. Limits are validated before BLAST whenever their size is known.

Uploads are untrusted UTF-8 nucleotide FASTA and optional CSV. ZIP, FASTQ,
gapped alignments, duplicate record identifiers within a file, and reserved
BLAST ID prefixes are rejected. Influenza needs supported segment tokens such
as `01-HA|sample`, `03-M|sample`, or `08-NS|sample`. Duplicate sample IDs across
different files remain separate records unless metadata links them.

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
