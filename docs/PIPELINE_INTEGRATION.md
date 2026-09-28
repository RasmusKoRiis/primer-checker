# Routine Nextflow integration

The routine wrappers enable PCR checks for influenza and PCR plus NGS checks
for SARS-CoV-2 and RSV. Bare pipeline invocations remain opt-in with
`--primer_check true`. PCR and NGS run as independent tasks, each with
`errorStrategy 'ignore'`. A technical failure in either task does not stop the
sequencing workflow or suppress the other assay's outputs.

## Wrapper settings

PCR locations are defined in the wrappers, not in pipeline defaults:

| Wrapper | Default PCR database directory |
| --- | --- |
| fluseq | `/mnt/tempdata/influensa_db/flu_seq_db/pcr-primers` |
| rsvseq | `/mnt/tempdata/rsv_db/pcr-primers` |
| sarsseq | `/mnt/tempdata/sars_db/pcr-primers` |

Use `-P /path/to/primers.json` or `-P /path/to/pcr-primers` to override.
A directory loads all its top-level `*.json` databases in sorted order;
duplicate organism/assay/primer combinations are rejected. Keep any referenced
assets alongside their JSON at the original relative paths. Do not put example
databases in the production database directory.

SARS reuses the sequencing scheme selected with `-p`. RSV reuses
`/mnt/tempdata/rsv_db/assets/primer_schemes`, selecting `RSVA/<scheme>` or
`RSVB/<scheme>` using the inferred subtype and the wrapper's `-p` setting.
For either wrapper, `-N /path/to/ngs-assets` overrides this location. For RSV,
the override is the root containing both subtype directories. FASTA-only SARS
runs need `-N` to enable an NGS check when no sequencing scheme is available.

Environment equivalents are `PRIMER_CHECK_PCR`, `PRIMER_CHECK_NGS_DIR`
(SARS/RSV), `PRIMER_CHECK_CONTAINER` and `PRIMER_CHECK_ENABLED` (default `true`).
`PRIMER_CHECK_ENABLED=false` disables the optional analysis. The pipeline
parameters are `--primer_check_pcr`, `--primer_check_ngs_dir` and
`--primer_check_container`.

NGS primer sequences come from a primer FASTA when available. Recognised names
are `SARS-CoV-2.primers.fasta`, `RSVA.primer.fasta`, `RSVB.primer.fasta`,
`primers.fasta` and `primer.fasta`. Otherwise an ARTIC-style BED must have actual
oligo sequences in column 7. The loader prefers `<virus>.primer.bed` over
`<virus>.scheme.bed`. It never reconstructs oligos from a reference genome:
some RSV BED intervals cover much more than the primer itself. BED coordinates
are omitted if the interval length does not match the actual oligo length.

## Latest version and first deployment

`Dockerfile.cli` packages the Python CLI, BLAST and report resources. The
`cli-container.yml` GitHub Actions workflow tests and publishes every push to
`main` as `ghcr.io/rasmuskoriis/primer-checker:latest` and an immutable commit tag.
Make the GHCR package readable by the analysis server (public visibility, or
authenticate Docker to GHCR).

Online Docker tasks specify `--pull=always`, and the primer-check tasks have
`cache false`. Thus they use the latest **published CLI image**, including on
`-resume`; they do not silently reuse an older image if the pull fails. A failed
image pull is handled by `errorStrategy 'ignore'`, like other technical errors.
The image contains its source commit in `REVISION`; reports record that commit,
the Python and BLAST versions, and SHA-256 hashes of the primer inputs.

Before the first production run, publish the primer-checker changes so the new
image exists, then deploy pipeline revisions containing the shared integration
and the updated wrappers. Merely editing a local checkout does not update the
remote pipeline fetched by a wrapper or the published container.

SARS `-o` offline mode uses the locally cached container and does not pull.
An offline run cannot guarantee the newest remote version. Pre-pull the image
while online. Other container engines should use an explicitly refreshed image;
the routine wrappers use Docker.

Local build for development:

```bash
docker build -f Dockerfile.cli -t primer-checker:dev \
  --build-arg PRIMER_CHECKER_REVISION="$(git rev-parse HEAD)" .
```

`--primer_check_source /path/to/primer-checker` is an explicit development
override that stages local Python code instead of using the code in the image.
It disables the automatic image pull; the runtime still needs Python and BLAST.

## Routing and outputs

The integration consumes declared consensus channels. Influenza runs after
header configuration and before whole-segment coverage filtering. Its subtype
file selects H1/H3/B PCR primers using the database's available subtype labels.
`H1N1`/`H3N2` can map to legacy `H1`/`H3`; `VIC`/`VICVIC` map to B/Victoria when
present, otherwise type B. Segment matching accepts fluseq's subtype suffixes
and treats M/MP as equivalent. RSV uses `meta.subtype`; unknown RSV is never
silently treated as RSV-A. SARS uses its ARTIC consensus, with the FASTA workflow
also supported. Human influenza FASTA mode is supported; avian workflows are
outside this integration.

Successful tasks publish into `<outdir>/primer_check/`:

```text
task_status.csv
pcr_primer_report.csv
pcr_primer_report.html
pcr_primer_report.status.csv
pcr_primer_report.provenance.json
ngs_primer_report.csv                 # SARS and RSV only
ngs_primer_report.html
ngs_primer_report.status.csv
ngs_primer_report.provenance.json
```

CSV rows include run and pipeline sample IDs alongside FASTA record IDs.
True BLAST no-hit results remain rows with `Hit_Status=no_hit`. Unknown subtype,
missing assay, missing consensus and missing matching influenza segment are
listed in the status CSV. Empty analyses still produce CSV headers and HTML.
Malformed databases, unavailable primer sequences, missing input files and
BLAST execution failures exit nonzero and are ignored by Nextflow; inspect the
Nextflow task logs when an assay's reports are absent.
`task_status.csv` is refreshed even when a task fails. On a resumed run, check
this file before using reports: older reports can remain in a reused output
directory after an ignored failure, and a `failed` assay has no current report.

The pipeline adapter marks a hit as `indeterminate` if its aligned consensus
contains ambiguous bases. A no-hit result on a sequence containing unknown or
ambiguous bases is also indeterminate: the expected binding site cannot be
established from that search. Identity and mismatch counts are blank for these
rows, and the HTML report excludes them from confirmed-hit statistics.
This analysis measures compatibility with the consensus, not amplicon depth.
The existing SARS mismatch/dashboard and depth analyses remain separate.

Normal uploads include the full `primer_check` directory. CSV-only uploads
include its CSV files in the existing analysis destination; HTML remains in
the run output and normal full-run upload.

## Synthetic verification and maintenance

From the primer-checker checkout, with Python, pytest, BLAST and Nextflow 24.10.2+
available:

```bash
python3 -m pytest -q tests/test_primer_checker.py tests/test_alignment_orientation.py \
  tests/test_influenza_selection.py tests/test_pipeline_primer_check.py \
  tests/test_nextflow_primer_check.py tests/test_report_indeterminate.py
```

Tests generate synthetic PCR databases, primer FASTAs/BEDs and consensuses.
They exercise SARS, RSV-A/B, H1/H3/B, M/MP, reverse primers, no hits, masked bases,
missing assays/segments, colliding staged filenames, ignored task failures and
resuming with an uncached primer-check task. Nextflow tests run locally without
containers or network access; they skip if Nextflow is not installed.

The canonical Nextflow files live in `integrations/nextflow`. Synchronise them
into the three independent pipeline repositories after changing the module:

```bash
python3 scripts/sync_nextflow_modules.py ../nf-core-sars ../nf-core-rsvseq ../nf-core-fluseq
python3 scripts/sync_nextflow_modules.py --check ../nf-core-sars ../nf-core-rsvseq ../nf-core-fluseq
```

Each pipeline receives the same standalone `tests/primer_check` harness.
