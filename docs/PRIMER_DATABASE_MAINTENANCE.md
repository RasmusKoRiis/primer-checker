# Maintaining a primer database

The repository ships **synthetic dummy data only**. The website and batch runner
use `primer_db/dummy_primers.json` by default. The invented sequences demonstrate
PCR, NGS and influenza routing; they are not a validated assay library. No remote
primer database is downloaded at runtime.

## Bring your own primers

Use **Upload a database** or **Build a database** on the website. Enter each actual
oligo in its 5′ → 3′ direction, including reverse primers. Download the resulting
JSON for reuse. Uploads are temporary and do not replace the installed database.

The recommended self-contained format is schema `1.0` with `schemes[]`, a
`database_version`, and a version for each scheme. Each primer has an `id`,
`name`, `sequence`, `role`, and `segment`. PCR and NGS lists use the optional
`assay_type` field. See the dummy JSON for a complete working example, or the
website's database format preview for a smaller template.

For influenza, use an organism such as `Influenza-A`, explicit segment labels,
and `subtype_tags`. Untagged primers apply to every subtype of that organism;
tagged primers apply to the exact selected tag. Header segment labels must match
the database. The synthetic `primer_db/dummy-influenza.fasta` illustrates all
eight standard Influenza-A segment labels; additional labels and subtypes work
when supplied by your database.

## Keep dummy and real data distinct

The supplied JSON includes `purpose: "synthetic-test-only"`, a `dummy-1.0`
version, and `DUMMY` in every assay and primer name. Keep those labels on test
copies. The website preserves the purpose marker when this file is uploaded,
and labels its result tables, HTML, CSV provenance and manifest accordingly.

When creating a real library, replace all synthetic sequences and their labels,
use your own versions and source information, and omit the synthetic purpose
marker. Omission does not imply scientific validation. Validate your library
and independently verify its suitability for your work.

The old real library, panel assets, archive and generated report have been
removed from the current branch. Existing Git history, other branches, old
releases and old deployments are separate copies; removing current files does
not remove those copies.

## CLI workflow

Convert a legacy organism → primer → sequence dictionary when needed:

```bash
python3 scripts/convert_legacy_primers.py \
  --input /path/to/legacy-primers.json \
  --output /path/to/my-primers.json \
  --database-version 1.0 \
  --source "Your laboratory"
python3 primer_checker.py --primers /path/to/my-primers.json --validate-primers
```

Record the exact database version used for an analysis and keep a copy in your
own storage. The CLI also supports schema `3.0` and BED/FASTA-backed panels;
`primer_db/assets/panel_database.example.json` is a small synthetic example.
Asset paths resolve relative to the database JSON. Public website uploads must
be self-contained.

To run the bundled demonstration:

```bash
python3 primer_checker.py \
  --primers primer_db/dummy_primers.json --virus Demo-virus --assay-type pcr \
  --fasta public/example.fasta --output result/dummy.csv --html-report result/dummy.html
```

Expected: two records, three primers, six hits, one `9:T>A` mismatch, and reverse
primer hits with decreasing subject coordinates. Assay and primer names retain
`DUMMY` in CLI CSV and HTML output, and the CLI prints a dummy-data notice.
