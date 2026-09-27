# DUMMY primer data — software testing only

All sequences in this directory are invented examples, not a real or validated
primer library. Organism, subtype and segment labels demonstrate software
selection rules only. Upload or build your own database for your analysis.

- `dummy_primers.json`: self-contained website default, version `dummy-1.0`.
  Three PCR oligos and two NGS oligos for `Demo-virus`, plus eight synthetic
  Influenza-A segment examples. Downloadable and re-uploadable on the website.
- `../public/example.fasta`: two artificial templates for the Demo-virus PCR/NGS
  examples. The second record has one forward-primer mismatch at position 9.
  The reverse-primer binding site is the reverse complement of the input oligo.
- `dummy-influenza.fasta`: matching artificial records for all eight segments.
- `assets/panel_database.example.json`: schema 3 CLI example using the synthetic
  BED and oligo FASTA in `assets/dummy/`.

The `purpose: "synthetic-test-only"` marker labels web uploads and result exports
as dummy data. Keep it when copying test data. For your own real library, replace
all examples, use your own metadata, and omit the synthetic purpose marker.

The former real database and assets are no longer shipped in this branch's
current files. This does not erase earlier Git commits or other copies.
