# Simple Primer Checker Test Data

This folder contains a minimal mixed-virus batch dataset for quick local checks.
It has one FASTA file each for SARS-CoV-2, RSV-A, and Influenza H3, plus a small legacy primer JSON and metadata CSV.

Run analysis plus HTML report:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder fixtures/simple_test_data \
  --primers fixtures/simple_test_data/simple_primers.json \
  --metadata-csv fixtures/simple_test_data/metadata.csv \
  --output result/simple_test_report.csv \
  --html-report result/simple_test_report.html
```

Run analysis only:

```bash
python3 scripts/run_primer_checker_batch.py \
  --input-folder fixtures/simple_test_data \
  --primers fixtures/simple_test_data/simple_primers.json \
  --metadata-csv fixtures/simple_test_data/metadata.csv \
  --output result/simple_test_report.csv \
  --analysis-only
```
