# Primer Database Maintenance Proposal

## Current problem

The current primer database is easy to read but too implicit:

```json
{
  "Influenza-B": {
    "triplex_InfB_F_NS": "TCCTCAAYTCACTCTTCGAGCG"
  }
}
```

This stores only primer name and sequence. Important metadata such as segment, role, organism, scheme version, source, and references must be inferred from names. That makes updates fragile.

## Recommended database shape

Keep supporting the legacy file, but add a normalized versioned format for future updates:

```json
{
  "schema_version": "1.0",
  "database_version": "2026-05-06",
  "schemes": [
    {
      "scheme_id": "fhi-influenza-b-triplex",
      "display_name": "FHI Influenza B triplex",
      "organism": "Influenza-B",
      "version": "2026-05-06",
      "status": "current",
      "source": "FHI",
      "references": [],
      "primers": [
        {
          "id": "triplex_InfB_F_NS",
          "name": "triplex_InfB_F_NS",
          "sequence": "TCCTCAAYTCACTCTTCGAGCG",
          "role": "forward_primer",
          "segment": "NS",
          "gene": "NS",
          "pool": "triplex",
          "strand": "plus",
          "notes": ""
        }
      ]
    }
  ]
}
```

This format is now supported by `primer_checker.py`. The legacy FHI primer file was converted into:

```text
primer_db/fhi_primers.normalized.json
```

## Easy update workflow

1. Edit one scheme file or an input spreadsheet/CSV, not the analysis code.
2. If starting from the old legacy JSON, run:

```bash
python3 scripts/convert_legacy_primers.py \
  --input old_primers.json \
  --output primer_db/fhi_primers.normalized.json \
  --database-version 2026-05-06 \
  --source FHI
```

3. Run a converter if using CSV in the future:

```bash
python3 scripts/convert_primers_csv.py --input new_scheme.csv --output primer_db/schemes/influenza-b/fhi-triplex/2026-05-06.json
```

4. Validate before use:

```bash
python3 primer_checker.py --primers primer_db/schemes/influenza-b/fhi-triplex/2026-05-06.json --validate-primers
```

5. Run analysis using that exact versioned file.
6. Keep old scheme versions in the repository for reproducibility.

## Why this is easier to maintain

- Adding a primer no longer depends on naming conventions.
- Segment and role are explicit, so NS/HA/M filtering is safer.
- Old versions can be retained and rerun.
- Validation can catch invalid IUPAC characters, missing metadata, duplicate IDs, and unsupported segments before analysis.
- Reports can eventually include `database_version`, `scheme_id`, and `scheme_version`.

## Suggested implementation steps

Completed:

- Added support for the normalized JSON shape in the existing loader.
- Kept legacy JSON support as compatibility mode.
- Added pure-Python validation for required normalized fields.
- Added a converter from the current legacy format to normalized JSON.
- Converted the current FHI primer database into `primer_db/fhi_primers.normalized.json`.

Remaining:

- Optionally add CSV-to-JSON conversion for easier non-programmer updates.
- Add a JSON Schema file if external validation tooling becomes useful.
- Add report columns for `database_version`, `scheme_id`, and `scheme_version` if downstream consumers want them.
