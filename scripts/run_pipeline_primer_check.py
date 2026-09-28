#!/usr/bin/env python3
"""Run one assay type for an explicitly routed Nextflow consensus manifest.

The manifest is a list of {sample_id, fasta: [paths], subtype, subtype_file}.
No sample/virus identity is inferred from filenames. PCR input is a JSON file
or a directory of JSON databases, with external assets kept beside each JSON.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
import primer_analysis as analysis  # noqa: E402
import primer_report as report  # noqa: E402


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_pcr_databases(path):
    # Nextflow may stage a single JSON as a symlink with a different basename.
    # Resolve it before loading relative panel assets from the same database.
    root = Path(path).resolve()
    files = sorted(root.glob('*.json')) if root.is_dir() else [root]
    if not files:
        raise ValueError(f'No PCR JSON databases found in {root}')
    combined = {}
    seen = set()
    provenance = []
    for database in files:
        records, validation = analysis.load_primer_records(str(database))
        analysis.print_validation_messages(validation)
        provenance.append({'file': str(database), 'sha256': file_digest(database)})
        for organism, primers in records.items():
            for primer in primers:
                if primer.assay_type != 'pcr':
                    continue
                key = (organism, primer.scheme_id, primer.name)
                if key in seen:
                    raise ValueError(f'Duplicate PCR primer across databases: {key}')
                seen.add(key)
                combined.setdefault(organism, []).append(primer)
    if not combined:
        raise ValueError(f'No PCR primers found in {root}')
    return combined, provenance


def first_file(directory, names):
    return next((directory / name for name in names if (directory / name).is_file()), None)


def load_ngs_panel(root, virus, scheme):
    directory = Path(root)
    short = virus.replace('-', '') if virus.startswith('RSV-') else virus
    if virus.startswith('RSV-'):
        if not scheme:
            raise ValueError('RSV NGS checks require the selected primer scheme')
        directory = directory / short / scheme
    if not directory.is_dir():
        raise FileNotFoundError(f'NGS primer directory does not exist: {directory}')
    fasta = first_file(directory, [f'{virus}.primers.fasta', f'{short}.primer.fasta',
                                   'primers.fasta', 'primer.fasta'])
    bed = first_file(directory, [f'{short}.primer.bed', 'primer.bed',
                                 f'{short}.scheme.bed', 'ncov-2019_midnight.scheme.bed'])
    if fasta is None and bed is None:
        raise ValueError(f'No primer FASTA or BED in {directory}')
    version = scheme or directory.name
    panel_id = f'{virus}:{version}'
    records = []
    # Prefer actual oligo sequences. Never manufacture primers from reference
    # coordinates: some RSV BED files describe much larger amplicon intervals.
    sequences = analysis.read_panel_fasta_records(str(fasta)) if fasta else {}
    entries = analysis.read_panel_bed_records(str(bed)) if bed else []
    metadata = {entry['name']: entry for entry in entries}
    if not fasta:
        if not entries or any(not entry['sequence'] for entry in entries):
            raise ValueError(f'{bed}: BED needs oligo sequences in column 7, or a primer FASTA')
        sequences = {entry['name']: entry['sequence'] for entry in entries}
        if len(sequences) != len(entries):
            raise ValueError(f'Duplicate primer names in {bed}')
    if not sequences:
        raise ValueError(f'No NGS primer sequences in {directory}')
    for name, sequence in sequences.items():
        if not sequence or set(sequence.upper()) - analysis.ALLOWED_SEQUENCE_CODES:
            raise ValueError(f'Invalid primer sequence: {name}')
        entry = metadata.get(name, {})
        # A matching name can still describe an amplicon rather than an oligo.
        # Retain pool/strand, but do not report those intervals as binding sites.
        start, end = entry.get('start', ''), entry.get('end', '')
        if start and end and int(end) - int(start) != len(sequence):
            start, end = '', ''
        records.append(analysis.PrimerRecord(
            organism=virus, name=name, sequence=sequence.upper(), assay_type='ngs',
            scheme_id=panel_id, scheme_version=version, assay_name=panel_id,
            pool=entry.get('pool', ''), strand=entry.get('strand', ''),
            start=start, end=end,
            reference_name=entry.get('chrom', ''),
        ))
    assets = [{'file': str(path), 'sha256': file_digest(path)} for path in (fasta, bed) if path]
    return records, assets


def select_target(virus, subtype, records):
    subtype = subtype.strip().upper()
    if virus == 'SARS-CoV-2':
        return virus, None
    if virus == 'RSV':
        return {'RSVA': ('RSV-A', None), 'RSV-A': ('RSV-A', None),
                'RSVB': ('RSV-B', None), 'RSV-B': ('RSV-B', None)}.get(subtype)
    choices = analysis.influenza_selections(records)
    if subtype in {'VIC', 'VICVIC', 'VICTORIA', 'B/VICTORIA'}:
        for choice in ('B/VICTORIA', 'B/VIC', 'B'):
            if choice in choices:
                return 'influenza', choice
        return None
    if subtype in choices:
        return 'influenza', subtype
    # Preserve the existing H1N1 -> H1 and H3N2 -> H3 database convention.
    family = re.fullmatch(r'(H\d+)N\d+', subtype)
    if family and family[1] in choices:
        return 'influenza', family[1]
    # An A database with only shared (untagged) primers is also valid.
    if re.fullmatch(r'H\d+(?:N\d+)?', subtype) and 'A' in choices:
        a_primers = analysis.build_influenza_subtype_records(records, 'A')
        if a_primers and all(not primer.subtype_tags for primer in a_primers):
            return 'influenza', 'A'
    return None


def run(args):
    samples = json.loads(Path(args.manifest).read_text())
    if not isinstance(samples, list):
        raise ValueError('The consensus manifest must be a list')
    records, assets = ({}, [])
    if args.assay_type == 'pcr':
        if not args.pcr_db:
            raise ValueError('PCR database is missing; set --primer_check_pcr in the pipeline')
        records, assets = load_pcr_databases(args.pcr_db)
    elif not args.ngs_dir:
        raise ValueError('NGS primer assets are missing; set --primer_check_ngs_dir')
    analysis.ensure_blastn_available()
    rows, statuses, panel_cache = [], [], {}
    used_assays = {}
    seen_samples = set()
    for sample in samples:
        sample_id = str(sample['sample_id'])
        if sample_id in seen_samples:
            raise ValueError(f'Duplicate sample in consensus manifest: {sample_id}')
        seen_samples.add(sample_id)
        subtype = sample.get('subtype') or ''
        if sample.get('subtype_file'):
            subtype = Path(sample['subtype_file']).read_text().strip()
        target = select_target(args.virus, subtype, records)
        base_status = {'Run_ID': args.run_id, 'Sample_ID': sample_id, 'Assay_Type': args.assay_type}
        if target is None:
            statuses.append({**base_status, 'Status': 'unclassified', 'Detail': subtype or 'No subtype'})
            continue
        virus, flu_type = target
        if args.assay_type == 'ngs':
            if virus not in panel_cache:
                panel_cache[virus], panel_assets = load_ngs_panel(args.ngs_dir, virus, args.ngs_scheme)
                assets.extend(panel_assets)
            selected = panel_cache[virus]
        else:
            try:
                virus, selected = analysis.select_primer_records(records, virus, flu_type, assay_type='pcr')
            except SystemExit as exc:
                statuses.append({**base_status, 'Status': 'no_assay', 'Detail': str(exc)})
                continue
        for primer in selected:
            used_assays[(primer.organism, primer.scheme_id)] = {
                'organism': primer.organism, 'assay_id': primer.scheme_id,
                'assay_version': primer.scheme_version, 'database_version': primer.database_version,
            }
        sample_rows = []
        files = sample.get('fasta', [])
        if not files:
            statuses.append({**base_status, 'Status': 'no_consensus', 'Detail': 'No consensus files'})
            continue
        for fasta in files:
            if not Path(fasta).is_file():
                raise FileNotFoundError(f'Consensus file missing: {fasta}')
            file_rows = analysis.process_fasta_file(
                str(fasta), virus, selected,
                execution=analysis.BlastExecution(strict_errors=True),
            )
            sequences = analysis.read_fasta_sequences(str(fasta))
            for row in file_rows:
                subject = sequences[row['Subject_Sequence_ID']]
                aligned_ambiguity = row['Hit_Status'] == 'hit' and any(
                    base not in 'ACGT-' for base in row['Subject_Alignment'].upper())
                incomplete_no_hit = row['Hit_Status'] == 'no_hit' and any(
                    base not in 'ACGT' for base in subject.upper())
                if aligned_ambiguity or incomplete_no_hit:
                    row.update(Hit_Status='indeterminate', Percent_Identity='', Mismatches='')
            sample_rows.extend(file_rows)
        present = {(row['Assay_ID'], row['Primer_Name']) for row in sample_rows}
        missing = [p.name for p in selected if (p.scheme_id, p.name) not in present]
        if missing:
            statuses.append({**base_status, 'Status': 'no_matching_segment', 'Detail': ', '.join(missing)})
        for row in sample_rows:
            row.update(Run_ID=args.run_id, Sample_ID=sample_id)
        rows.extend(sample_rows)
        if sample_rows:
            statuses.append({**base_status, 'Status': 'analysed', 'Detail': f'{len(sample_rows)} comparisons'})
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with Path(f'{prefix}.csv').open('w', newline='') as handle:
        analysis.write_csv_rows(rows, handle, extra_fieldnames=('Run_ID', 'Sample_ID'))
    # A legitimate empty analysis must still satisfy the workflow's outputs.
    report.write_html_report(rows, f'{prefix}.html')
    with Path(f'{prefix}.status.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['Run_ID', 'Sample_ID', 'Assay_Type', 'Status', 'Detail'])
        writer.writeheader()
        writer.writerows(statuses)
    revision = REPO_ROOT / 'REVISION'
    provenance = {
        'run_id': args.run_id, 'assay_type': args.assay_type,
        'primer_checker_revision': revision.read_text().strip() if revision.exists() else 'local-source',
        'python': sys.version.split()[0],
        'blast': subprocess.check_output([analysis.resolve_blastn(), '-version'], text=True).splitlines()[0],
        'primer_assets': assets, 'manifest_sha256': file_digest(args.manifest),
        'assays': list(used_assays.values()),
        'comparisons': len(rows),
    }
    Path(f'{prefix}.provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--virus', required=True, choices=['SARS-CoV-2', 'RSV', 'Influenza'])
    parser.add_argument('--assay-type', required=True, choices=['pcr', 'ngs'])
    parser.add_argument('--pcr-db')
    parser.add_argument('--ngs-dir')
    parser.add_argument('--ngs-scheme', default='')
    parser.add_argument('--run-id', default='Unknown')
    parser.add_argument('--output-prefix', required=True)
    args = parser.parse_args()
    try:
        run(args)
    except (OSError, ValueError, analysis.BlastError) as exc:
        parser.exit(1, f'Primer check failed: {exc}\n')


if __name__ == '__main__':
    main()
