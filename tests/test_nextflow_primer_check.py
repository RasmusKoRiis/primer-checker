"""Run the actual shared Nextflow subworkflow on synthetic consensus data."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from test_pipeline_primer_check import ROOT, OLIGO, csv_rows, write_database, write_fasta, write_mixed_database

pytestmark = pytest.mark.skipif(not shutil.which('nextflow'), reason='Nextflow is not installed')
HARNESS = ROOT / 'integrations/nextflow/tests/primer_check'


def execute(tmp_path, virus, manifest, database, ngs, resume=False, final_fasta=None):
    path = tmp_path / 'manifest.json'
    path.write_text(json.dumps(manifest))
    if final_fasta is None:
        final_fasta = tmp_path / 'final.fasta'
        final_fasta.write_text(''.join(Path(fasta).read_text() for sample in manifest for fasta in sample['fasta']))
    command = ['nextflow', '-C', str(HARNESS / 'nextflow.config'), 'run', str(HARNESS / 'main.nf'),
               '--test_manifest', str(path), '--test_fasta', str(final_fasta), '--test_virus', virus,
               '--primer_check_pcr', str(database), '--primer_check_ngs_dir', str(ngs),
               '--primer_check_source', str(ROOT), '--outdir', str(tmp_path / 'results'),
               '-work-dir', str(tmp_path / 'work'), '-ansi-log', 'false']
    if resume:
        command += ['-resume']
    result = subprocess.run(command, cwd=tmp_path, env={**os.environ, 'NXF_OFFLINE': 'true'},
                            text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / 'results/continued.txt').is_file()
    return result


@pytest.mark.parametrize('virus,targets', [
    ('SARS-CoV-2', [('SC2', '', 'SC2')]),
    ('RSV', [('R1', 'RSVA', 'R1'), ('R2', 'RSVB', 'R2'), ('RU', 'UNKNOWN', 'RU')]),
    ('Influenza', [('F1', 'H1N1', 'F1|01-HA-H1N1'), ('F3', 'H3N2', 'F3|03-MP-H3N2'),
                   ('FB', 'VICVIC', 'FB|01-HA-VICVIC')]),
])
def test_nextflow_routes_and_publishes(tmp_path, virus, targets):
    database = write_database(tmp_path / 'pcr database')
    ngs = tmp_path / 'ngs'
    write_fasta(ngs / 'primer.fasta', 'NGS_LEFT')
    for subtype in ('RSVA', 'RSVB'):
        write_fasta(ngs / subtype / 'TEST/primer.fasta', subtype + '_LEFT')
    manifest = []
    for sample, subtype, header in targets:
        # Deliberately collide basenames and include spaces: Nextflow must
        # preserve sample identity when staging inputs in numbered directories.
        fasta = write_fasta(tmp_path / sample / 'same name.fa', header)
        row = {'sample_id': sample, 'subtype': subtype, 'fasta': [fasta]}
        if virus == 'Influenza':
            subtype_file = tmp_path / sample / 'subtype.txt'
            subtype_file.write_text(subtype)
            row['subtype_file'] = str(subtype_file)
            del row['subtype']
        manifest.append(row)
    execute(tmp_path, virus, manifest, database, ngs)
    output = tmp_path / 'results/primer_check'
    rows = csv_rows(output / 'pcr_primer_report.csv')
    assert {row['Sample_ID'] for row in rows} == {s[0] for s in targets if s[1] != 'UNKNOWN'}
    assert all(row['Mismatches'] == '0' for row in rows)
    assert (output / 'pcr_primer_report.html').is_file()
    assert all(row['Status'] == 'complete' for row in csv_rows(output / 'task_status.csv'))
    assert (output / 'ngs_primer_report.csv').is_file() == (virus != 'Influenza')
    if virus == 'SARS-CoV-2':
        result = execute(tmp_path, virus, manifest, database, ngs, resume=True)
        # Mutable latest must not reuse an earlier analysis from the task cache.
        assert result.stdout.count('Submitted process > PRIMER_CHECK_RUN:PRIMER_CHECK') == 2


def test_ignored_failure_keeps_other_assay_and_pipeline(tmp_path):
    ngs = tmp_path / 'ngs'
    write_fasta(ngs / 'primer.fasta', 'NGS_LEFT')
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample')
    result = execute(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}],
                     tmp_path / 'missing-database', ngs)
    output = tmp_path / 'results/primer_check'
    assert (output / 'ngs_primer_report.csv').is_file()
    assert not (output / 'pcr_primer_report.csv').exists()
    assert 'ignored' in (result.stdout + result.stderr).lower()
    assert {row['Assay_Type']: row['Status'] for row in csv_rows(output / 'task_status.csv')} == {
        'pcr': 'failed', 'ngs': 'complete'}


def test_unified_pcr_database_and_separate_ngs_assets(tmp_path):
    database = write_mixed_database(tmp_path / 'pcr' / 'unified.json', 'virus')
    ngs = tmp_path / 'sequencing scheme'
    write_fasta(ngs / 'primer.fasta', 'NGS_LEFT')
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample')
    execute(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}],
            database.parent, ngs)
    output = tmp_path / 'results/primer_check'
    assert {row['Assay_Type']: row['Status'] for row in csv_rows(output / 'task_status.csv')} == {
        'pcr': 'complete', 'ngs': 'complete'}
    assert csv_rows(output / 'pcr_primer_report.csv')[0]['Primer_Name'] == 'F'
    assert csv_rows(output / 'ngs_primer_report.csv')[0]['Primer_Name'] == 'NGS_LEFT'


@pytest.mark.parametrize('virus,subtype,header', [
    ('SARS-CoV-2', '', 'kept'), ('RSV', 'RSVA', 'kept'),
    ('Influenza', 'H3N2', 'kept|01-HA-H3N2'),
])
def test_nextflow_checks_only_final_export(tmp_path, virus, subtype, header):
    database = write_database(tmp_path / 'pcr', unified=True)
    ngs = tmp_path / 'ngs'
    write_fasta(ngs / ('RSVA/TEST/primer.fasta' if virus == 'RSV' else 'primer.fasta'), 'NGS_LEFT')
    original = write_fasta(tmp_path / 'original.fa', header)
    filtered = write_fasta(tmp_path / 'filtered.fa', 'filtered')
    final_fasta = write_fasta(tmp_path / 'export.fasta', header, OLIGO[:18] + 'A' + OLIGO[19:])
    execute(tmp_path, virus, [
        {'sample_id': 'kept', 'subtype': subtype, 'fasta': [original]},
        {'sample_id': 'filtered', 'subtype': subtype, 'fasta': [filtered]},
    ], database, ngs, final_fasta=final_fasta)
    output = tmp_path / 'results/primer_check'
    for assay in (['pcr'] if virus == 'Influenza' else ['pcr', 'ngs']):
        rows = csv_rows(output / f'{assay}_primer_report.csv')
        assert len(rows) == 1
        assert rows[0]['Sample_ID'] == 'kept'
        assert rows[0]['Mismatches'] == '1'
        assert rows[0]['Fasta_File'] == 'export.fasta'
        assert {row['Sample_ID']: row['Status'] for row in csv_rows(output / f'{assay}_primer_report.status.csv')} == {
            'kept': 'analysed', 'filtered': 'no_consensus'}
