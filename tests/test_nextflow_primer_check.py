"""Run the actual shared Nextflow subworkflow on synthetic consensus data."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from test_pipeline_primer_check import ROOT, csv_rows, write_database, write_fasta

pytestmark = pytest.mark.skipif(not shutil.which('nextflow'), reason='Nextflow is not installed')
HARNESS = ROOT / 'integrations/nextflow/tests/primer_check'


def execute(tmp_path, virus, manifest, database, ngs, resume=False):
    path = tmp_path / 'manifest.json'
    path.write_text(json.dumps(manifest))
    command = ['nextflow', '-C', str(HARNESS / 'nextflow.config'), 'run', str(HARNESS / 'main.nf'),
               '--test_manifest', str(path), '--test_virus', virus,
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
