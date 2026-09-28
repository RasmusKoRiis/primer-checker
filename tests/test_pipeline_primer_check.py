"""Synthetic integration cases: real BLAST, no patient data or network required."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import shutil
import sys

import pytest

import primer_analysis as analysis

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pipeline_check', ROOT / 'scripts/run_pipeline_primer_check.py')
pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline)
OLIGO = 'GTCAGACATCGATGCTACGTCAGGATCGTACCTAGCTGAC'


def write_fasta(path, header, sequence=OLIGO):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f'>{header}\n{sequence}\n')
    return str(path)


def write_database(path):
    path.mkdir(parents=True, exist_ok=True)
    data = {
        'SARS-CoV-2': {'SARS_F': OLIGO},
        'RSV-A': {'RSVA_F': OLIGO},
        'RSV-B': {'RSVB_F': OLIGO},
        'Influenza-A': {'H1_F_HA': OLIGO, 'H3_F_HA': OLIGO, 'H3_F_M': OLIGO},
        'Influenza-B': {'B_F_HA': OLIGO},
    }
    (path / 'primers.json').write_text(json.dumps(data))
    return path


def invoke(tmp_path, virus, samples, assay='pcr', database=None, ngs=None, extra_env=None):
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(samples))
    command = [sys.executable, str(ROOT / 'scripts/run_pipeline_primer_check.py'),
               '--manifest', str(manifest), '--virus', virus, '--assay-type', assay,
               '--run-id', 'SYNTHETIC', '--output-prefix', str(tmp_path / 'result')]
    if database:
        command += ['--pcr-db', str(database)]
    if ngs:
        command += ['--ngs-dir', str(ngs), '--ngs-scheme', 'TEST']
    env = None
    if extra_env:
        import os
        env = {**os.environ, **extra_env}
    return subprocess.run(command, capture_output=True, text=True, env=env)


def csv_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


@pytest.mark.parametrize('virus,subtype,header,expected', [
    ('SARS-CoV-2', '', 'sample', 'SARS_F'),
    ('RSV', 'RSVA', 'sample', 'RSVA_F'),
    ('RSV', 'RSVB', 'sample', 'RSVB_F'),
    ('Influenza', 'H1N1', 'sample|01-HA-H1N1', 'H1_F_HA'),
    ('Influenza', 'H3N2', 'sample|01-HA-H3N2', 'H3_F_HA'),
    ('Influenza', 'H3N2', 'sample|03-MP-H3N2', 'H3_F_M'),
    ('Influenza', 'VICVIC', 'sample|01-HA-VICVIC', 'B_F_HA'),
])
def test_pcr_routes_real_blast(tmp_path, virus, subtype, header, expected):
    db = write_database(tmp_path / 'pcr with spaces')
    fasta = write_fasta(tmp_path / 'sample.fa', header)
    result = invoke(tmp_path, virus, [{'sample_id': 'sample', 'subtype': subtype, 'fasta': [fasta]}], database=db)
    assert result.returncode == 0, result.stderr
    rows = csv_rows(tmp_path / 'result.csv')
    assert {row['Primer_Name'] for row in rows} == {expected}
    assert rows[0]['Mismatches'] == '0'
    assert rows[0]['Sample_ID'] == 'sample'
    assert rows[0]['Run_ID'] == 'SYNTHETIC'
    assert (tmp_path / 'result.html').is_file()


@pytest.mark.parametrize('sequence,status', [(OLIGO, 'hit'), ('C' * 90, 'no_hit'),
                                           ('N' * 90, 'indeterminate'),
                                           (OLIGO[:18] + 'N' + OLIGO[19:], 'indeterminate')])
def test_hit_no_hit_and_masked_consensus(tmp_path, sequence, status):
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample', sequence)
    result = invoke(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}],
                    database=write_database(tmp_path / 'pcr'))
    assert result.returncode == 0, result.stderr
    row = csv_rows(tmp_path / 'result.csv')[0]
    assert row['Hit_Status'] == status
    if status == 'indeterminate':
        assert row['Percent_Identity'] == row['Mismatches'] == ''
        assert 'status_indeterminate' in (tmp_path / 'result.html').read_text()


@pytest.mark.parametrize('virus,subtype', [('SARS-CoV-2', ''), ('RSV', 'RSVA'), ('RSV', 'RSVB')])
def test_ngs_uses_actual_fasta_not_long_bed_intervals(tmp_path, virus, subtype):
    ngs = tmp_path / 'ngs'
    scheme = ngs / subtype / 'TEST' if subtype else ngs
    write_fasta(scheme / 'primer.fasta', 'oligo_LEFT', OLIGO)
    short = subtype or 'SARS-CoV-2'
    (scheme / f'{short}.primer.bed').write_text(f'ref\t0\t1000\toligo_LEFT\t1\t+\n')
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample')
    result = invoke(tmp_path, virus, [{'sample_id': 'sample', 'subtype': subtype, 'fasta': [fasta]}], 'ngs', ngs=ngs)
    assert result.returncode == 0, result.stderr
    row = csv_rows(tmp_path / 'result.csv')[0]
    assert row['Primer_Sequence'] == OLIGO
    assert row['Primer_Start'] == row['Primer_End'] == ''
    assert row['Assay_Type'] == 'ngs'
    assert row['Assay_ID'] == f'{"RSV-" + subtype[-1] if subtype else virus}:TEST'


def test_bed_sequences_and_reverse_primer(tmp_path):
    ngs = tmp_path / 'ngs'
    ngs.mkdir()
    reverse = analysis.reverse_complement(OLIGO)
    (ngs / 'SARS-CoV-2.primer.bed').write_text(f'ref\t0\t{len(OLIGO)}\toligo_RIGHT\t2\t-\t{reverse}\n')
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample')
    result = invoke(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}], 'ngs', ngs=ngs)
    assert result.returncode == 0, result.stderr
    assert csv_rows(tmp_path / 'result.csv')[0]['Mismatches'] == '0'


@pytest.mark.parametrize('virus,subtype,header,status', [
    ('RSV', 'UNKNOWN', 'sample', 'unclassified'),
    ('Influenza', 'H3N2', 'sample|06-NP-H3N2', 'no_matching_segment'),
])
def test_unmatched_inputs_still_create_reports(tmp_path, virus, subtype, header, status):
    fasta = write_fasta(tmp_path / 'sample.fa', header)
    result = invoke(tmp_path, virus, [{'sample_id': 'sample', 'subtype': subtype, 'fasta': [fasta]}],
                    database=write_database(tmp_path / 'pcr'))
    assert result.returncode == 0, result.stderr
    assert csv_rows(tmp_path / 'result.csv') == []
    assert csv_rows(tmp_path / 'result.status.csv')[0]['Status'] == status
    assert (tmp_path / 'result.html').is_file()


def test_no_assay_is_reported(tmp_path):
    db = write_database(tmp_path / 'pcr')
    (db / 'primers.json').write_text(json.dumps({'RSV-A': {'F': OLIGO}}))
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample')
    result = invoke(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}], database=db)
    assert result.returncode == 0, result.stderr
    assert csv_rows(tmp_path / 'result.status.csv')[0]['Status'] == 'no_assay'


def test_technical_blast_failure_is_not_no_hit(tmp_path):
    blast = tmp_path / 'blastn'
    blast.write_text('#!/bin/sh\necho "simulated BLAST failure" >&2\nexit 2\n')
    blast.chmod(0o755)
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample')
    result = invoke(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}],
                    database=write_database(tmp_path / 'pcr'), extra_env={'BLASTN_PATH': str(blast)})
    assert result.returncode != 0
    assert 'BLAST analysis failed' in result.stderr
    assert not (tmp_path / 'result.csv').exists()


def test_bed_without_sequences_fails_instead_of_inventing_primers(tmp_path):
    (tmp_path / 'SARS-CoV-2.primer.bed').write_text('ref\t0\t900\tp\t1\t+\n')
    with pytest.raises(ValueError, match='column 7'):
        pipeline.load_ngs_panel(tmp_path, 'SARS-CoV-2', 'TEST')


def test_empty_csv_and_write_failure(tmp_path):
    analysis.write_csv_report([], str(tmp_path / 'empty.csv'))
    assert csv_rows(tmp_path / 'empty.csv') == []
    with pytest.raises(OSError):
        analysis.write_csv_report([], str(tmp_path / 'missing' / 'result.csv'))


def test_missing_consensus_fails(tmp_path):
    result = invoke(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': ['missing.fa']}],
                    database=write_database(tmp_path / 'pcr'))
    assert result.returncode != 0
    assert 'Consensus file missing' in result.stderr


def test_known_mismatch_and_empty_run(tmp_path):
    database = write_database(tmp_path / 'pcr')
    sequence = OLIGO[:18] + ('A' if OLIGO[18] != 'A' else 'C') + OLIGO[19:]
    fasta = write_fasta(tmp_path / 'sample.fa', 'sample', sequence)
    result = invoke(tmp_path, 'SARS-CoV-2', [{'sample_id': 'sample', 'fasta': [fasta]}], database=database)
    assert result.returncode == 0, result.stderr
    row = csv_rows(tmp_path / 'result.csv')[0]
    assert row['Mismatches'] == '1'
    assert row['Mismatch_Positions'] == '19'
    result = invoke(tmp_path, 'SARS-CoV-2', [], database=database)
    assert result.returncode == 0, result.stderr
    assert csv_rows(tmp_path / 'result.csv') == []
    assert (tmp_path / 'result.html').is_file()


def test_staged_json_preserves_relative_assets(tmp_path):
    database = tmp_path / 'database'
    shutil.copytree(ROOT / 'primer_db/assets', database)
    source = database / 'panel_database.example.json'
    data = json.loads(source.read_text())
    data['viruses'][0]['pcr'] = {'schemes': [{
        'scheme_id': 'synthetic-pcr', 'display_name': 'Synthetic PCR', 'version': 'TEST',
        'assay_type': 'pcr', 'primers': [{'id': 'F', 'name': 'F', 'sequence': OLIGO,
                                        'role': 'forward_primer', 'segment': '', 'subtype_tags': []}],
    }]}
    source.write_text(json.dumps(data))
    staged = tmp_path / 'pcr_database'
    staged.symlink_to(source)
    records, _ = pipeline.load_pcr_databases(staged)
    assert records['Demo-virus'][0].sequence == OLIGO
