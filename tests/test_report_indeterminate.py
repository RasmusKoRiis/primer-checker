"""Execute the report's JavaScript for uncertain primer-binding evidence."""
import re
import shutil
import subprocess

import pytest

import primer_report

pytestmark = pytest.mark.skipif(not shutil.which('node'), reason='Node.js is not installed')


def test_indeterminate_rows_are_not_confirmed_hits(tmp_path):
    html = primer_report.build_html_report([])
    script = re.findall(r'<script>(.*?)</script>', html, re.S)[-1]
    script_file = tmp_path / 'report.js'
    script_file.write_text(script)
    subprocess.run(['node', '--check', str(script_file)], check=True, capture_output=True, text=True)
    names = ['parseNumber', 'isNoHit', 'statusLabel', 'calculateMismatchDistribution',
             'calculatePrimerStats', 'calculateRisk']
    functions = '\n'.join(re.search(r'    function ' + name + r'\(.*?\n    }', script, re.S)[0]
                          for name in names)
    thresholds = re.search(r'    const RISK_THRESHOLDS = .*?\n    };', script, re.S)[0]
    test = functions + '\n' + thresholds + '''
      const assert = require('node:assert/strict');
      const t = key => key;
      const roleForRow = () => 'primer';
      const calculateMismatchPositionCounts = () => ({counts: []});
      const parseMismatchPositions = () => [];
      const terminalMismatchPosition = () => false;
      const unknown = {Hit_Status: 'indeterminate', Mismatches: '', Percent_Identity: ''};
      const confirmed = {Hit_Status: 'hit', Mismatches: 0, Percent_Identity: 100};
      const noHit = {Hit_Status: 'no_hit', Mismatches: '', Percent_Identity: 'No hit'};
      assert.equal(statusLabel(unknown), 'status_indeterminate');
      let stats = calculatePrimerStats([unknown]);
      assert.equal(stats.hits, 0);
      assert.equal(stats.noHits, 0);
      assert.equal(stats.indeterminate, 1);
      assert.equal(calculateRisk(stats), 'Indeterminate');
      stats = calculatePrimerStats([unknown, confirmed]);
      assert.equal(stats.hits, 1);
      assert.equal(stats.avgPercentIdentity, 100);
      stats = calculatePrimerStats([unknown, confirmed, noHit]);
      assert.equal(stats.noHitRate, 0.5);
      assert.equal(stats.distribution.indeterminate, 1);
      assert.equal(stats.distribution['0'], 1);
    '''
    subprocess.run(['node', '-e', test], check=True, capture_output=True, text=True)
