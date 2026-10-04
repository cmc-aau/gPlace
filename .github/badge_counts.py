"""Badge values for the README, written as GitHub Actions step outputs (key=value lines):
the unit-test counts from `python -m unittest` output, and the pplacer-parity count from the
report of `gplace test --models` (test sets compared with the stock pplacer; the
--fix_gamma_underflow set is compared with gPlace's own validated output and not counted).

Usage: python .github/badge_counts.py <unittest log> <gplace_test_report.tsv> >> "$GITHUB_OUTPUT"
"""
import csv
import os
import re
import sys


def unit_tests(path):
    if not os.path.exists(path):
        return 'not run', 'lightgrey'
    log = open(path).read()
    ran = re.search(r'^Ran (\d+) tests?', log, re.M)
    if not ran:
        return 'did not finish', 'red'
    n = int(ran.group(1))
    count = lambda key: sum(int(v) for v in re.findall(rf'\b{key}=(\d+)', log))
    bad, skipped = count('failures') + count('errors'), count('skipped')
    msg = f'{n - bad - skipped} passed'
    if skipped:
        msg += f', {skipped} skipped'
    if bad:
        msg += f', {bad} failed'
    return msg, 'red' if bad else 'brightgreen'


def parity(path):
    if not os.path.exists(path):
        return 'not run', 'lightgrey'
    with open(path, newline='') as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t') if '--fix_gamma_underflow' not in r['test_set']]
    ok = [r for r in rows if r['identical'] == r['genomes'] and r['same_tree'] == 'True']
    genomes = sum(int(r['genomes']) for r in ok)
    msg = f'{len(ok)}/{len(rows)} test sets identical ({genomes} genomes)'
    return msg, 'brightgreen' if rows and len(ok) == len(rows) else 'red'


if __name__ == '__main__':
    for key, (msg, color) in (('tests', unit_tests(sys.argv[1])), ('parity', parity(sys.argv[2]))):
        print(f'{key}_message={msg}')
        print(f'{key}_color={color}')
