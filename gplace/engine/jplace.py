###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""jplace (version 3) files written the way pplacer v1.1.alpha19 writes them.

pplacer's Json module prints floats with ``%.12g`` and integers with ``%d``; with gPlace's
bit-identical doubles every number in the file is textually identical to pplacer's output.
Whitespace and key order are JSON-insignificant and not reproduced byte-for-byte. Like pplacer,
non-finite values are written as ``-nan``/``nan``/``-inf``/``inf`` (strictly not valid JSON);
:func:`read_jplace` accepts them.
"""

import json
import math
import re

from gplace.config.output import JPLACE_FIELDS

_JSON_ESCAPE = re.compile(r'([\\"/\b\f\n\r\t])')
_ESC = {'\\': '\\\\', '"': '\\"', '/': '\\/', '\b': '\\b', '\f': '\\f', '\n': '\\n', '\r': '\\r', '\t': '\\t'}


def _qs(s):
    """Json.quote from pplacer (also escapes '/')."""
    return '"' + _JSON_ESCAPE.sub(lambda m: _ESC[m.group(1)], s) + '"'


def _f(v):
    """OCaml's %.12g via glibc printf: a NaN with the sign bit set prints as -nan (pplacer's
    like-weight ratios after a -inf likelihood are such NaNs)."""
    if v != v:
        return '-nan' if math.copysign(1.0, v) < 0 else 'nan'
    return '%.12g' % v


def write_jplace(path, tree, results, invocation, records=None, classify=None):
    """Write a jplace file.

    Parameters
    ----------
    path : str
        Output file.
    tree : str
        Reference tree with {edge} numbers (Placer.jplace_tree()).
    results : list
        [(name, placements, warnings)], placements as [distal, edge, lwr, likelihood, pendant].
    invocation : str
        Command line recorded in the metadata.
    records : list, optional
        [[name, ...], ...]: names sharing a record (pplacer writes identical query sequences as one
        record, Placer.records); default one record per result.
    classify : callable, optional
        edge number -> tax_id or None (taxonomy.Classifier): adds pplacer's "classification"
        field, first in the sorted field list.
    """
    pl_of = {name: pl for name, pl, _ in results}
    if records is None:
        records = [[name] for name, _, _ in results]
    fields = (['classification'] if classify is not None else []) + JPLACE_FIELDS
    with open(path, 'w') as fh:
        fh.write('{"tree": ' + _qs(tree) + ', "placements":\n  [\n')
        for i, names in enumerate(records):
            rows = []
            for p in pl_of[names[0]]:
                row = '%s, %d, %s, %s, %s' % (_f(p[0]), int(p[1]), _f(p[2]), _f(p[3]), _f(p[4]))
                if classify is not None:
                    tax = classify(int(p[1]))
                    row = ('null' if tax is None else _qs(tax)) + ', ' + row
                rows.append('[' + row + ']')
            nm = ', '.join('[%s, 1]' % _qs(n) for n in names)
            fh.write('  {"p": [%s], "nm": [%s]}%s\n' % (', '.join(rows), nm, ',' if i + 1 < len(records) else ''))
        fh.write('  ], "metadata": {"invocation": ' + _qs(invocation) + '}, "version": 3, "fields":\n')
        fh.write('  [' + ', '.join(_qs(f) for f in fields) + ']\n}\n')


_NONFINITE = [(re.compile(r'(?<![\w."])-?nan(?![\w"])'), 'NaN'),
              (re.compile(r'(?<![\w."])-inf(?![\w"])'), '-Infinity'),
              (re.compile(r'(?<![\w."])inf(?![\w"])'), 'Infinity')]


def read_jplace(path):
    """Parse a jplace file (tolerating pplacer's nan/inf) into
    (tree, {name: [ {field: value, ...}, ... ]}, fields)."""
    text = open(path).read()
    for rx, rep in _NONFINITE:
        text = rx.sub(rep, text)
    js = json.loads(text)
    fields = js['fields']
    out = {}
    for pl in js['placements']:
        names = [n[0] for n in pl['nm']] if 'nm' in pl else list(pl['n'])
        rows = [dict(zip(fields, p)) for p in pl['p']]
        for n in names:
            out[n] = rows
    return js.get('tree', ''), out, fields


def _records(path):
    """The name lists of a file's records (pplacer writes identical queries as one record)."""
    text = open(path).read()
    for rx, rep in _NONFINITE:
        text = rx.sub(rep, text)
    return sorted(tuple(n[0] for n in pl['nm']) if 'nm' in pl else tuple(pl['n'])
                  for pl in json.loads(text)['placements'])


def compare_jplace(path_a, path_b):
    """Compare two jplace files placement by placement.

    Returns a dict with counts of identical best edges, identical full edge lists, identical
    printed values (every field equal) and the largest absolute difference per field for the
    best placement, plus up to 25 examples of differing genomes."""
    tree_a, a, fields_a = read_jplace(path_a)
    tree_b, b, fields_b = read_jplace(path_b)
    names = sorted(set(a) | set(b))
    rep = {'genomes': len(names), 'only_a': 0, 'only_b': 0, 'best_edge': 0, 'edge_list': 0,
           'identical': 0, 'same_tree': tree_a == tree_b,
           'same_records': _records(path_a) == _records(path_b), 'same_fields': fields_a == fields_b,
           'max_delta': {k: 0.0 for k in ('likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length')},
           'examples': []}
    for n in names:
        if n not in a:
            rep['only_b'] += 1; continue
        if n not in b:
            rep['only_a'] += 1; continue
        pa, pb = a[n], b[n]
        best = pa[0]['edge_num'] == pb[0]['edge_num']
        rep['best_edge'] += best
        same_list = [r['edge_num'] for r in pa] == [r['edge_num'] for r in pb]
        rep['edge_list'] += same_list
        rep['identical'] += pa == pb
        if best:
            for k in rep['max_delta']:
                rep['max_delta'][k] = max(rep['max_delta'][k], abs(pa[0][k] - pb[0][k]))
        if not same_list and len(rep['examples']) < 25:
            rep['examples'].append((n, [(r['edge_num'], r['like_weight_ratio']) for r in pa],
                                    [(r['edge_num'], r['like_weight_ratio']) for r in pb]))
    return rep
