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
non-finite values are written as ``nan``/``inf`` (strictly not valid JSON); :func:`read_jplace`
accepts them.
"""

import json
import re

from gplace.config.output import JPLACE_FIELDS

_JSON_ESCAPE = re.compile(r'([\\"/\b\f\n\r\t])')
_ESC = {'\\': '\\\\', '"': '\\"', '/': '\\/', '\b': '\\b', '\f': '\\f', '\n': '\\n', '\r': '\\r', '\t': '\\t'}


def _qs(s):
    """Json.quote from pplacer (also escapes '/')."""
    return '"' + _JSON_ESCAPE.sub(lambda m: _ESC[m.group(1)], s) + '"'


def _f(v):
    return '%.12g' % v


def write_jplace(path, tree, results, invocation):
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
    """
    with open(path, 'w') as fh:
        fh.write('{"tree": ' + _qs(tree) + ', "placements":\n  [\n')
        for i, (name, pl, _) in enumerate(results):
            ps = ', '.join('[%s, %d, %s, %s, %s]' % (_f(p[0]), int(p[1]), _f(p[2]), _f(p[3]), _f(p[4]))
                           for p in pl)
            fh.write('  {"p": [%s], "nm": [[%s, 1]]}%s\n' % (ps, _qs(name), ',' if i + 1 < len(results) else ''))
        fh.write('  ], "metadata": {"invocation": ' + _qs(invocation) + '}, "version": 3, "fields":\n')
        fh.write('  [' + ', '.join(_qs(f) for f in JPLACE_FIELDS) + ']\n}\n')


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


def compare_jplace(path_a, path_b):
    """Compare two jplace files placement by placement.

    Returns a dict with counts of identical best edges, identical full edge lists, identical
    printed values (every field equal) and the largest absolute difference per field for the
    best placement, plus up to 25 examples of differing genomes."""
    tree_a, a, _ = read_jplace(path_a)
    tree_b, b, _ = read_jplace(path_b)
    names = sorted(set(a) | set(b))
    rep = {'genomes': len(names), 'only_a': 0, 'only_b': 0, 'best_edge': 0, 'edge_list': 0,
           'identical': 0, 'same_tree': tree_a == tree_b,
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
