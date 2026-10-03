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

"""``guppy tog`` (pplacer v1.1.alpha19, guppy_tog.ml): a tree with every placed query added as a
pendant edge at its best placement, written byte for byte as guppy writes it (Newick or, with
``--xml``, phyloXML).

What guppy does, and so what this does:

* each jplace is read as pplacer reads it (Placerun_io.of_json_file): a query's placements are
  sorted by like-weight ratio, highest first (stable), and queries without placements are dropped;
* the "best" placement is the FIRST of that list, whatever the criterion: pplacer's
  Pquery.opt_best_something compares ``criterion h`` (the first) instead of ``criterion x``, so
  ``--pp`` only requires posterior probabilities to be present (Placement.No_PP otherwise);
* every query becomes a leaf (its name, branch length = pendant length) attached at its distal
  length on its edge; a query with several names becomes a node (pendant length) with one
  zero-length leaf per name (Gtree.add_subtrees_by_map). Queries on one edge are attached in
  increasing distal length (stable sort), each one cutting the remaining edge: the lower part
  gets (distal - previous distal), the new node above it the rest (0 if within 1e-10 below 0);
* the new node's children are (the edge below, the query); branch lengths are printed with %g,
  labels quoted as pplacer quotes them, ``--node-numbers`` turns the reference tree's labels into
  ``@label`` and puts its node numbers where bootstraps go (``[n]``).
"""

import csv
import gzip
import json
import os
import re

from gplace.engine.newick import _tokens, parse_newick, quote_label
from gplace.exceptions import GPlaceExit

BL_TOL = 1e-10                       # Gtree.bl_arithmetic_tol


class PlaceTree(object):
    """A tree under construction: per node its children, bark (bl, label, edge_label, red decor)
    and whether it has bark at all (pplacer prints nothing for a node without bark)."""

    def __init__(self):
        self.children, self.bl, self.label, self.edge_label, self.red, self.has_bark = [], [], [], [], [], []

    def add(self, children, bl=None, label=None, edge_label=None, red=False, has_bark=True):
        for a, v in ((self.children, children), (self.bl, bl), (self.label, label),
                     (self.edge_label, edge_label), (self.red, red), (self.has_bark, has_bark)):
            a.append(v)
        return len(self.children) - 1


# ----------------------------------------------------------------------------- reading ---

def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def placerun_name(path):
    """MaybeZipped.chop_extension (Filename.basename path): drop .gz, then the last extension."""
    base = os.path.basename(path)
    if base.endswith('.gz'):
        base = base[:-3]
    root, ext = os.path.splitext(base)
    return root if ext else base


def read_placerun(path):
    """Placerun_io.of_any_file for a .jplace / .json file: (name, reference tree, number map,
    pqueries); a pquery is (names, placements sorted by like-weight ratio)."""
    stem = path[:-3] if path.endswith('.gz') else path
    if stem.endswith('.place'):
        raise GPlaceExit(f'{path}: the old .place format is not supported; use a .jplace file.')
    if not (stem.endswith('.jplace') or stem.endswith('.json')):
        raise GPlaceExit(f'unfamiliar suffix on {path}')
    with _open(path) as fh:
        j = json.load(fh)
    if 'version' not in j:
        raise GPlaceExit(f"{path}: invalid placerun: no 'version' field")
    if j['version'] not in (2, 3):
        raise GPlaceExit(f"{path}: jplace version {j['version']} is not supported (versions 2 and 3 are).")
    tree = parse_newick(j['tree'])
    # the {n} numbers of a jplace tree, in text (= post-order) order, are pplacer's node ids
    toks = list(_tokens(j['tree'].strip()))
    nums = [int(toks[k + 1][1]) for k in range(len(toks) - 2) if toks[k][0] == '{' and toks[k + 2][0] == '}']
    if len(nums) != tree.n_nodes:
        raise GPlaceExit(f'{path}: the jplace tree does not number every edge.')
    id_of = {n: v for v, n in enumerate(nums)}
    fields = j['fields']
    pqueries = []
    for rec in j['placements']:
        if 'n' in rec:
            n = rec['n']
            names = [n] if isinstance(n, str) else list(n)
        else:
            names = [nm[0] for nm in rec['nm']]
        places = [dict(zip(fields, p)) for p in rec['p']]
        places.sort(key=lambda p: p['like_weight_ratio'], reverse=True)       # stable, highest first
        pqueries.append((names, places))
    if not pqueries:
        raise GPlaceExit(f'{path} has no placements!')
    return placerun_name(path), tree, id_of, pqueries


def split_placerun(arg):
    """Placerun_io.of_split_file for ``placefile:groups.csv`` (rows: query name, group): one
    placerun per group, named after it, as pplacer builds them: groups in decreasing order, and in
    each group the queries, and a query's names, in reverse order."""
    m = re.match(r'^(.+):(.+?)$', arg, re.S)
    if not m:
        raise GPlaceExit('csv file provided with no jplace to split')
    _, tree, id_of, pqueries = read_placerun(m.group(1))
    group_of = {}
    with open(m.group(2), newline='') as fh:
        for row in csv.reader(fh):
            if len(row) != 2:
                raise GPlaceExit('malformed row in csv file')
            group_of[row[0]] = row[1]
    split = {}
    for names, places in pqueries:
        groups = {}
        for nm in names:
            if nm in group_of:
                groups.setdefault(group_of[nm], []).insert(0, nm)
        for g in sorted(groups):
            split.setdefault(g, []).insert(0, (groups[g], places))
    return [(g, tree, id_of, split[g]) for g in sorted(split, reverse=True)]


def placeruns(arg):
    """Placerun_io.maybe_of_split_file: a placefile, or ``placefile:groups.csv``."""
    if arg.endswith('.csv'):
        return split_placerun(arg)
    if arg.endswith('.biom'):
        raise GPlaceExit(f'{arg}: BIOM input is not supported.')
    return [read_placerun(arg)]


# ---------------------------------------------------------------------------- the tree ---

def _has_bark(tree, v):
    return tree.bl[v] is not None or tree.label[v] is not None or tree.edge_label[v] is not None


def tog_tree(tree, id_of, pqueries, numbered=False):
    """Guppy_tog.tog_tree on the reference tree; returns (PlaceTree, root)."""
    by_loc = {}
    for names, places in pqueries:
        if not places:
            continue                                     # unplaced queries are left out
        best = places[0]
        by_loc.setdefault(id_of[best['edge_num']], []).append((best['distal_length'], names, best['pendant_length']))
    out = PlaceTree()
    mapped = {}
    for v in tree.postorder():
        bark = _has_bark(tree, v)
        label, edge_label = tree.label[v], tree.edge_label[v]
        if numbered and bark:                            # Newick_bark.to_numbered: '@label', or '' if none
            label, edge_label = ('@' + label if label is not None else ''), str(v)
        top = out.add([mapped.pop(c) for c in tree.children[v]], tree.bl[v], label, edge_label, False, bark)
        pos = 0.0
        for new_pos, names, pendant in sorted(by_loc.get(v, []), key=lambda a: a[0]):
            where = new_pos - pos
            new_top_bl = out.bl[top] - where
            if new_top_bl < 0.0:
                if abs(new_top_bl) < BL_TOL:
                    new_top_bl = 0.0
                else:
                    raise GPlaceExit('Internal edge became %g when attaching subtrees to %d' % (new_top_bl, v))
            out.bl[top] = where
            out.has_bark[top] = True
            if len(names) == 1:
                sub = out.add([], pendant, names[0], None, True)
            else:
                leaves = [out.add([], 0.0, nm, None, True) for nm in names]
                sub = out.add(leaves, pendant, None, None, True)
            top = out.add([top, sub], new_top_bl, None, None, False)
            pos = new_pos
        mapped[v] = top
    return out, mapped[tree.root]


# ------------------------------------------------------------------------------ output ---

def _bark_newick(t, v):
    """Newick_bark.to_newick_string (no node numbers: tog writes plain Newick)."""
    if not t.has_bark[v]:
        return ''
    s = quote_label(t.label[v]) if t.label[v] is not None else ''
    if t.bl[v] is not None:
        s += ':%g' % t.bl[v]
    if t.edge_label[v] is not None:
        s += '[%s]' % quote_label(t.edge_label[v])
    return s


def _postorder(t, root):
    order, stack = [], [(root, False)]
    while stack:
        v, done = stack.pop()
        if done or not t.children[v]:
            order.append(v)
            continue
        stack.append((v, True))
        for c in reversed(t.children[v]):
            stack.append((c, False))
    return order


def to_newick(t, root):
    out = {}
    for v in _postorder(t, root):
        kids = t.children[v]
        out[v] = ('(' + ','.join(out.pop(c) for c in kids) + ')' if kids else '') + _bark_newick(t, v)
    return out[root] + ';'


def _ocaml_float(s):
    """OCaml float_of_string, as Newick_bark.maybe_float uses it (None if not a float): every '_' is
    dropped, then hexadecimal ([+-]0x...) or C strtod, which skips leading blanks but must consume
    the rest of the string."""
    if s is None:
        return None
    s = s.replace('_', '')
    if not s.strip() or s != s.rstrip():
        return None
    try:
        if re.match(r'\s*[+-]?0[xX]', s):
            return float.fromhex(s.strip())
        return float(s)
    except ValueError:
        return None


def _xml_escape(s):
    """Character data as xmlm writes it."""
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')


def _clade_xml(t, v, parts):
    parts.append('<clade>')
    if t.has_bark[v]:
        is_leaf = not t.children[v]
        conf = _ocaml_float(t.edge_label[v])              # Newick_bark.get_confidence_name_opt
        name = t.label[v]
        if conf is None and not is_leaf:
            c = _ocaml_float(t.label[v])
            if c is not None:
                conf, name = c, None
        if name is not None:                              # xmlm writes an empty element as <name/>
            parts.append('<name>%s</name>' % _xml_escape(name) if name else '<name/>')
        if t.bl[v] is not None:
            parts.append('<branch_length>%g</branch_length>' % t.bl[v])
        if conf is not None:
            parts.append('<confidence type="confidence">%g</confidence>' % conf)
        if t.red[v]:
            parts.append('<color><red>255</red><green>0</green><blue>0</blue></color>')
    for c in t.children[v]:
        _clade_xml(t, c, parts)
    parts.append('</clade>')


def to_phyloxml(trees):
    """Phyloxml.gtrees_to_file (xmlm output: no indentation, no final newline)."""
    import sys
    old = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old, 100000))
    try:
        parts = ['<?xml version="1.0" encoding="UTF-8"?>\n',
                 '<phyloxml xmlns="http://www.phyloxml.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                 'xsi:schemaLocation="http://www.phyloxml.org/1.10/phyloxml.xsd">']
        for t, root in trees:
            parts.append('<phylogeny rooted="true">')
            _clade_xml(t, root, parts)
            parts.append('</phylogeny>')
        parts.append('</phyloxml>')
    finally:
        sys.setrecursionlimit(old)
    return ''.join(parts)


# ----------------------------------------------------------------------------- command ---

def out_target(out_file=None, out_dir=None, prefix=None):
    """Guppy_cmdobjs.output_cmd#out_file_or_dir: ('file', path) or ('dir', dir, prefix)."""
    if out_file is None:
        return ('dir', out_dir if out_dir is not None else '.', prefix or '')
    name = (prefix or '') + out_file
    return ('file', os.path.join(out_dir, name) if out_dir is not None else name)


def tog(placefiles, out_file=None, out_dir=None, prefix=None, pp=False, xml=False, node_numbers=False):
    """guppy tog [--pp] [--xml] [--node-numbers] [-o FILE | --out-dir DIR] [--prefix P] placefile[s]
    (a placefile may be ``x.jplace:groups.csv``). Returns the list of files written."""
    named = []
    for name, tree, id_of, pqueries in [pr for arg in placefiles for pr in placeruns(arg)]:
        if pp:
            for names, places in pqueries:
                if places and places[0].get('post_prob') is None:
                    raise GPlaceExit('no posterior probability present (pplacer: Fatal error: '
                                     'exception Placement.No_PP)')
        named.append((name, tog_tree(tree, id_of, pqueries, numbered=node_numbers)))
    ext = '.xml' if xml else '.tre'
    render = (lambda trees: to_phyloxml(trees)) if xml else \
        (lambda trees: ''.join(to_newick(t, r) + '\n' for t, r in trees))
    target = out_target(out_file, out_dir, prefix)
    written = []
    if target[0] == 'file':
        jobs = [(target[1], [tr for _, tr in named])]
    else:
        base = os.path.join(target[1], target[2])          # Filename.concat dir prefix
        jobs = [(base + name + '.tog' + ext, [tr]) for name, tr in named]
    for fname, trees in jobs:
        with open(fname, 'w') as fh:
            fh.write(render(trees))
        written.append(fname)
    return written
