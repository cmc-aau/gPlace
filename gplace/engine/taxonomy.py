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

"""Taxonomic classification of placements, as pplacer v1.1.alpha19 adds it when the reference
package has a taxonomy and seq_info (``Refpkg.tax_equipped``): every placement gets a
``classification`` field, the tax_id of its edge (or null).

Two ways to classify an edge, as in pplacer:
  * edge painting (default; ``edge_painting.ml``): for the most specific rank first, the edge
    gets a taxon if exactly one taxon of that rank is *cut* by the edge (``Convex``), otherwise
    the next rank up; NoTax if no rank has exactly one;
  * MRCA (``--mrca-class``; ``tax_map.ml``, ``classif_map.ml``): the taxonomic MRCA of the
    reference sequences below the edge.

Files (``tax_taxonomy.ml``, ``tax_seqinfo.ml``): the taxonomy is pplacer's tax table CSV,
``tax_id,parent_id,rank,tax_name`` then one column per rank holding the lineage's tax_ids; the
header's rank columns name the ranks and a tax_id's rank is its lineage column. Parents come from
consecutive non-empty lineage columns (not from parent_id). seq_info needs ``seqname`` and
``tax_id`` columns. Empty and ``NA`` cells are missing values.
"""

import csv
import json
import os

from gplace.exceptions import RefpkgMalformed

NOTAX = None


def tax_equipped(refpkg_dir):
    files = json.load(open(os.path.join(refpkg_dir, 'CONTENTS.json'))).get('files', {})
    return 'taxonomy' in files and 'seq_info' in files


def _entry(s):
    return None if s in ('', 'NA') else s


def _rows(path):
    with open(path, newline='') as fh:
        return [row for row in csv.reader(fh)]


class Taxonomy(object):
    """Tax_taxonomy.t: ancestor map, rank index of every tax_id, rank names."""

    def __init__(self, path):
        rows = _rows(path)
        if not rows:
            raise RefpkgMalformed(f'Empty taxonomy: {path}')
        if any(len(r) != len(rows[0]) for r in rows):
            raise RefpkgMalformed(f'Taxonomy array not rectangular: {path}')
        if len(rows[0]) < 4:
            raise RefpkgMalformed(f'This tax line didn\'t fit expectations: {rows[0]}')
        self.rank_names = rows[0][4:]
        if any(_entry(r) is None for r in self.rank_names):
            raise RefpkgMalformed('NA in taxon rank name line!')
        self.parent, self.rank = {}, {}
        for r in rows[1:]:
            lin = [(i, _entry(x)) for i, x in enumerate(r[4:])]
            lin = [(i, x) for i, x in lin if x is not None]
            for (i, x), (_, y) in zip(lin, lin[1:]):          # y's parent is x
                self._check(self.parent, y, x, 'parent')
                self._check(self.rank, x, i, 'rank')
            if lin:
                self._check(self.rank, lin[-1][1], lin[-1][0], 'rank')

    @staticmethod
    def _check(m, k, v, what):
        if k in m and m[k] != v:
            if what == 'parent':
                raise RefpkgMalformed(f'Tax table broken: tax_id {k} had established parent {m[k]} '
                                      f'but {v} is claiming to be the parent.')
            raise RefpkgMalformed(f'Tax table broken: {k} had established rank {m[k]} but was '
                                  f'also found at rank {v}.')
        m[k] = v

    def get_rank(self, ti):
        if ti is NOTAX:
            return 0
        if ti not in self.rank:
            raise RefpkgMalformed(f'Tax_taxonomy.get_tax_rank not known: {ti}')
        return self.rank[ti]

    def lineage(self, ti):
        """Root first, ending with ti; [] for NoTax."""
        if ti is NOTAX:
            return []
        out = [ti]
        while out[-1] in self.parent:
            out.append(self.parent[out[-1]])
        return out[::-1]

    def mrca(self, a, b):
        """Tax_taxonomy.mrca: walk the deeper-ranked taxon up until both meet."""
        a0, b0 = a, b
        while True:
            if b is NOTAX:
                return a
            if a is NOTAX:
                return b
            if a == b:
                return a
            ra, rb = self.get_rank(a), self.get_rank(b)
            prox, dist = (a, b) if ra < rb else (b, a)            # sort_by_rank
            if dist not in self.parent:
                raise RefpkgMalformed(f'No MRCA for tax_ids {a0} and {b0} in the taxonomy.')
            a, b = self.parent[dist], prox


def read_seqinfo(path):
    """Tax_seqinfo.of_csv: seqname -> tax_id (None for empty/NA)."""
    rows = _rows(path)
    if not rows:
        raise RefpkgMalformed('empty refdata in seqinfo')
    head = rows[0]
    for col in ('seqname', 'tax_id'):
        if col not in head:
            raise RefpkgMalformed(f'couldn\'t find {col} column in seqinfo')
    si, ti = head.index('seqname'), head.index('tax_id')
    out = {}
    for r in rows[1:]:
        d = dict(zip(head, r))
        name, tid = r[si], _entry(d.get('tax_id', ''))
        if name in out:
            raise RefpkgMalformed(f'Tax_refdata.of_csv: contradictory line for {name}')
        out[name] = tid
    return out


def _postorder(children, root):
    order, stack = [], [(root, False)]
    while stack:
        v, done = stack.pop()
        if done or not children[v]:
            order.append(v)
            continue
        stack.append((v, True))
        for c in reversed(children[v]):
            stack.append((c, False))
    return order


def _cutsets(children, root, colors):
    """Convex.build_sizemim_and_cutsetim (cut sets only): node -> set of colors cut by the edge
    above it. colors: leaf node -> color (a tax_id at one rank)."""
    below = {}
    for v in _postorder(children, root):
        if not children[v]:
            below[v] = {colors[v]} if v in colors else set()
        else:
            s = set()
            for c in children[v]:
                s |= below[c]
            below[v] = s
    cut = dict(below)
    stack = [(root, frozenset())]
    while stack:
        v, terminated = stack.pop()
        ch = children[v]
        if not ch:
            continue
        sets = [cut[c] for c in ch]
        between = set()                                   # colors in >= 2 of the subtrees
        for i in range(len(sets)):
            for j in range(i + 1, len(sets)):
                between |= sets[i] & sets[j]
        term = terminated | between
        for c, s in zip(ch, sets):
            s2 = s & term
            if s2 != s:                                   # else nothing below c is cut further
                cut[c] = s2
                stack.append((c, term))
    return cut


def edge_painting(tax, seqinfo, label, children, root):
    """Edge_painting.of_refpkg: (node -> tax_id or NoTax, ranks present among the reference
    sequences in the tree)."""
    leaf_of = {label[v]: v for v in range(len(children)) if not children[v] and label[v] is not None}
    rankmap = {}                                          # rank -> {leaf node: tax_id}
    for name in sorted(seqinfo):                          # StringMap order
        node = leaf_of.get(name)
        if node is None:
            continue
        for t in tax.lineage(seqinfo[name]):
            rankmap.setdefault(tax.get_rank(t), {})[node] = t
    if not rankmap:
        return {}, set()
    cut = {r: _cutsets(children, root, colors) for r, colors in rankmap.items()}
    highest = max(rankmap)
    paint = {}
    for v in range(len(children)):
        r = highest
        while True:
            if r not in cut:
                if r <= 0:
                    raise RefpkgMalformed('Edge painting reached below the top rank of the taxonomy.')
                r -= 1
                continue
            cs = cut[r][v]
            if len(cs) == 1:
                paint[v] = next(iter(cs)); break
            if r == 0:
                paint[v] = NOTAX; break
            r -= 1
    return paint, set(rankmap)


def mrca_map(tax, seqinfo, label, children, root):
    """Classif_map with --mrca-class: every node maps to the taxonomic MRCA of the reference
    sequences below it (Tax_map.fill_out; all_mrcas fills every node with that value)."""
    full = {}
    for v in _postorder(children, root):
        if not children[v]:
            if label[v] is None or label[v] not in seqinfo:
                raise RefpkgMalformed(f'tax_id for {label[v]} not found!')
            full[v] = seqinfo[label[v]]
        else:
            ts = [full[c] for c in children[v]]
            m = ts[0]
            for t in ts[1:]:
                m = tax.mrca(m, t)
            full[v] = m
    return full


class Classifier(object):
    """Edge number -> tax_id (str) or None, for a tax-equipped reference package."""

    def __init__(self, refpkg_dir, tree, mrca_class=False):
        files = json.load(open(os.path.join(refpkg_dir, 'CONTENTS.json')))['files']
        tax = Taxonomy(os.path.join(refpkg_dir, files['taxonomy']))
        seqinfo = read_seqinfo(os.path.join(refpkg_dir, files['seq_info']))
        for name, t in seqinfo.items():
            if t is not NOTAX and t not in tax.rank:
                raise RefpkgMalformed(f'Tax_taxonomy.get_tax_rank not known: {t} (seq_info row {name})')
        self.method = 'MRCA' if mrca_class else 'edge painting'
        self.missing_ranks = []
        if mrca_class:
            self.map = mrca_map(tax, seqinfo, tree.label, tree.children, tree.root)
        else:
            self.map, present = edge_painting(tax, seqinfo, tree.label, tree.children, tree.root)
            self.missing_ranks = [r for i, r in enumerate(tax.rank_names) if i not in present]

    def __call__(self, edge):
        return self.map.get(edge, NOTAX)
