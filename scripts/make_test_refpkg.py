#!/usr/bin/env python3
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

"""Build the small reference package bundled with gPlace for `gplace test`.

The mini refpkg is one clade (~250 genomes) of a GTDB-Tk R226 class-level refpkg,
with the refpkg's own FastTree WAG+CAT model. Queries are GlobDB genomes whose
GTDB-Tk placement falls inside that clade, plus a few that do not. The expected
output is produced by running GTDB-Tk's stock pplacer (v1.1.alpha19) on it:

    pplacer -m WAG -j 1 -c mini.refpkg -o expected_pplacer.jplace queries.fasta

Usage:
    make_test_refpkg.py <source refpkg dir> <GTDB-Tk export dir with jplace + msa> <out dir>
                        [--leaves 250] [--queries 16] [--others 4] [--seed 1]
"""

import argparse
import hashlib
import json
import os
import random
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gplace.engine.newick import parse_newick, quote_label  # noqa: E402


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                out.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf)))
    return out


def subtree_newick(t, v, is_root=True):
    """Newick for the clade rooted at v (no node numbers; the clade root keeps its label but
    drops its branch length)."""
    parts = []

    def rec(u, top):
        s = ''
        if t.children[u]:
            s += '(' + ','.join(rec(c, False) for c in t.children[u]) + ')'
        if t.label[u] is not None:
            s += quote_label(t.label[u])
        if not top and t.bl[u] is not None:
            s += ':' + repr(t.bl[u])
        return s

    old = sys.getrecursionlimit(); sys.setrecursionlimit(100000)
    try:
        return rec(v, is_root) + ';'
    finally:
        sys.setrecursionlimit(old)


def md5(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('refpkg'); ap.add_argument('export'); ap.add_argument('out')
    ap.add_argument('--leaves', type=int, default=250)
    ap.add_argument('--queries', type=int, default=16)
    ap.add_argument('--others', type=int, default=4)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()

    contents = json.load(open(os.path.join(a.refpkg, 'CONTENTS.json')))['files']
    t = parse_newick(open(os.path.join(a.refpkg, contents['tree'])).read())
    nleaf = [0] * t.n_nodes
    for v in t.postorder():
        nleaf[v] = 1 if not t.children[v] else sum(nleaf[c] for c in t.children[v])
    clade = min((v for v in range(t.n_nodes) if t.children[v] and v != t.root),
                key=lambda v: (abs(nleaf[v] - a.leaves), v))
    members = set()
    stack = [clade]
    while stack:
        u = stack.pop(); members.add(u); stack.extend(t.children[u])
    leaves = {t.label[u] for u in members if not t.children[u]}
    print(f'clade node {clade} ({t.label[clade]}), {len(leaves)} leaves')

    pkg = os.path.join(a.out, 'mini.refpkg')
    os.makedirs(pkg, exist_ok=True)
    open(os.path.join(pkg, 'mini_tree.tree'), 'w').write(subtree_newick(t, clade) + '\n')
    with open(os.path.join(pkg, 'mini_msa.fa'), 'w') as fh:
        for n, s in read_fasta(os.path.join(a.refpkg, contents['aln_fasta'])):
            if n in leaves:
                fh.write(f'>{n}\n{s}\n')
    shutil.copy(os.path.join(a.refpkg, contents['phylo_model']), os.path.join(pkg, 'phylo_model.json'))
    files = {'phylo_model': 'phylo_model.json', 'aln_fasta': 'mini_msa.fa', 'tree': 'mini_tree.tree'}
    json.dump({'metadata': {'format_version': '1.1',
                            'locus': f'gPlace test clade of {os.path.basename(os.path.normpath(a.refpkg))}'},
               'files': files,
               'md5': {k: md5(os.path.join(pkg, f)) for k, f in files.items()},
               'log': ['gPlace scripts/make_test_refpkg.py'], 'rollback': None, 'rollforward': None},
              open(os.path.join(pkg, 'CONTENTS.json'), 'w'), indent=4)

    # queries: best GTDB-Tk placement inside the clade, plus a few outside
    inside, outside, seqs = [], [], {}
    for root, _, fnames in sorted(os.walk(a.export)):
        if 'user_msa_file.fasta' in fnames and 'pplacer.class_level.bac120.json' in fnames:
            js = json.load(open(os.path.join(root, 'pplacer.class_level.bac120.json')))
            iE = js['fields'].index('edge_num')
            best = {pl['nm'][0][0] if 'nm' in pl else pl['n'][0]: pl['p'][0][iE] for pl in js['placements']}
            for n, s in read_fasta(os.path.join(root, 'user_msa_file.fasta')):
                if n in best:
                    seqs[n] = s
                    (inside if best[n] in members else outside).append(n)
    rnd = random.Random(a.seed)
    chosen = sorted(rnd.sample(sorted(inside), min(a.queries, len(inside)))) + \
        sorted(rnd.sample(sorted(outside), min(a.others, len(outside))))
    with open(os.path.join(a.out, 'queries.fasta'), 'w') as fh:
        for n in chosen:
            fh.write(f'>{n}\n{seqs[n]}\n')
    print(f'{len(chosen)} queries ({min(a.queries, len(inside))} placed inside the clade by GTDB-Tk)')


if __name__ == '__main__':
    main()
