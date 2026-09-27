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

"""Write a gplace place_batch batchfile for GTDB-Tk class-level alignments.

Finds every <...>/tree_<N>/<msa_name> below --msa_root and pairs it with
<refpkg_root>/gtdbtk.package.<N>.refpkg. The output jplace is written next to the alignment
(--out_name), or into the same relative path below --out_root.

Example:
    make_gtdbtk_batchfile.py --msa_root runs --refpkg_root release226/split/class_level/pplacer > batch.tsv
"""

import argparse
import os
import re
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--msa_root', required=True, help='directory searched for tree_<N>/<msa_name>')
    ap.add_argument('--refpkg_root', required=True, help='directory with gtdbtk.package.<N>.refpkg')
    ap.add_argument('--msa_name', default='user_msa_file.fasta', help='alignment file name')
    ap.add_argument('--out_name', default='pplacer.class_level.bac120.json', help='output jplace file name')
    ap.add_argument('--out_root', default=None, help='write outputs below this directory instead of next to the alignments')
    a = ap.parse_args()
    n = 0
    for root, _, files in sorted(os.walk(a.msa_root)):
        m = re.search(r'tree_(\d+)$', root)
        if not m or a.msa_name not in files:
            continue
        refpkg = os.path.join(a.refpkg_root, f'gtdbtk.package.{m.group(1)}.refpkg')
        if not os.path.isfile(os.path.join(refpkg, 'CONTENTS.json')):
            sys.exit(f'No reference package for {root}: {refpkg}')
        out_dir = root if a.out_root is None else os.path.join(a.out_root, os.path.relpath(root, a.msa_root))
        print(f'{os.path.abspath(os.path.join(root, a.msa_name))}\t{os.path.abspath(refpkg)}\t'
              f'{os.path.abspath(os.path.join(out_dir, a.out_name))}')
        n += 1
    print(f'{n} alignments', file=sys.stderr)


if __name__ == '__main__':
    main()
