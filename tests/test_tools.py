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

"""gplace tog and gplace taxit against outputs of the real guppy tog (pplacer v1.1.alpha19) and
taxit (taxtastic v1.0.1), stored in tests/data/tools."""

import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from gplace.engine import testsets
from gplace.exceptions import GPlaceExit
from gplace.tools import taxit
from gplace.tools.tog import tog
from gplace.fileio import read_json

TOOLS = os.path.join(testsets.TEST_DATA, 'tools')
MINI = testsets.MINI
PARITY = testsets.PARITY


class TestTog(unittest.TestCase):
    """Newick, phyloXML and --node-numbers output byte for byte as guppy tog wrote it: queries with
    several names, ties, distal 0 / full edge, quoting, XML escaping, float-like labels."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_identical_to_guppy(self):
        for name in ('dups', 'edge_cases', 'records', 'names', 'floatlabels'):
            jp = os.path.join(TOOLS, f'{name}.jplace')
            for args, suffix in (({}, '.tog.tre'), ({'xml': True}, '.tog.xml'), ({'node_numbers': True}, '.tog_nn.tre')):
                out = os.path.join(self.tmp, name + suffix)
                tog([jp], out_file=out, **args)
                with open(out, 'rb') as a, open(os.path.join(TOOLS, name + suffix), 'rb') as b:
                    self.assertEqual(a.read(), b.read(), f'{name}{suffix}')

    def test_output_naming(self):
        jp = os.path.join(TOOLS, 'dups.jplace')
        os.makedirs(os.path.join(self.tmp, 'd'))
        self.assertEqual(tog([jp], out_dir=os.path.join(self.tmp, 'd'), prefix='p_'),
                         [os.path.join(self.tmp, 'd', 'p_dups.tog.tre')])
        self.assertEqual(tog([jp], out_dir=self.tmp, xml=True), [os.path.join(self.tmp, 'dups.tog.xml')])
        out = os.path.join(self.tmp, 'two.tre')
        tog([jp, jp], out_file=out)
        with open(out) as fh:
            self.assertEqual(len(fh.read().splitlines()), 2)       # one tree per placefile

    def test_pp_requires_posterior_probabilities(self):
        with self.assertRaises(GPlaceExit) as cm:
            tog([os.path.join(TOOLS, 'dups.jplace')], out_file=os.path.join(self.tmp, 'x'), pp=True)
        self.assertIn('Placement.No_PP', str(cm.exception))

    def test_split_by_csv(self):
        jp = os.path.join(TOOLS, 'dups.jplace')
        names = [r['nm'][0][0] for r in read_json(jp)['placements']][:6]
        csv_path = os.path.join(self.tmp, 'g.csv')
        with open(csv_path, 'w') as fh:
            fh.write(''.join(f'{n},{"b" if i % 2 else "a"}\n' for i, n in enumerate(names)))
        written = tog([f'{jp}:{csv_path}'], out_dir=self.tmp)
        self.assertEqual([os.path.basename(w) for w in written], ['b.tog.tre', 'a.tog.tre'])   # pplacer's order


class TestTaxit(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_phylo_models_identical_to_taxtastic(self):
        """RAxML, RAxML-NG, IQ-TREE, FastTree and PhyML (DNA and protein) stats files parse to the
        JSON taxtastic writes."""
        expected = read_json(os.path.join(TOOLS, 'expected_phylo_models.json'))
        self.assertEqual(len(expected), 16)
        for f, e in expected.items():
            path = os.path.join(TOOLS, f)
            self.assertEqual(taxit.guess_stats_type(path), e['stats_type'], f)
            self.assertEqual(json.dumps(taxit.phylo_model_of(None, path, e['frequency_type']), indent=4),
                             e['phylo_model'], f)

    def test_create_check_info(self):
        files = read_json(os.path.join(MINI, 'CONTENTS.json'))['files']
        pkg = os.path.join(self.tmp, 'p.refpkg')
        st = taxit.create(pkg, 'test locus', author='me', reroot=False,
                          aln_fasta=os.path.join(MINI, files['aln_fasta']), tree=os.path.join(MINI, files['tree']),
                          tree_stats=os.path.join(TOOLS, 'stats_ft_wag.log'),
                          seq_info=os.path.join(PARITY, 'mini_full.seq_info.csv'),
                          taxonomy=os.path.join(PARITY, 'mini_full.taxonomy.csv'))
        self.assertEqual(st, 0)
        c = read_json(os.path.join(pkg, 'CONTENTS.json'))
        self.assertEqual(list(c), ['metadata', 'files', 'md5', 'log', 'rollback', 'rollforward'])
        self.assertEqual(list(c['files']), ['phylo_model', 'aln_fasta', 'seq_info', 'taxonomy', 'tree', 'tree_stats'])
        self.assertEqual(c['log'], ['Stripped refpkg (removed 0 files)', 'Loaded initial files into empty refpkg'])
        self.assertEqual(c['metadata']['author'], 'me')
        model = read_json(os.path.join(pkg, c['files']['phylo_model']))
        self.assertEqual((model['subs_model'], model['ras_model']), ('WAG', 'Price-CAT'))
        out = io.StringIO()
        with redirect_stdout(out):
            self.assertEqual(taxit.check(pkg), 1)              # no aln_sto: as taxit says
        self.assertEqual(out.getvalue(), 'RefPkg has no key aln_sto\n')
        out = io.StringIO()
        with redirect_stdout(out):
            taxit.info(pkg)
        self.assertTrue(out.getvalue().startswith('number of sequences: 247\npackage components\n'))
        self.assertEqual(taxit.create(pkg, 'x'), 1)            # exists, no --clobber

    def test_md5_mismatch_is_invalid(self):
        pkg = os.path.join(self.tmp, 'p.refpkg')
        taxit.create(pkg, 'x', aln_fasta=os.path.join(MINI, 'mini_msa.fa'), reroot=False)
        with open(os.path.join(pkg, 'mini_msa.fa'), 'a') as fh:
            fh.write(' ')
        with self.assertRaises(ValueError) as cm:
            taxit.check(pkg)
        self.assertIn('did not match its MD5 sum', str(cm.exception))


if __name__ == '__main__':
    unittest.main()
