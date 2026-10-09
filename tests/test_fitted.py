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

"""Trees fitted with IQ-TREE and RAxML-NG (tests/data/fits): the model pplacer places with, what
it drops (warnings), why it cannot use some fits (errors), and placements identical to stock pplacer
on the package taxit builds from the same files."""

import os
import shutil
import tempfile
import unittest

from gplace.engine import core, fitted, testsets
from gplace.engine.factory import new_placer, tiled_placer
from gplace.engine.jplace import write_jplace, compare_jplace
from gplace.engine.reference import Reference
from gplace.exceptions import RefpkgMalformed

FITS = os.path.join(testsets.TEST_DATA, 'fits')
PROT_ALN = os.path.join(testsets.MINI, 'mini_msa.fa')
DNA_ALN = os.path.join(testsets.MODELS, 'like_gtr.refpkg', 'aln.fasta')
QUERIES = os.path.join(testsets.TEST_DATA, 'queries.fasta')
GTR_Q = os.path.join(testsets.MODELS, 'like_gtr_queries.fasta')


def fit(name):
    return os.path.join(FITS, name)


class TestReadFit(unittest.TestCase):

    def test_iqtree_and_raxmlng(self):
        cases = {
            'iq_prot_LG_F_G4.iqtree': ('LG+F+G4', 'LG', 'empirical', 'gamma', 4),
            'iq_prot_LG.iqtree': ('LG', 'LG', 'model', 'none', None),
            'iq_prot_LG_I_G4.iqtree': ('LG+I+G4', 'LG', 'model', 'invariant+gamma', 4),
            'iq_prot_LG_R4.iqtree': ('LG+R4', 'LG', 'model', 'free-rate', 4),
            'iq_prot_Q_pfam_G4.iqtree': ('Q.pfam+G4', 'Q.pfam', 'model', 'gamma', 4),
            'rx_prot_LG_G4_FC.raxml.log': ('LG+FC+G4m', 'LG', 'empirical', 'gamma', 4),
            'rx_prot_LG_G4_FO.raxml.log': ('LG+FO+G4m', 'LG', 'ML', 'gamma', 4),
            'rx_dna_GTR_G4.raxml.log': ('GTR+FO+G4m', 'GTR', 'ML', 'gamma', 4),
        }
        for name, (model, matrix, freqs, rates, n) in cases.items():
            with self.subTest(fit=name):
                f = fitted.read_fit(fit(name))
                self.assertEqual((f['model'], f['matrix'], f['freqs'], f['rates'], f['n_cats']),
                                 (model, matrix, freqs, rates, n))
                self.assertFalse(f['median'])
        self.assertIsNone(fitted.read_fit(os.path.join(testsets.MINI, 'mini_tree.tree')))


class TestModelFromFits(unittest.TestCase):

    def setUp(self):
        fitted._warned.clear()                  # warnings are logged once per file

    def model(self, stats, aln=PROT_ALN, tree=None):
        return Reference(tree=tree or fit('iq_prot_LG_F_G4.treefile'), ref_aln=aln, stats=fit(stats)).phylo_model()

    def test_what_pplacer_drops(self):
        with self.assertLogs('timestamp', level='WARNING') as logs:
            pm = self.model('rx_prot_LG_G4_FC.raxml.log')
        self.assertEqual((pm['subs_model'], pm['ras_model'], pm['empirical_frequencies']), ('LG', 'gamma', False))
        self.assertIn("frequencies counted from the alignment (pplacer uses LG's own", logs.output[0])
        with self.assertLogs('timestamp', level='WARNING') as logs:
            self.model('rx_prot_LG_G4_FO.raxml.log')
        self.assertIn('ML-estimated frequencies (pplacer counts them from the reference alignment)', logs.output[0])

    def test_nothing_dropped(self):
        f = fitted.read_fit(fit('iq_prot_LG_F_G4.iqtree'))
        pm = self.model('iq_prot_LG_F_G4.iqtree')
        self.assertEqual(fitted.dropped(f, pm), [])

    def test_fits_pplacer_cannot_use(self):
        for name, why in (('iq_prot_LG.iqtree', 'no rate variation'), ('iq_prot_LG_I_G4.iqtree', 'invariant-sites'),
                          ('iq_prot_LG_R4.iqtree', 'free-rate'), ('iq_prot_Q_pfam_G4.iqtree', 'pplacer knows WAG'),
                          ('iq_dna_HKY_F_G4.iqtree', 'special case of GTR')):
            with self.subTest(fit=name), self.assertRaises(RefpkgMalformed) as cm:
                self.model(name, aln=DNA_ALN if '_dna_' in name else PROT_ALN)
            self.assertIn(why, str(cm.exception))
            self.assertIn('refit with', str(cm.exception))


class TestPlacedAsPplacer(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_fit_files_directly(self):
        """--tree/--ref_aln/--stats with the fit files: identical to pplacer on taxit's package."""
        for name, tree, aln, queries, tile in (
                ('iq_prot_LG_F_G4', 'iq_prot_LG_F_G4.treefile', PROT_ALN, QUERIES, 60),
                ('rx_dna_GTR_G4', 'rx_dna_GTR_G4.raxml.bestTree', DNA_ALN, GTR_Q, 20)):
            stats = fit(name + ('.iqtree' if name.startswith('iq') else '.raxml.log'))
            for tiled in (False, True):
                with self.subTest(fit=name, tiled=tiled):
                    ref = Reference(tree=fit(tree), ref_aln=aln, stats=stats)
                    p = new_placer(ref, core.read_fasta(queries), 2)
                    if tiled:
                        p = tiled_placer(p, tile)
                    results = [r for r in p.expand(p.run()) if r[1]]
                    out = os.path.join(self.tmp, 'out.jplace')
                    write_jplace(out, p.jplace_tree(), results, 'unittest', p.records([r[0] for r in results]))
                    r = compare_jplace(out, fit(name + '.pplacer.jplace'))
                    self.assertTrue(r['same_tree'])
                    self.assertTrue(r['same_records'])
                    self.assertEqual(r['identical'], r['genomes'])
