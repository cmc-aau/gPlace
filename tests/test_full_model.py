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

"""place --full_model: IQ-TREE's complete fitted model (tests/data/fits/full_*): a matrix pplacer lacks
(Q.pfam+G4), weighted rate categories (HKY+F+I+R2), mixtures (EX2, LG4X) and median gamma rates (LG+G4). The fits are read as IQ-TREE
fitted them, and the placement likelihoods equal IQ-TREE's own (full_model_expected.json). RAxML-NG's
complete model (tests/data/fits/rx_full_*): read as RAxML-NG reads its .raxml.bestModel; the placement
likelihoods equal RAxML-NG's own (rx_full_expected.json)."""

import os
import shutil
import tempfile
import unittest

from gplace.engine import core, iqtree, raxmlng, testsets
from gplace.engine.factory import new_placer
from gplace.engine.jplace import write_jplace
from gplace.engine.reference import Reference
from gplace.fileio import read_json

FITS = os.path.join(testsets.TEST_DATA, 'fits')
PROT = (os.path.join(testsets.MINI, 'mini_msa.fa'), os.path.join(testsets.TEST_DATA, 'queries.fasta'))
DNA = (os.path.join(testsets.MODELS, 'like_gtr.refpkg', 'aln.fasta'), os.path.join(testsets.MODELS, 'like_gtr_queries.fasta'))
CASES = {'prot_Q_pfam_G4': PROT, 'dna_HKY_F_I_R2': DNA, 'mix_EX2': PROT, 'mix_LG4X': PROT, 'prot_LG_G4_median': PROT}
ENGINE = os.path.exists(core.engine_path())


def fit_files(tag):
    return os.path.join(FITS, f'full_{tag}.ckp.gz'), os.path.join(FITS, f'full_{tag}.iqtree')


@unittest.skipUnless(ENGINE, 'engine not built')
class TestFullModel(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_read_fits(self):
        f = iqtree.read_fit(*fit_files('prot_Q_pfam_G4'))
        self.assertEqual((f['matrix'], f['rates'], f['n_cats'], f['alphabet']), ('Q.pfam', 'gamma', 4, core.PROTEIN))
        self.assertAlmostEqual(sum(f['freqs']), 1.0, places=15)
        f = iqtree.read_fit(*fit_files('dna_HKY_F_I_R2'))
        self.assertEqual((f['rates'], f['n_cats'], f['categories']), ('invariant+free-rate', 2, 3))
        self.assertEqual(f['R'][0][2], f['R'][1][3])                      # HKY: both transitions share kappa
        f = iqtree.read_fit(*fit_files('mix_EX2'))
        self.assertEqual(([c['matrix'] for c in f['components']], f['paired']), (['BurEX2', 'ExpEX2'], False))
        f = iqtree.read_fit(*fit_files('mix_LG4X'))
        self.assertEqual((len(f['components']), f['paired'], f['categories']), (4, True, 4))
        f = iqtree.read_fit(*fit_files('prot_LG_G4_median'))
        self.assertEqual((f['median'], f['n_cats']), (True, 4))                # --gamma-median
        self.assertFalse(iqtree.read_fit(*fit_files('prot_Q_pfam_G4'))['median'])

    def test_gamma_rates(self):
        """IQ-TREE's discrete gamma (Yang's algorithm, RateGamma::computeRates): bit-identical to its C++ code."""
        from gplace.engine import yang_gamma
        for n, median, p, alpha, rates in read_json(os.path.join(FITS, 'full_model_expected.json'))['iqtree_gamma']:
            with self.subTest(n=n, median=median, p_invar=p, alpha=alpha):
                self.assertEqual(yang_gamma.iqtree_rates(n, alpha, median, 1.0 / (1.0 - p) if p else 1.0), rates)

    def test_not_supported_yet(self):
        """Fits gPlace cannot place with yet stop with a message instead of being read in part."""
        from gplace.exceptions import RefpkgMalformed
        with self.assertRaises(RefpkgMalformed) as cm:
            iqtree.read_fit(os.path.join(FITS, 'full_partitioned_dna.ckp.gz'), os.path.join(FITS, 'full_partitioned_dna.iqtree'))
        self.assertIn('partitioned models', str(cm.exception))

    def test_rate_categories_average_one(self):
        """Weights sum to 1 and the mean rate is 1, as in IQ-TREE (+I: a rate-0 category)."""
        for tag in CASES:
            with self.subTest(fit=tag):
                f = iqtree.read_fit(*fit_files(tag))
                if f.get('components'):
                    rates, weights, _ = iqtree.mixture_categories(f)
                else:
                    rates, weights = iqtree.categories(f)
                self.assertAlmostEqual(sum(weights), 1.0, places=9)
                self.assertAlmostEqual(sum(w * r for w, r in zip(weights, rates)), 1.0, places=9)

    def test_placements_match_iqtree(self):
        expected = read_json(os.path.join(FITS, 'full_model_expected.json'))['fits']
        for tag, (aln, queries) in CASES.items():
            with self.subTest(fit=tag):
                ref = Reference(tree=os.path.join(FITS, f'full_{tag}.treefile'), ref_aln=aln,
                                stats=os.path.join(FITS, f'full_{tag}.iqtree'), full_model='auto')
                p = new_placer(ref, core.read_fasta(queries), 2)     # every query: pplacer's pre-mask uses them all
                results = [r for r in p.expand(p.run()) if r[1]]
                out = os.path.join(self.tmp, f'{tag}.jplace')
                write_jplace(out, p.jplace_tree(), results, 'unittest', p.records([r[0] for r in results]))
                J = read_json(out)
                f = J['fields']
                best = {(rec['nm'][0][0] if 'nm' in rec else rec['n'][0]): rec['p'][0] for rec in J['placements']}
                for name, e in expected[tag].items():
                    edge, like = best[name][f.index('edge_num')], best[name][f.index('likelihood')]
                    self.assertEqual(edge, e['edge_num'])
                    self.assertAlmostEqual(like, e['likelihood'], delta=1e-6 * max(1.0, abs(like)) + 1e-6)
                    self.assertAlmostEqual(like, e['iqtree_likelihood'], delta=1e-5)


RX_CASES = {'dna_HKY_FC_I_G4': DNA, 'dna_GTR_I_R3': DNA, 'prot_LG4X': PROT, 'prot_LG_G4a': PROT, 'prot_WAG_FC_IC_G4': PROT}


def rx_files(tag):
    return os.path.join(FITS, f'rx_full_{tag}.raxml.bestModel'), os.path.join(FITS, f'rx_full_{tag}.raxml.log')


@unittest.skipUnless(ENGINE, 'engine not built')
class TestFullModelRAxMLNG(unittest.TestCase):
    """place --full_model with RAxML-NG 2.0.3 fits (tests/data/fits/rx_full_*): the .raxml.bestModel read as
    RAxML-NG reads it; placement likelihoods equal RAxML-NG's own (rx_full_expected.json)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def model_file(self, text):
        path = os.path.join(self.tmp, 'run.raxml.bestModel')
        with open(path, 'w') as fh:
            fh.write(text + '\n')
        return path

    def test_read_fits(self):
        f = raxmlng.read_fit(*rx_files('dna_HKY_FC_I_G4'))
        self.assertEqual((f['alphabet'], f['freqs_mode'], f['p_invar'], f['rates'], f['categories']),
                         (core.DNA, 'empirical', 0.195313, 'invariant+gamma', 5))
        R = f['matrices'][0][0]
        self.assertEqual((R[0][2], R[1][3]), (4.166402, 4.166402))     # HKY: A-G and C-T share kappa
        self.assertEqual((R[0][1], R[0][3], R[1][2], R[2][3]), (1.0, 1.0, 1.0, 1.0))
        f = raxmlng.read_fit(*rx_files('prot_LG4X'))
        self.assertEqual(([c['matrix'] for c in f['components']], f['paired'], f['categories']),
                         (['LG4X1', 'LG4X2', 'LG4X3', 'LG4X4'], True, 4))
        weights, rates = f['free']
        self.assertAlmostEqual(sum(weights), 1.0, places=15)
        self.assertAlmostEqual(sum(w * r for w, r in zip(weights, rates)), 1.0, places=15)
        f = raxmlng.read_fit(*rx_files('prot_LG_G4a'))
        self.assertEqual((f['median'], f['n_cats']), (True, 4))
        f = raxmlng.read_fit(self.model_file('DNA010020{1/3.5/7}+FU{0.3/0.2/0.2/0.3}+G4m{0.5}, noname = 1-216'))
        self.assertEqual((f['matrices'][0][0][0][2], f['matrices'][0][0][1][3], f['matrices'][0][0][0][1]), (3.5, 7.0, 1.0))

    def test_gamma_rates(self):
        """RAxML-NG's discrete gamma (Yang's algorithm, coraxlib): bit-identical to its C code."""
        for n, median, alpha, rates in read_json(os.path.join(FITS, 'rx_full_expected.json'))['gamma']:
            with self.subTest(n=n, median=median, alpha=alpha):
                self.assertEqual(raxmlng.gamma_rates(n, alpha, median), rates)

    def test_counted_frequencies(self):
        """+FC / +IC: counted in the reference alignment as RAxML-NG counts them (its log, 6 decimals)."""
        f = raxmlng.read_fit(*rx_files('prot_WAG_FC_IC_G4'))
        ref = [(n, s.upper()) for n, s in core.read_fasta(PROT[0])]
        freqs, pinv = raxmlng.alignment_stats(ref, core.PROTEIN)
        self.assertLess(max(abs(a - b) for a, b in zip(freqs, f['log_freqs'])), 5.0001e-7)
        self.assertLess(abs(pinv - f['log_pinv']), 5.0001e-7)

    def test_not_supported_yet(self):
        """Fits gPlace cannot place with yet stop with a message instead of being read in part."""
        from gplace.exceptions import RefpkgMalformed
        cases = {'GTR{1/2/1/1/2/1}+FU{0.25/0.25/0.25/0.25}+G4m{0.5}, p1 = 1-100\n'
                 'GTR{1/2/1/1/2/1}+FU{0.25/0.25/0.25/0.25}+G4m{0.5}, p2 = 101-216': 'partitioned models',
                 'GTR{1/2/1/1/2/1}+FU{0.25/0.25/0.25/0.25}+G4m{0.5}+ASC_LEWIS, noname = 1-216': '+ASC_LEWIS',
                 'BIN+G4m{0.5}, noname = 1-100': 'binary',
                 'LG4M+IU{0.1}+G4m{0.5}, noname = 1-5036': 'invariable sites (+I) in a mixture',
                 'GTR{1/2/1/1/2/1}+FU{0.25/0.25/0.25/0.25}+G4m{0.5}, noname = 1-100\\3': 'whole alignment'}
        for text, message in cases.items():
            with self.subTest(model=text):
                with self.assertRaises(RefpkgMalformed) as cm:
                    raxmlng.read_fit(self.model_file(text))
                self.assertIn(message, str(cm.exception))
        f = raxmlng.read_fit(self.model_file('LG+G4m{0.5}, noname = 1-100'))
        with self.assertRaises(RefpkgMalformed) as cm:
            raxmlng.complete(f, [(n, s.upper()) for n, s in core.read_fasta(PROT[0])])
        self.assertIn('covers columns 1-100', str(cm.exception))

    def test_placements_match_raxmlng(self):
        expected = read_json(os.path.join(FITS, 'rx_full_expected.json'))['fits']
        for tag, (aln, queries) in RX_CASES.items():
            with self.subTest(fit=tag):
                ref = Reference(tree=os.path.join(FITS, f'rx_full_{tag}.raxml.bestTree'), ref_aln=aln,
                                stats=os.path.join(FITS, f'rx_full_{tag}.raxml.log'), full_model='auto')
                p = new_placer(ref, core.read_fasta(queries), 2)
                results = [r for r in p.expand(p.run()) if r[1]]
                out = os.path.join(self.tmp, f'{tag}.jplace')
                write_jplace(out, p.jplace_tree(), results, 'unittest', p.records([r[0] for r in results]))
                J = read_json(out)
                f = J['fields']
                best = {(rec['nm'][0][0] if 'nm' in rec else rec['n'][0]): rec['p'][0] for rec in J['placements']}
                for name, e in expected[tag].items():
                    edge, like = best[name][f.index('edge_num')], best[name][f.index('likelihood')]
                    self.assertEqual(edge, e['edge_num'])
                    self.assertAlmostEqual(like, e['likelihood'], delta=1e-6 * max(1.0, abs(like)) + 1e-6)
                    self.assertAlmostEqual(like, e['raxmlng_likelihood'], delta=1e-6)
