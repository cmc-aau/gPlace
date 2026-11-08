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

"""Posterior probabilities (pplacer's -p, --unif-prior, --inform-prior, --prior-lower, --pp-rel-err)
against the output of the stock pplacer v1.1.alpha19 run with the same options
(tests/data/posterior), untiled and tiled; the -p option order; the priors' means."""

import os
import shutil
import tempfile
import unittest

from gplace.cli import get_main_parser
from gplace.engine import core, testsets
from gplace.engine.core import Search
from gplace.engine.factory import new_placer, tiled_placer
from gplace.engine.jplace import write_jplace, compare_jplace
from gplace.engine.newick import parse_newick
from gplace.engine.placer import tree_length, prior_mean_map

POST = os.path.join(testsets.TEST_DATA, 'posterior')
MINI = testsets.MINI
MINI_Q = os.path.join(POST, 'mini_q4.fasta')
PP = dict(calc_pp=True, keep_at_most=20, keep_factor=0.001, max_strikes=20)    # what -p sets
try:
    core.lib()
    ENGINES = True
except Exception:
    ENGINES = False


@unittest.skipUnless(ENGINES, 'engine not built')
class TestPosterior(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def assertPlacedAsPplacer(self, ref, queries, name, **search):
        for tiled in (False, True):
            with self.subTest(case=name, tiled=tiled):
                p = new_placer(ref, core.read_alignment(queries), 2, search=Search(**search))
                if tiled:
                    p = tiled_placer(p, 40)
                results = [r for r in p.expand(p.run()) if r[1]]
                out = os.path.join(self.tmp, 'out.jplace')
                write_jplace(out, p.jplace_tree(), results, 'unittest', p.records([r[0] for r in results]))
                r = compare_jplace(out, os.path.join(POST, f'{name}.pplacer.jplace'))
                self.assertTrue(r['same_tree'])
                self.assertTrue(r['same_fields'])
                self.assertEqual(r['identical'], r['genomes'])

    def test_default_prior(self):
        self.assertPlacedAsPplacer(MINI, MINI_Q, 'mini_p', **PP)

    def test_uniform_prior(self):
        self.assertPlacedAsPplacer(MINI, MINI_Q, 'mini_p_unif', uniform_prior=True, **PP)

    def test_informative_prior(self):
        self.assertPlacedAsPplacer(MINI, MINI_Q, 'mini_p_inform', informative_prior=True, prior_lower=0.1, **PP)

    def test_integration_not_converging(self):
        """--pp-rel-err 1e-6: inner integrations fail inside the outer one (NaN, as ocamlgsl's
        callback returns), the bounds are halved ten times and the likelihood is used instead."""
        self.assertPlacedAsPplacer(MINI, MINI_Q, 'mini_p_err1e-6', pp_rel_err=1e-6, **PP)

    def test_options_after_p(self):
        """--keep-at-most 3 -p --max-strikes 6: -p overrides the earlier option, the later one -p."""
        self.assertPlacedAsPplacer(MINI, MINI_Q, 'mini_keep3_p_strikes6',
                                   **dict(PP, max_strikes=6))

    def test_short_edges_and_zero_prior(self):
        """Edges shorter than 1e-5 (cut =~ 0: the pendant length only) and, with --inform-prior on a
        2e-6 leaf edge, a prior that is 0 at the ML pendant length (upper limit by halving)."""
        self.assertPlacedAsPplacer(os.path.join(POST, 'mini_short.refpkg'), os.path.join(POST, 'short_q6.fasta'),
                                   'mini_short_inform', informative_prior=True, **PP)

    def test_gamma_model(self):
        self.assertPlacedAsPplacer(os.path.join(testsets.MODELS, 'like_wag.refpkg'),
                                   os.path.join(POST, 'wag_g_q6.fasta'), 'wag_g_p', **PP)


class TestPosteriorOptions(unittest.TestCase):

    def parse(self, *args):
        return get_main_parser().parse_args(['place', '--refpkg', 'r', '--msa_file', 'q', '--out_file', 'o'] + list(args))

    def test_p_sets_search_where_it_appears(self):
        o = self.parse('--keep_at_most', '3', '-p', '--max_strikes', '6')
        self.assertEqual((o.calc_pp, o.keep_at_most, o.keep_factor, o.max_strikes), (True, 20, 0.001, 6))
        o = self.parse('--pp', '--keep-factor', '0.5')
        self.assertEqual((o.keep_at_most, o.keep_factor, o.max_strikes), (20, 0.5, 20))
        o = self.parse()
        self.assertEqual((o.calc_pp, o.keep_at_most, o.keep_factor, o.max_strikes), (False, 7, 0.01, 6))

    def test_prior_means(self):
        """Gtree.tree_length (root edge excluded) and Core.prior_mean_map on a small tree."""
        t = parse_newick('((A:1,B:2):0.5,(C:3,(D:4,E:5):1):2,F:0.25);')
        bl = t.branch_lengths_for_placement()
        self.assertEqual(tree_length(t.children, bl, t.root), 1 + 2 + 0.5 + 3 + 4 + 5 + 1 + 2 + 0.25)
        means = prior_mean_map(t.children, bl, t.root, 0.1)
        ab = t.children[t.children[t.root][0]]
        self.assertEqual(means[ab[0]], 0.1 + 1.)                        # a leaf: prior_lower + its length
        self.assertEqual(means[t.children[t.root][0]], ((0.1 + 0.5 + 1.) + (0.1 + 0.5 + 2.)) / 2.)
        self.assertNotIn(t.root, means)


if __name__ == '__main__':
    unittest.main()
