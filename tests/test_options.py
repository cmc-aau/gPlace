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

"""pplacer's search settings, --no-pre-mask, Stockholm and gzipped query files, reference files
without a reference package (-t/-r/-s and statistics files) and model refinement, against the
output of the stock pplacer v1.1.alpha19 run with the same options (tests/data/options)."""

import gzip
import json
import math
import os
import shutil
import tempfile
import unittest

from gplace.engine import core, testsets
from gplace.engine.core import Search
from gplace.engine.factory import new_placer, tiled_placer
from gplace.engine.jplace import write_jplace, compare_jplace
from gplace.engine.reference import Reference, ModelSettings
from gplace.engine.stats import parse_stats, StatsParsingError
from gplace.exceptions import GPlaceException

OPT = os.path.join(testsets.TEST_DATA, 'options')
MINI = testsets.MINI
MODELS = testsets.MODELS
QUERIES = os.path.join(testsets.TEST_DATA, 'queries.fasta')
GTR_Q = os.path.join(MODELS, 'like_gtr_queries.fasta')
WAG_Q = os.path.join(MODELS, 'like_wag_queries.fasta')
HIV = dict(tree=os.path.join(OPT, 'hiv.phy_phyml_tree.txt'), ref_aln=os.path.join(OPT, 'hiv.fasta'))
ACTIN = dict(tree=os.path.join(OPT, 'actin.phy_phyml_tree.txt'), ref_aln=os.path.join(OPT, 'actin.fasta'))


def expected(name):
    return os.path.join(OPT, f'{name}.pplacer.jplace')


class _Placing(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def assertPlacedAsPplacer(self, ref, queries, name, tile_edges=None, **kw):
        """Place untiled (and tiled with tile_edges) and compare with pplacer's jplace `name`."""
        for tiled in ((False, True) if tile_edges else (False,)):
            with self.subTest(case=name, tiled=tiled):
                p = new_placer(ref, queries, 2, **kw)
                if tiled:
                    p = tiled_placer(p, tile_edges)
                results = [r for r in p.expand(p.run()) if r[1]]     # unplaced: omitted, as pplacer
                out = os.path.join(self.tmp, 'out.jplace')
                write_jplace(out, p.jplace_tree(), results, 'unittest', p.records([r[0] for r in results]))
                r = compare_jplace(out, expected(name))
                self.assertTrue(r['same_tree'])
                self.assertTrue(r['same_records'])
                self.assertEqual(r['identical'], r['genomes'])


class TestSearchSettings(_Placing):

    def test_defaults_and_checks(self):
        s = Search()
        self.assertEqual((s.ml_tolerance, s.start_pend, s.max_pend, s.max_strikes, s.strike_box, s.max_pitches,
                          s.keep_at_most, s.keep_factor, s.no_pre_mask), (0.01, 0.1, 2.0, 6, 3.0, 40, 7, 0.01, False))
        self.assertEqual(s.changed(), [])
        self.assertEqual(Search(max_strikes=0).t_max_pitches, math.inf)          # 0: every edge
        self.assertEqual(Search(max_strikes=3).changed(), [('max_strikes', 3)])
        for kw, msg in ((dict(start_pend=0), 'strictly positive'), (dict(start_pend=-1), 'strictly positive'),
                        (dict(start_pend=2), 'strictly less than maximum'),
                        (dict(start_pend=0.5, max_pend=0.4), 'strictly less than maximum')):
            with self.assertRaises(GPlaceException) as cm:
                Search(**kw)
            self.assertIn(msg, str(cm.exception))

    def test_identical_to_pplacer(self):
        mini_q = core.read_fasta(QUERIES)
        for name, ref, queries, kw, tile in (
                ('mini.combo', MINI, mini_q, dict(ml_tolerance=0.001, start_pend=0.05, max_pend=5, strike_box=6,
                                                  max_pitches=80), 60),
                ('mini.p_like', MINI, mini_q, dict(max_strikes=20, keep_at_most=20, keep_factor=0.001), None),
                ('mini.npm_dups', MINI, core.read_fasta(os.path.join(testsets.PARITY, 'mini_queries_dups.fasta')),
                 dict(no_pre_mask=True), 60),
                ('gtr_g.ms0', os.path.join(MODELS, 'like_gtr.refpkg'), core.read_fasta(GTR_Q), dict(max_strikes=0), 20),
                ('wag_g.kam1', os.path.join(MODELS, 'like_wag.refpkg'), core.read_fasta(WAG_Q), dict(keep_at_most=1),
                 None),
                ('dna_cat.kf0', os.path.join(MODELS, 'dna_cat.refpkg'), core.read_fasta(GTR_Q), dict(keep_factor=0),
                 20)):
            self.assertPlacedAsPplacer(ref, queries, name, tile, search=Search(**kw))

    def test_nothing_kept(self):
        """keep_factor > 1 keeps no placement: pplacer omits those queries."""
        p = new_placer(MINI, core.read_fasta(QUERIES)[:3], 2, search=Search(keep_factor=2))
        self.assertTrue(all(not pl for _, pl, _ in p.expand(p.run())))

    def test_no_pre_mask_needs_a_residue(self):
        q = core.read_fasta(QUERIES)[:2]
        q.append(('EMPTY', '-' * len(q[0][1])))
        with self.assertRaises(GPlaceException) as cm:
            new_placer(MINI, q, 1, search=Search(no_pre_mask=True))
        self.assertIn("sequence 'EMPTY' has no informative sites.", str(cm.exception))


class TestAlignmentFiles(_Placing):

    def test_stockholm_reader(self):
        path = os.path.join(self.tmp, 'a.sto')
        open(path, 'w').write('# STOCKHOLM 1.0\n#=GF ID x\n\nb  AC.-\na  ac--  \n#=GC SS_cons ....\n\n'
                              'b  GT\na  gt\n//\n\n')
        self.assertEqual(core.read_stockholm(path), [('b', 'AC--GT'), ('a', 'ac--gt')])   # descending names
        self.assertEqual(core.read_alignment(path), core.read_stockholm(path))
        for text, msg in (('# STOCKHOLM 1.0\na AC\nb A\n//\n', "doesn't match the previous sequence lengths"),
                          ('# STOCKHOLM 1.0\na AC\n', "didn't reach footer by EOF"),
                          ('a AC\n# STOCKHOLM 1.0\n//\n', 'something before header'),
                          ('# STOCKHOLM 1.0\na AC\n//\nb AC\n', 'something after footer'),
                          ('# STOCKHOLM 1.0\na AC x\n//\n', 'line 2')):
            open(path, 'w').write(text)
            with self.assertRaises(GPlaceException) as cm:
                core.read_stockholm(path)
            self.assertIn(msg, str(cm.exception))

    def test_fasta_reader(self):
        path = os.path.join(self.tmp, 'a.fasta.gz')
        with gzip.open(path, 'wt') as fh:
            fh.write(';comment\n>a desc\nAC\n;more\nGT\n>b\nAA\n')
        self.assertEqual(core.read_alignment(path), [('a', 'ACGT'), ('b', 'AA')])
        open(os.path.join(self.tmp, 'b.fa'), 'w').write('AC\n>a\nAC\n')
        with self.assertRaises(GPlaceException):
            core.read_fasta(os.path.join(self.tmp, 'b.fa'))

    def test_stockholm_queries_identical(self):
        """The bundled queries as gzipped interleaved Stockholm (lowercase, '.', markup) place as
        the FASTA file does in pplacer."""
        q = core.read_fasta(QUERIES)
        lines = ['# STOCKHOLM 1.0', '#=GF ID q', '']
        for b in range(0, len(q[0][1]), 1000):
            lines += [f'{n}   {s[b:b + 1000].lower().replace("-", ".")}' for n, s in q] + ['']
        path = os.path.join(self.tmp, 'q.sto.gz')
        with gzip.open(path, 'wt') as fh:
            fh.write('\n'.join(lines + ['//']) + '\n')
        sto = core.read_alignment(path)
        self.assertEqual(sorted(sto), sorted((n, s.lower()) for n, s in q))
        p = new_placer(MINI, sto, 2)
        results = p.expand(p.run())
        out = os.path.join(self.tmp, 'out.jplace')
        write_jplace(out, p.jplace_tree(), results, 'unittest', p.records([r[0] for r in results]))
        r = compare_jplace(out, os.path.join(testsets.TEST_DATA, 'expected_pplacer.jplace'))
        self.assertEqual(r['identical'], r['genomes'])

    def test_refpkg_aln_sto(self):
        """A reference package with aln_sto (Stockholm) instead of aln_fasta, empirical
        frequencies (which depend on the reference order: pplacer's descending names)."""
        pkg = os.path.join(self.tmp, 'sto.refpkg')
        os.makedirs(pkg)
        files = json.load(open(os.path.join(MINI, 'CONTENTS.json')))['files']
        ref = core.read_fasta(os.path.join(MINI, files['aln_fasta']))
        with open(os.path.join(pkg, 'ref.sto'), 'w') as fh:
            fh.write('# STOCKHOLM 1.0\n' + ''.join(f'{n} {s}\n' for n, s in ref) + '//\n')
        shutil.copy(os.path.join(MINI, files['tree']), os.path.join(pkg, 'tree.nwk'))
        shutil.copy(os.path.join(MODELS, 'mini_WAG_emp.phylo_model.json'), os.path.join(pkg, 'pm.json'))
        json.dump({'files': {'aln_sto': 'ref.sto', 'tree': 'tree.nwk', 'phylo_model': 'pm.json'},
                   'metadata': {'format_version': '1.1'}}, open(os.path.join(pkg, 'CONTENTS.json'), 'w'))
        self.assertPlacedAsPplacer(pkg, core.read_fasta(QUERIES), 'refpkg_aln_sto', 60)


class TestReferenceFiles(_Placing):
    """pplacer -t/-r/-s: the model of a statistics file (always gamma)."""

    def test_stats_files(self):
        self.assertEqual(parse_stats(os.path.join(OPT, 'hiv.phy_phyml_stats.txt')),
                         ('GTR', 4, 0.771, [1.17728, 3.08642, 0.60698, 0.92672, 4.74366, 1.0]))
        self.assertEqual(parse_stats(os.path.join(OPT, 'nogamma_stats.txt'), 3, 2.0)[1:3], (3, 2.0))
        self.assertEqual(parse_stats(os.path.join(OPT, 'RAxML_info.gtr_fe'))[:3], ('GTR', 4, 0.771))
        self.assertEqual(parse_stats(os.path.join(OPT, 'RAxML_info.wag723')), ('WAG', 4, 0.289, None))
        self.assertEqual(parse_stats(os.path.join(OPT, 'actin_model.json')), ('WAG', 4, 0.289, None))
        for f, msg in (('RAxML8_info.gtr', "couldn't find alpha line"), ('phyml33_stats.txt', "header didn't match")):
            with self.assertRaises(StatsParsingError) as cm:
                parse_stats(os.path.join(OPT, f))
            self.assertIn(msg, str(cm.exception))

    def test_identical_to_pplacer(self):
        gtr_q, wag_q = core.read_fasta(GTR_Q), core.read_fasta(WAG_Q)
        for name, files, queries, settings, tile in (
                ('phyml_gtr', dict(HIV, stats=os.path.join(OPT, 'hiv.phy_phyml_stats.txt')), gtr_q, None, 20),
                ('phyml_nogamma_g', dict(HIV, stats=os.path.join(OPT, 'nogamma_stats.txt')), gtr_q,
                 ModelSettings(gamma_cats=3, gamma_alpha=2.0), None),
                ('raxml_gtr_fe', dict(HIV, stats=os.path.join(OPT, 'RAxML_info.gtr_fe')), gtr_q, None, None),
                ('raxml_wag', dict(ACTIN, stats=os.path.join(OPT, 'RAxML_info.wag723')), wag_q, None, 20),
                ('json_wag_modelfreqs', dict(ACTIN, stats=os.path.join(OPT, 'actin_model.json')), wag_q,
                 ModelSettings(model_freqs=True), None)):
            self.assertPlacedAsPplacer(Reference(settings=settings, **files), queries, name, tile)

    def test_references_in_the_query_file(self):
        hiv = core.read_fasta(HIV['ref_aln'])
        ref = Reference(tree=HIV['tree'], stats=os.path.join(OPT, 'hiv.phy_phyml_stats.txt'))
        self.assertPlacedAsPplacer(ref, core.read_fasta(GTR_Q) + hiv, 'refs_in_query_no_r')
        with self.assertRaises(GPlaceException) as cm:
            new_placer(ref, core.read_fasta(GTR_Q), 1)
        self.assertIn('Please specify a reference alignment', str(cm.exception))


class TestRefinement(_Placing):
    """Gcat_model.refine: new site categories, when the query file holds the reference sequences
    or with --always-refine (gamma models: no change)."""

    def test_identical_to_pplacer(self):
        files = json.load(open(os.path.join(MINI, 'CONTENTS.json')))['files']
        mini_ref = core.read_fasta(os.path.join(MINI, files['aln_fasta']))
        mini_q, gtr_q = core.read_fasta(QUERIES), core.read_fasta(GTR_Q)
        self.assertPlacedAsPplacer(MINI, mini_q + mini_ref, 'refine_mini_refs', 60)
        self.assertPlacedAsPplacer(MINI, mini_q, 'always_refine_mini', None, always_refine=True)
        self.assertPlacedAsPplacer(os.path.join(MODELS, 'dna_cat.refpkg'), gtr_q, 'always_refine_dna_cat', 20,
                                   always_refine=True)
        p = new_placer(os.path.join(MODELS, 'like_gtr.refpkg'), gtr_q[:2], 1, always_refine=True)
        self.assertFalse(p.refined)                                               # Gmix: no-op


if __name__ == '__main__':
    unittest.main()
