import math
import os
import tempfile
import unittest

from gplace.engine import core, models, testsets
from gplace.exceptions import RefpkgMalformed
from gplace.fileio import read_json, write_json

ENGINES = all(os.path.exists(core.engine_path(ns)) for ns in (20, 4))


class TestModelSelection(unittest.TestCase):
    """pplacer's rules for choosing the model of a refpkg (no engine needed)."""

    def test_alphabet_and_width(self):
        self.assertIs(models.alphabet_of({'subs_model': 'LG'}), core.PROTEIN)
        self.assertIs(models.alphabet_of({'subs_model': 'GTR'}), core.DNA)
        self.assertEqual(models.width_of({'subs_model': 'WAG', 'ras_model': 'Price-CAT'}), 20)
        self.assertEqual(models.width_of({'subs_model': 'GTR', 'ras_model': 'Price-CAT'}), 4)
        self.assertEqual(models.width_of({'subs_model': 'JTT', 'ras_model': 'gamma',
                                          'gamma': {'alpha': 0.7, 'n_cats': 8}}), 160)
        with self.assertRaises(RefpkgMalformed):
            models.alphabet_of({'subs_model': 'HKY'})

    def test_dna_alphabet_is_pplacer_nuc_map(self):
        m = core.DNA.char_map
        self.assertEqual(list(m['R']), [1, 0, 1, 0])          # G or A
        self.assertEqual(list(m['B']), [0, 1, 1, 1])          # G or T or C
        for c in 'N-?X':
            self.assertEqual(list(m[c]), [1, 1, 1, 1])
        self.assertNotIn('U', m)                              # pplacer rejects U as well
        self.assertEqual(core.DNA.gap_code, len(core.DNA.codes))

    def test_gtr_matrix(self):
        R = models.gtr_exchangeabilities({'ac': 1, 'ag': 2, 'at': 3, 'cg': 4, 'ct': 5, 'gt': 6})
        self.assertEqual(R, [[0, 1, 2, 3], [1, 0, 4, 5], [2, 4, 0, 6], [3, 5, 6, 0]])
        with self.assertRaises(RefpkgMalformed):
            models.gtr_exchangeabilities({'ac': 1})

    def test_protein_tables(self):
        for name in models.PROTEIN_MODELS:
            R, pi = core.load_protein_model(name)
            self.assertEqual((len(R), len(pi)), (20, 20))
            self.assertTrue(all(R[i][j] == R[j][i] for i in range(20) for j in range(20)))
            self.assertAlmostEqual(sum(pi), 1.0, places=5)

    def test_empirical_frequencies_are_sequential(self):
        """emper_freq: ambiguity codes add fractions, so the sum is order-dependent; check against
        a plain loop over the characters in alignment order."""
        aln = [('a', 'ACGTX-N'), ('b', 'RYMK?SWHBVDN'), ('c', 'AAAAX')]
        got = models.emper_freq(aln, core.DNA)
        tot = [0.0] * 4
        for _, seq in aln:
            for ch in seq:
                if ch in '-?':
                    continue
                v = list(core.DNA.char_map[ch]); n = 0.0
                for x in v:
                    n = n + x
                for i in range(4):
                    tot[i] = tot[i] + v[i] * (1.0 / n)
        n = 0.0
        for x in tot:
            n = n + x
        self.assertEqual(got, [x * (1.0 / n) for x in tot])
        with self.assertRaises(RefpkgMalformed):
            models.emper_freq([('a', 'ACGU')], core.DNA)


@unittest.skipUnless(ENGINES, 'engines not built (run: gplace build_engine)')
class TestModelsEngine(unittest.TestCase):

    def test_discrete_gamma(self):
        from gplace.engine.selftest import check_discrete_gamma
        self.assertEqual(check_discrete_gamma(), [])
        self.assertEqual(models.discrete_gamma(4, 0.05)[0], 0.0)       # pplacer's epsilon rule
        r = models.discrete_gamma(8, 0.7)
        self.assertAlmostEqual(sum(r) / 8, 1.0, places=12)

    def test_discrete_gamma_many_categories(self):
        """21 or more categories: GSL calls log1p (glibc 2.19's, built in), not an aborting trap."""
        for n, a in [(21, 2.0), (24, 0.6), (32, 0.3), (100, 1.0)]:
            r = models.discrete_gamma(n, a)
            self.assertEqual(len(r), n)
            self.assertAlmostEqual(sum(r) / n, 1.0, places=10)

    def test_discrete_gamma_gsl_error_stops_like_pplacer(self):
        """pplacer turns GSL errors into exceptions: alpha 1e6 stops it with
        Gsl_error.Gsl_exn(6, "inverse failed to converge"); gPlace stops with the same message."""
        from gplace.exceptions import RefpkgMalformed
        with self.assertRaises(RefpkgMalformed) as cm:
            models.discrete_gamma(4, 1e6)
        self.assertIn('Gsl_error.Gsl_exn(6, "inverse failed to converge")', str(cm.exception))
        models.discrete_gamma(4, 0.5)                                    # the next model is unaffected

    def test_dna_and_gamma_kernels(self):
        from gplace.engine.selftest import check_kernels, check_gamma_kernels
        for ns in (20, 4):
            self.assertEqual(check_kernels(ns=ns)[1], [])
            self.assertEqual(check_gamma_kernels(ns=ns)[1], [])

    def test_model_test_sets_identical(self):
        """Every bundled model test set, untiled and tiled, identical to its expected output (the
        stock pplacer v1.1.alpha19). Sets whose untiled placement does not fit in the memory
        available (many gamma categories on a small machine) are checked tiled only."""
        from gplace.engine import resources, taxonomy
        from gplace.engine.factory import new_placer, tiled_placer
        from gplace.engine.jplace import write_jplace, compare_jplace
        with tempfile.TemporaryDirectory() as d:
            for ts in testsets.model_sets(d):
                for tiled in (False, True):
                    with self.subTest(test_set=ts.name, tiled=tiled):
                        p = new_placer(ts.refpkg, core.read_fasta(ts.queries), 2, ts.fix_gamma, ts.discard)
                        avail = resources.available_memory()
                        if not tiled and avail and p.memory_estimate(2) > avail:
                            self.skipTest(f'untiled needs {p.memory_estimate(2) / 2 ** 30:.1f} GiB, '
                                          f'{avail / 2 ** 30:.1f} GiB available (checked tiled)')
                        classify = (taxonomy.Classifier(ts.refpkg, p.tree, ts.mrca_class)
                                    if taxonomy.tax_equipped(ts.refpkg) else None)
                        if tiled:
                            p = tiled_placer(p, ts.tile_edges)
                        out = os.path.join(d, 'out.jplace')
                        results = p.expand(p.run())
                        write_jplace(out, p.jplace_tree(), results, 'unittest',
                                     p.records([r[0] for r in results]), classify)
                        r = compare_jplace(out, ts.expected)
                        self.assertTrue(r['same_tree'])
                        self.assertTrue(r['same_records'])
                        self.assertTrue(r['same_fields'])
                        self.assertEqual(r['identical'], r['genomes'])

    def test_gamma_fix_changes_only_the_last_digits(self):
        """At alpha 0.05 the blocked rescaling loses precision: the fix changes values slightly,
        never the edges."""
        from gplace.engine.jplace import compare_jplace
        ref = os.path.join(testsets.MODELS, 'mini_WAG_g0.05x4.pplacer.jplace')
        fix = os.path.join(testsets.MODELS, 'mini_WAG_g0.05x4.fix.jplace')
        r = compare_jplace(ref, fix)
        self.assertEqual(r['edge_list'], r['genomes'])
        self.assertLess(r['identical'], r['genomes'])
        self.assertLess(r['max_delta']['likelihood'], 1e-3)

    def test_unsupported_models_rejected(self):
        from gplace.engine.factory import new_placer
        with tempfile.TemporaryDirectory() as d:
            ts = testsets.model_sets(d)[0]
            contents = read_json(os.path.join(ts.refpkg, 'CONTENTS.json'))
            for bad in ({'subs_model': 'HKY', 'ras_model': 'Price-CAT'},
                        {'subs_model': 'WAG', 'ras_model': 'free'},
                        {'subs_model': 'GTR', 'ras_model': 'gamma', 'gamma': {'alpha': 1, 'n_cats': 4}}):
                pkg = os.path.join(d, 'bad.refpkg'); os.makedirs(pkg, exist_ok=True)
                write_json(os.path.join(pkg, 'pm.json'), bad)
                files = dict(contents['files'], phylo_model=os.path.join(pkg, 'pm.json'))
                write_json(os.path.join(pkg, 'CONTENTS.json'), {'files': files})
                with self.assertRaises(RefpkgMalformed):
                    new_placer(pkg, core.read_fasta(ts.queries)[:2], 1)

    def test_cache_is_deterministic_with_many_threads(self):
        """Building the likelihood cache with many threads on a small tree, where tasks finish
        while the first ones are still being submitted, gives the same vectors every time (a node
        must never be submitted twice)."""
        import hashlib
        from gplace.engine.factory import new_placer
        with tempfile.TemporaryDirectory() as d:
            ts = [t for t in testsets.model_sets(d) if 'CAT (FastTree)' in t.name][0]
            p = new_placer(ts.refpkg, core.read_fasta(ts.queries), 32)
            seen = set()
            for _ in range(25):
                p.cache()
                seen.add(hashlib.md5(p.darr.tobytes() + p.parr.tobytes()).hexdigest())
            self.assertEqual(len(seen), 1)


class TestNonFinite(unittest.TestCase):

    @unittest.skipUnless(ENGINES, 'engine not built')
    def test_negative_nan_written_as_pplacer(self):
        """pplacer (x86-64, glibc printf) writes the like-weight ratio of a -inf likelihood as -nan:
        the NaN from -inf - (-inf) has its sign bit set on x86-64. ARM64 makes a positive NaN, so
        gPlace sets the sign itself; the result must be -nan on every platform."""
        from gplace.engine.jplace import _f
        from gplace.engine.placer import weighted_placements
        for lls in ([-math.inf, -math.inf], [-100.0, -math.inf]):
            pl = weighted_placements(core.lib(), [(e, (ll, 0.1, 0.2)) for e, ll in enumerate(lls)])
            lwr = {p[1]: p[2] for p in pl}
            self.assertTrue(math.isnan(lwr[1]) and math.copysign(1.0, lwr[1]) < 0)
            self.assertEqual(_f(lwr[1]), '-nan')
        self.assertEqual(lwr[0], 1.0)                              # the finite one keeps weight 1
        self.assertEqual(_f(math.copysign(float('nan'), 1.0)), 'nan')
        self.assertEqual(_f(-math.inf), '-inf')

if __name__ == '__main__':
    unittest.main()
