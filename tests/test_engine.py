import os
import tempfile
import unittest

from gplace.engine import core

ENGINE = os.path.exists(core.engine_path())
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'gplace', 'tests', 'data')


@unittest.skipUnless(ENGINE, 'engine not built (run: gplace build_engine)')
class TestEngine(unittest.TestCase):

    def test_libm_matches_pplacer_binary(self):
        from gplace.engine.selftest import check_libm
        n, bad = check_libm()
        self.assertGreater(n, 0)
        self.assertEqual(bad, 0)

    def test_wag_eigen_matches_pplacer_binary(self):
        from gplace.engine.selftest import check_wag_eigen
        self.assertEqual(check_wag_eigen(), (True, True, True))

    def test_simd_kernels_match_scalar(self):
        from gplace.engine.selftest import check_kernels
        _, fails = check_kernels()
        self.assertEqual(fails, [])

    def test_mini_placement_identical_to_pplacer(self):
        from gplace.engine.placer import Placer
        from gplace.engine.jplace import write_jplace, compare_jplace
        placer = Placer(os.path.join(DATA, 'mini.refpkg'), core.read_fasta(os.path.join(DATA, 'queries.fasta')), 2)
        results = placer.run()
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'gplace.jplace')
            write_jplace(out, placer.jplace_tree(), results, 'unittest')
            r = compare_jplace(out, os.path.join(DATA, 'expected_pplacer.jplace'))
        self.assertTrue(r['same_tree'])
        self.assertEqual(r['identical'], 20)

    def test_ranking_without_a_score_table(self):
        """--ranked_edges: the untiled prefilter keeps only each query's first ranked edges, scoring
        the queries in batches; a query that needs more is ranked in full when it gets there. With
        one stored edge per query (every search goes past it) and batches of 3 queries, placements
        are still pplacer's; memory no longer grows with queries x edges. Default: the score table."""
        from gplace.engine import resources
        from gplace.engine.placer import Placer
        from gplace.engine.jplace import write_jplace, compare_jplace
        placer = Placer(os.path.join(DATA, 'mini.refpkg'), core.read_fasta(os.path.join(DATA, 'queries.fasta')), 2)
        placer.cache()
        placer.prefilter()
        self.assertIsNone(placer.rank)                                   # default: score table
        self.assertEqual(placer.scores.shape, (20, placer.nnodes - 1))
        placer.prefilter(n_cand=1, query_batch=3)
        self.assertEqual(placer.rank.shape, (20, 1))
        results = placer.expand(placer.place_all())
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'gplace.jplace')
            write_jplace(out, placer.jplace_tree(), results, 'unittest')
            r = compare_jplace(out, os.path.join(DATA, 'expected_pplacer.jplace'))
        self.assertEqual(r['identical'], 20)
        est = lambda q, k: resources.estimate_memory(273290, 136644, 5036, q, 136646, 5036, 40, n_cand=k)
        self.assertLess(est(1_000_000, 164) - est(100_000, 164), 30 * 2 ** 30)
        self.assertGreater(est(1_000_000, None) - est(100_000, None), 1500 * 2 ** 30)   # the table


if __name__ == '__main__':
    unittest.main()


@unittest.skipUnless(ENGINE, 'engine not built')
class TestEngineImports(unittest.TestCase):

    def test_only_exact_libm_functions_imported(self):
        """The wheel runs on any distribution's libm; only functions that are exact everywhere
        (sqrt, frexp, ldexp, fmod, round, ...) may come from it. exp/log/pow/sin/hypot/log1p are built in,
        also for GSL's own calls (build.py LIBM_CR); other libm functions GSL references are traps."""
        import shutil
        nm = shutil.which('nm')
        if nm is None:
            self.skipTest('nm not available')
        from gplace.config.common import ENGINE_VARIANTS
        from gplace.engine.build import undefined_symbols
        allowed_extra = {'avx2': set(), 'neon': set(), 'portable': {'fma'}}   # fma: correctly rounded everywhere
        found = 0
        for variant in ENGINE_VARIANTS:
            imported = set()
            for ns in (20, 4):
                path = core.engine_path(ns, variant)
                if not os.path.exists(path):
                    continue
                found += 1
                imported |= set(undefined_symbols(path))
            self._check(imported, allowed_extra[variant])
        self.assertGreater(found, 0)

    def _check(self, imported, allowed_extra):
        libm = {'exp', 'log', 'pow', 'hypot', 'expm1', 'log1p', 'log2', 'log10', 'exp2', 'cbrt',
                'sin', 'cos', 'tan', 'atan', 'atan2', 'lgamma', 'tgamma', 'erf', 'erfc', 'fma', 'sqrt',
                'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh', 'asin', 'acos', 'sincos',
                'frexp', 'ldexp', 'fabs', 'floor', 'ceil', 'round', 'trunc', 'fmod', 'modf',
                '__sincos_stret', '__sincospi_stret', '__exp10'}         # macOS forms of sincos and exp10
        allowed = {'sqrt', 'frexp', 'ldexp', 'fabs', 'floor', 'ceil', 'trunc', 'modf', 'fmod', 'round'}
        self.assertEqual(imported & libm - allowed - allowed_extra, set())


class TestGslDownload(unittest.TestCase):

    def test_next_mirror_when_one_fails(self):
        """An unreachable server and a file with the wrong SHA-256 are skipped; the first good copy
        is kept, and no partial file is left behind when every mirror fails."""
        import hashlib
        from gplace.engine.build import download_gsl
        from gplace.exceptions import GPlaceExit
        with tempfile.TemporaryDirectory() as d:
            good, bad = os.path.join(d, 'good.tar.gz'), os.path.join(d, 'bad.tar.gz')
            with open(good, 'wb') as fh:
                fh.write(b'gsl 1.16')
            with open(bad, 'wb') as fh:
                fh.write(b'an error page')
            sha = hashlib.sha256(b'gsl 1.16').hexdigest()
            out = os.path.join(d, 'gsl-1.16.tar.gz')
            urls = ('http://127.0.0.1:9/gsl-1.16.tar.gz', 'file://' + bad, 'file://' + good)
            self.assertEqual(download_gsl(out, urls=urls, sha256=sha, timeout=5, rounds=1), out)
            with open(out, 'rb') as fh:
                self.assertEqual(fh.read(), b'gsl 1.16')
            os.remove(out)
            with self.assertRaises(GPlaceExit) as ctx:
                download_gsl(out, urls=urls[:2], sha256=sha, timeout=5, rounds=2)
            self.assertIn('GPLACE_GSL_TARBALL', str(ctx.exception))
            self.assertEqual(str(ctx.exception).count('127.0.0.1'), 2)
            self.assertFalse(os.path.exists(out) or os.path.exists(out + '.part'))
