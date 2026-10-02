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


if __name__ == '__main__':
    unittest.main()


@unittest.skipUnless(ENGINE, 'engine not built')
class TestEngineImports(unittest.TestCase):

    def test_only_exact_libm_functions_imported(self):
        """The wheel runs on any distribution's libm; only functions that are exact everywhere
        (sqrt, frexp, ldexp, fmod, round, ...) may come from it. exp/log/pow/sin/hypot are built in,
        also for GSL's own calls (build.py LIBM_CR); other libm functions GSL references are traps."""
        import shutil
        import subprocess
        nm = shutil.which('nm')
        if nm is None:
            self.skipTest('nm not available')
        imported = set()
        for ns in (20, 4):
            out = subprocess.run([nm, '-D', '--undefined-only', core.engine_path(ns)],
                                 capture_output=True, text=True).stdout
            imported |= {line.split()[-1].split('@')[0] for line in out.splitlines() if line.strip()}
        libm = {'exp', 'log', 'pow', 'hypot', 'expm1', 'log1p', 'log2', 'log10', 'exp2', 'cbrt',
                'sin', 'cos', 'tan', 'atan', 'atan2', 'lgamma', 'tgamma', 'erf', 'erfc', 'fma', 'sqrt',
                'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh', 'asin', 'acos', 'sincos',
                'frexp', 'ldexp', 'fabs', 'floor', 'ceil', 'round', 'trunc', 'fmod', 'modf'}
        allowed = {'sqrt', 'frexp', 'ldexp', 'fabs', 'floor', 'ceil', 'trunc', 'modf', 'fmod', 'round'}
        self.assertEqual(imported & libm - allowed, set())
