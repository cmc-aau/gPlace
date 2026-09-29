import os
import tempfile
import unittest
from unittest.mock import patch

from gplace.engine import resources


class TestResources(unittest.TestCase):

    def test_parse_cpus(self):
        self.assertEqual(resources.parse_cpus('auto'), 'auto')
        self.assertEqual(resources.parse_cpus('AUTO'), 'auto')
        self.assertEqual(resources.parse_cpus('12'), 12)
        with self.assertRaises(ValueError):
            resources.parse_cpus('0')

    def test_candidates_and_heuristic(self):
        self.assertEqual(resources.default_candidates(64), [8, 16, 24, 32, 40, 48, 64])
        self.assertEqual(resources.default_candidates(1), [1])
        self.assertEqual(resources.heuristic_cpus(64), 40)

    def test_memory_model_tree2(self):
        """GTDB-Tk R226 class-level tree 2, 22,061 queries, 40 workers: measured 43.8 GiB peak."""
        m = resources.estimate_memory(43130, 21564, 5036, 22061, 21566, 5036, 40) / 2 ** 30
        self.assertGreater(m, 43.8)
        self.assertLess(m, 43.8 * 1.10)

    def test_memory_grows_with_queries_and_workers(self):
        a = resources.estimate_memory(1000, 499, 5000, 100, 501, 5036, 4)
        self.assertLess(a, resources.estimate_memory(1000, 499, 5000, 1000, 501, 5036, 4))
        self.assertLess(a, resources.estimate_memory(1000, 499, 5000, 100, 501, 5036, 16))

    def test_tiled_memory_model_full_tree(self):
        """GTDB R226 full bacterial tree, 24 tiles of <= 13,889 edges (6,944 slots), 40 workers:
        measured 14.9 GiB (4,040 queries) and 19.9 GiB (127,343 queries)."""
        for q, measured in ((4040, 14.9), (127343, 19.9)):
            m = resources.estimate_memory_tiled(6944, 13889, 24, 5036, q, 136646, 5036, 40) / 2 ** 30
            self.assertGreater(m, measured)
            self.assertLess(m, measured * 1.10)

    def test_tiled_memory_grows_with_tiles_and_queries(self):
        a = resources.estimate_memory_tiled(1000, 2000, 10, 5000, 100, 1001, 5036, 4)
        self.assertLess(a, resources.estimate_memory_tiled(2000, 4000, 10, 5000, 100, 1001, 5036, 4))
        self.assertLess(a, resources.estimate_memory_tiled(1000, 2000, 10, 5000, 10000, 1001, 5036, 4))

    def test_parse_memory(self):
        G = 2 ** 30
        self.assertEqual(resources.parse_memory('64G'), 64 * G)
        self.assertEqual(resources.parse_memory('64gb'), 64 * G)
        self.assertEqual(resources.parse_memory('64GiB'), 64 * G)
        self.assertEqual(resources.parse_memory('64'), 64 * G)
        self.assertEqual(resources.parse_memory('1.5T'), int(1.5 * 1024 * G))
        self.assertEqual(resources.parse_memory('500M'), 500 * 2 ** 20)
        self.assertEqual(resources.parse_memory('4096K'), 4 * 2 ** 20)
        for bad in ('abc', '', '-5G', '0', '64X'):
            with self.assertRaises(ValueError):
                resources.parse_memory(bad)

    def test_available_memory(self):
        m = resources.available_memory()
        self.assertTrue(m is None or m > 0)

    def test_saved_tuning_roundtrip(self):
        with tempfile.TemporaryDirectory() as d, patch.dict(os.environ, {'GPLACE_CACHE_DIR': d}):
            self.assertIsNone(resources.load_tuned_cpus())
            resources.save_tuned_cpus(7, {'refpkg': 'x'})
            self.assertEqual(resources.load_tuned_cpus(), 7)

    def test_slurm_time(self):
        self.assertEqual(resources.slurm_time(3725), '01:02:05')
        self.assertEqual(resources.slurm_time(90000), '1-01:00:00')


if __name__ == '__main__':
    unittest.main()
