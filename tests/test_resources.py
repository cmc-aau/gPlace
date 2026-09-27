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
