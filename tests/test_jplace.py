import os
import tempfile
import unittest

from gplace.engine.jplace import write_jplace, read_jplace, compare_jplace

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'gplace', 'tests', 'data')


class TestJplace(unittest.TestCase):

    def test_number_format_is_pplacer(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.jplace')
            write_jplace(path, '(A:1{0},B:1{1}){2};',
                         [('q1', [[0.0758007809160999, 11, 1.0, -13036489.620231, 0.0629997563923999]], [])],
                         'test')
            text = open(path).read()
            self.assertIn('[0.0758007809161, 11, 1, -13036489.6202, 0.0629997563924]', text)
            _, pl, fields = read_jplace(path)
            self.assertEqual(fields, ['distal_length', 'edge_num', 'like_weight_ratio', 'likelihood', 'pendant_length'])
            self.assertEqual(pl['q1'][0]['edge_num'], 11)

    def test_nonfinite_values_are_read(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.jplace')
            write_jplace(path, '(A:1{0},B:1{1}){2};', [('q1', [[0.1, 0, float('nan'), -float('inf'), 0.2]], [])], 't')
            _, pl, _ = read_jplace(path)
            self.assertTrue(pl['q1'][0]['likelihood'] == -float('inf'))

    def test_compare_identical(self):
        exp = os.path.join(DATA, 'expected_pplacer.jplace')
        r = compare_jplace(exp, exp)
        self.assertEqual(r['identical'], r['genomes'])
        self.assertEqual(r['genomes'], 20)


if __name__ == '__main__':
    unittest.main()
