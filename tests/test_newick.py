import os
import unittest

from gplace.engine.newick import parse_newick, to_newick, quote_label
from gplace.fileio import read_text

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'gplace', 'tests', 'data')


class TestNewick(unittest.TestCase):

    def test_postorder_ids(self):
        t = parse_newick('((A:0.1,B:0.2)x:0.3,C:0.4)r;')
        # post-order: A=0, B=1, x=2, C=3, r=4 (these are the jplace edge numbers)
        self.assertEqual(t.label, ['A', 'B', 'x', 'C', 'r'])
        self.assertEqual(t.children[2], [0, 1])
        self.assertEqual(t.root, 4)
        self.assertIsNone(t.bl[4])

    def test_pplacer_format(self):
        t = parse_newick("((A:0.6845600000000001,B:0.2)'100.0:g__X;s__Y':0.3,C:1e-05)d__Bacteria;")
        self.assertEqual(to_newick(t),
                         "((A:0.68456{0},B:0.2{1})'100.0:g__X;s__Y':0.3{2},C:1e-05{3})d__Bacteria{4};")

    def test_quoting(self):
        self.assertEqual(quote_label('GB_GCA_000001.1'), 'GB_GCA_000001.1')
        self.assertEqual(quote_label('a b'), "'a b'")
        self.assertEqual(quote_label("it's"), "'it\\'s'")
        self.assertEqual(parse_newick("('it\\'s':1,b:2);").label[0], "it's")

    def test_root_branch_length_for_placement(self):
        t = parse_newick('(A:0.1,B:0.2);')
        self.assertEqual(t.branch_lengths_for_placement()[t.root], 0.0)

    def test_expected_tree_string(self):
        """The bundled pplacer output's tree string is reproduced byte for byte."""
        import json
        tree_txt = read_text(os.path.join(DATA, 'mini.refpkg', 'mini_tree.tree'))
        with open(os.path.join(DATA, 'expected_pplacer.jplace')) as fh:
            expected = json.load(fh)['tree']
        self.assertEqual(to_newick(parse_newick(tree_txt)), expected)


if __name__ == '__main__':
    unittest.main()
