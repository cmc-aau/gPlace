import os
import tempfile
import unittest

from gplace.engine import core, testsets
from gplace.engine.newick import parse_newick
from gplace.engine.taxonomy import Taxonomy, mrca_map, edge_painting, read_seqinfo
from gplace.exceptions import AlignmentMalformed, RefpkgMalformed

ENGINE = os.path.exists(core.engine_path())
PARITY = testsets.PARITY


def write(d, name, text):
    path = os.path.join(d, name)
    with open(path, 'w') as fh:
        fh.write(text)
    return path


TAXTABLE = ('"tax_id","parent_id","rank","tax_name","root","phylum","genus"\n'
            '"root","root","root","root","root","",""\n'
            '"p1","root","phylum","P1","root","p1",""\n'
            '"g1","p1","genus","G1","root","p1","g1"\n'
            '"g2","p1","genus","G2","root","p1","g2"\n'
            '"g3","root","genus","G3","root","","g3"\n')


class TestTaxonomy(unittest.TestCase):

    def test_taxtable_lineage_and_mrca(self):
        with tempfile.TemporaryDirectory() as d:
            t = Taxonomy(write(d, 't.csv', TAXTABLE))
        self.assertEqual(t.rank_names, ['root', 'phylum', 'genus'])
        self.assertEqual(t.lineage('g1'), ['root', 'p1', 'g1'])
        self.assertEqual(t.lineage('g3'), ['root', 'g3'])          # rank gap: parent from the lineage
        self.assertEqual(t.mrca('g1', 'g2'), 'p1')
        self.assertEqual(t.mrca('g1', 'g3'), 'root')
        self.assertEqual(t.mrca('g1', None), 'g1')

    def test_broken_taxtable(self):
        bad = TAXTABLE + '"g1","root","genus","G1","root","","g1"\n'      # g1 under two parents
        with tempfile.TemporaryDirectory() as d, self.assertRaises(RefpkgMalformed):
            Taxonomy(write(d, 't.csv', bad))

    def test_painting_and_mrca_on_a_small_tree(self):
        """((a:1,b:1)x:1,(c:1,d:1)y:1)r;  a,b in g1; c in g2; d in g3."""
        tree = parse_newick('((a:1,b:1):1,(c:1,d:1):1);')
        with tempfile.TemporaryDirectory() as d:
            tax = Taxonomy(write(d, 't.csv', TAXTABLE))
            si = read_seqinfo(write(d, 's.csv', 'seqname,tax_id\na,g1\nb,g1\nc,g2\nd,g3\n'))
        node = {tree.label[v]: v for v in tree.leaves()}
        par = {c: v for v, ch in enumerate(tree.children) for c in ch}
        x, y = par[node['a']], par[node['c']]
        mrca = mrca_map(tax, si, tree.label, tree.children, tree.root)
        self.assertEqual((mrca[x], mrca[y], mrca[tree.root]), ('g1', 'root', 'root'))
        paint, ranks = edge_painting(tax, si, tree.label, tree.children, tree.root)
        # an edge gets a taxon it *cuts* (present on both sides): inside g1 -> g1; the stem of the g1
        # clade and the edge to the lone g2 leaf lie inside p1 only
        self.assertEqual(paint[node['a']], 'g1')
        self.assertEqual(paint[x], 'p1')
        self.assertEqual(paint[node['c']], 'p1')
        self.assertEqual(ranks, {0, 1, 2})


@unittest.skipUnless(ENGINE, 'engine not built')
class TestDuplicatesAndOverlap(unittest.TestCase):

    def placer(self, queries, **kw):
        from gplace.engine.factory import new_placer
        return new_placer(testsets.MINI, core.read_fasta(queries), 2, **kw)

    def test_duplicates_share_a_record_in_input_order(self):
        p = self.placer(os.path.join(PARITY, 'mini_queries_dups.fasta'))
        self.assertEqual((len(p.all_names), len(p.query_names)), (24, 21))
        recs = p.records(p.all_names)
        self.assertEqual(len(recs), 21)
        self.assertIn(['GEMOTU_05958', 'GEMOTU_05958_DUP_exact', 'GEMOTU_05958_DUP_exact2'], recs)
        self.assertIn(['GEMOTU_17643', 'GEMOTU_17643_DUP_maskedout'], recs)   # differs only where masked
        self.assertNotIn('GOMCOTU_14244_GAPVAR', [n for r in recs if len(r) > 1 for n in r])  # '-' vs '?'
        # per output file: duplicates in different files are different records (place_batch)
        a, b = ['GEMOTU_05958', 'GEMOTU_17643'], ['GEMOTU_05958_DUP_exact', 'GEMOTU_17643_DUP_maskedout']
        self.assertEqual(p.records(a), [[n] for n in a])
        self.assertEqual(p.records(b), [[n] for n in b])

    def test_nonoverlapping_query(self):
        q = os.path.join(PARITY, 'mini_queries_nonoverlap.fasta')
        with self.assertRaisesRegex(AlignmentMalformed, "Sequence NO_OVERLAP_QUERY doesn't overlap any "
                                                        "reference sequences"):
            self.placer(q)
        p = self.placer(q, discard_nonoverlapped=True)
        self.assertEqual(p.discarded, ['NO_OVERLAP_QUERY'])
        self.assertNotIn('NO_OVERLAP_QUERY', p.all_names)

    def test_mrca_class_needs_every_leaf_in_seq_info(self):
        from gplace.engine.taxonomy import Classifier
        with tempfile.TemporaryDirectory() as d:
            pkg = testsets._mini_variant('tax_partial', d, 'partial')
            p = self.placer(os.path.join(PARITY, 'mini_queries_dups.fasta'))
            Classifier(pkg, p.tree)                              # edge painting copes with gaps
            with self.assertRaisesRegex(RefpkgMalformed, 'tax_id for .* not found!'):
                Classifier(pkg, p.tree, mrca_class=True)          # as pplacer --mrca-class


if __name__ == '__main__':
    unittest.main()
