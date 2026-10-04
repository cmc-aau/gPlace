import os
import tempfile
import unittest

import numpy as np

from gplace.engine import core, resources
from gplace.engine.newick import parse_newick
from gplace.engine.tiled import tile_layout, plan_tiling, TiledPlacer, MIN_TILE_EDGES, PLAN_FRACTION
from gplace.exceptions import GPlaceExit
from gplace.fileio import read_text

ENGINE = os.path.exists(core.engine_path())
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'gplace', 'tests', 'data')
MINI = os.path.join(DATA, 'mini.refpkg')


def mini_tree():
    return parse_newick(read_text(os.path.join(MINI, 'mini_tree.tree')))


class TestTileLayout(unittest.TestCase):

    def check(self, t, N, bounded=True):
        tiles = tile_layout(t.children, t.postorder(), t.root, N)
        roots = {r: i for i, (r, _, _) in enumerate(tiles)}
        self.assertEqual(tiles[-1][0], t.root)
        seen = [0] * t.n_nodes
        for i, (r, size, slots) in enumerate(tiles):
            self.assertTrue(t.children[r], 'a leaf roots a tile')
            stack, n, n_int, n_below = [r], 0, 0, 0
            while stack:
                u = stack.pop()
                seen[u] += 1; n += 1; n_int += bool(t.children[u])
                for c in t.children[u]:
                    if c in roots:
                        n_below += 1
                        self.assertLess(roots[c], i, 'tiles are not bottom-up')
                    else:
                        stack.append(c)
            self.assertEqual((n, n_int + n_below), (size, slots))
            if bounded and N >= 3:
                self.assertLessEqual(size, N)
        self.assertEqual(seen, [1] * t.n_nodes)
        return tiles

    def test_partition_bounded_bottom_up(self):
        t = mini_tree()
        for N in (2, 3, 7, 40, 100, 492, 493, 10 ** 6):
            self.check(t, N)
        self.assertEqual(len(self.check(t, t.n_nodes)), 1)

    def test_multifurcation(self):
        """x has 5 leaf children, which cannot root tiles: its tile has 6 nodes whatever N."""
        t = parse_newick('((a:1,b:1,c:1,d:1,e:1)x:1,(f:1,(g:1,h:1)y:1)z:1,i:1)r;')
        for N in (2, 3, 4, 5, 6, 20):
            self.check(t, N, bounded=N >= 6)


def balanced_tree(n_leaves):
    """Newick tree with n_leaves leaves, balanced."""
    count = [0]

    def sub(n):
        if n == 1:
            count[0] += 1
            return f'L{count[0]}:0.1'
        return f'({sub(n // 2)},{sub(n - n // 2)}):0.1'
    return parse_newick(sub(n_leaves)[:-4] + ';')


class TestPlan(unittest.TestCase):
    """20,000 reference genomes (39,999 nodes), 5,036 sites, 5,000 queries, 4 workers."""

    @classmethod
    def setUpClass(cls):
        cls.tree = balanced_tree(20000)

    def plan(self, **kw):
        return plan_tiling(self.tree, 5036, 5000, 20000, 5036, 4, **kw)

    def test_untiled_by_default_and_when_it_fits(self):
        p = self.plan()
        self.assertIsNone(p.tile_edges)
        self.assertEqual(p.memory, p.untiled_memory)
        self.assertIsNone(self.plan(max_memory=p.untiled_memory / PLAN_FRACTION + 1).tile_edges)
        self.assertIsNotNone(self.plan(max_memory=p.untiled_memory).tile_edges)    # 5% headroom

    def test_budget_picks_largest_tiles_that_fit(self):
        untiled = self.plan().untiled_memory
        budget = untiled // 3
        p = self.plan(max_memory=budget)
        self.assertIsNotNone(p.tile_edges)
        self.assertLessEqual(p.memory, budget * PLAN_FRACTION)
        self.assertLessEqual(p.max_edges, p.tile_edges)
        self.assertGreater(self.plan(tile_edges=int(p.tile_edges * 1.2)).memory, budget * PLAN_FRACTION)
        smaller = self.plan(max_memory=budget // 2)
        self.assertLess(smaller.tile_edges, p.tile_edges)

    def test_budget_too_small(self):
        with self.assertRaises(GPlaceExit):
            self.plan(max_memory=10 ** 6)

    def test_explicit_tile_edges(self):
        p = self.plan(tile_edges=5000)
        self.assertEqual(p.tile_edges, 5000)
        self.assertLessEqual(p.max_edges, 5000)
        self.assertGreater(p.n_tiles, 7)
        self.assertLess(p.memory, p.untiled_memory)
        with self.assertRaises(GPlaceExit):
            self.plan(tile_edges=1)

    def test_minimum_tile_size(self):
        tiny = self.plan(tile_edges=MIN_TILE_EDGES).memory
        with self.assertRaises(GPlaceExit):
            self.plan(max_memory=tiny - 1)


class TestMerge(unittest.TestCase):

    def test_merge_is_rank_order_with_ties(self):
        """Merging tiles keeps each query's first K edges by (score desc, edge asc), ties included."""
        rng = np.random.default_rng(3)
        n, K = 50, 6
        cs = np.full((n, K), -np.inf); ce = np.full((n, K), np.iinfo(np.int64).max, np.int64)
        all_s, all_e = [], []
        e0 = 0
        for E in (3, 11, 1, 25, 8):
            edges = np.arange(e0, e0 + E, dtype=np.int64)[::-1].copy(); e0 += E
            sc = rng.integers(-5, 0, (n, E)).astype(float)             # many ties
            sc[rng.random((n, E)) < 0.05] = -np.inf
            TiledPlacer._merge(cs, ce, 0, edges, sc)
            all_s.append(sc); all_e.append(np.broadcast_to(edges, (n, E)))
        S, E_ = np.concatenate(all_s, 1), np.concatenate(all_e, 1)
        order = np.lexsort((ce, -cs), axis=1)
        got = np.take_along_axis(ce, order, 1)
        for r in range(n):
            want = E_[r][np.lexsort((E_[r], -S[r]))][:K]
            self.assertEqual(got[r].tolist(), want.tolist())


@unittest.skipUnless(ENGINE, 'engine not built (run: gplace build_engine)')
class TestTiledEngine(unittest.TestCase):

    def run_mini(self, **kw):
        from gplace.engine.jplace import write_jplace, compare_jplace
        tp = TiledPlacer(MINI, core.read_fasta(os.path.join(DATA, 'queries.fasta')), 2, **kw)
        self.assertGreater(len(tp.tiles), 1)
        results = tp.run()
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'tiled.jplace')
            write_jplace(out, tp.jplace_tree(), results, 'unittest')
            r = compare_jplace(out, os.path.join(DATA, 'expected_pplacer.jplace'))
        self.assertTrue(r['same_tree'])
        self.assertEqual(r['identical'], 20)

    def test_tiled_identical_to_pplacer(self):
        self.run_mini(tile_edges=40)

    def test_rerank_and_small_batches_identical(self):
        """Queries running out of ranked candidates are re-ranked; tiny query batches."""
        self.run_mini(tile_edges=60, n_cand=2, round_size=1, query_batch=7)

    def test_benchmark_timings(self):
        """The benchmark's probe and work recording leave the placements identical; the time
        model grows with the number of queries and adds up its parts."""
        from gplace.engine.jplace import write_jplace, compare_jplace
        tp = TiledPlacer(MINI, core.read_fasta(os.path.join(DATA, 'queries.fasta')), 2, tile_edges=60)
        tp.record_work = True
        tp.cache()
        tp.prefilter(probe=5)
        results = tp.place_all()
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'tiled.jplace')
            write_jplace(out, tp.jplace_tree(), results, 'unittest')
            self.assertEqual(compare_jplace(out, os.path.join(DATA, 'expected_pplacer.jplace'))['identical'], 20)
        ti, n_ev, rates = tp.eval_rates([1, 2], min_evals=50)
        self.assertIn(ti, tp.work_log)
        self.assertGreaterEqual(n_ev, 50)
        self.assertEqual([w for w, _, _ in rates], [1, 2])
        self.assertTrue(all(dt > 0 and eps > 0 for _, dt, eps in rates))
        tp.release()
        eps = rates[-1][2]
        total, boundary, pref, place = tp.job_seconds(1000, eps)
        self.assertAlmostEqual(total, boundary + pref + place)
        self.assertGreater(tp.job_seconds(100000, eps)[0], total)
        self.assertGreater(tp.evals_per_query(), 1)


if __name__ == '__main__':
    unittest.main()
