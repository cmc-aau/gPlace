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

    def test_engine_rebuild_matches_python(self):
        """The engine's cache build and tile rebuilds (Placer._rebuild: a block of sites per thread;
        one chunk of transition matrices, or chunks of one or two nodes) hold the same bits as a
        Python rebuild node by node with the engine's kernels (protein and DNA CAT, gamma, gamma with
        --fix_gamma_underflow, an IQ-TREE mixture), untiled and tiled."""
        import shutil
        from gplace.engine import placer as placer_mod, testsets
        from gplace.engine.factory import new_placer, tiled_placer
        from gplace.engine.reference import Reference
        models, fits = os.path.join(DATA, 'models'), os.path.join(DATA, 'fits')
        q_mini, q_dna = os.path.join(DATA, 'queries.fasta'), os.path.join(models, 'like_gtr_queries.fasta')
        tmp = tempfile.mkdtemp()
        mixture = Reference(tree=os.path.join(fits, 'full_mix_LG4X.treefile'), ref_aln=os.path.join(MINI, 'mini_msa.fa'),
                            stats=os.path.join(fits, 'full_mix_LG4X.iqtree'), full_model='auto')
        cases = (('protein CAT', MINI, q_mini, 40, False), ('DNA CAT', os.path.join(models, 'dna_cat.refpkg'), q_dna, 15, False),
                 ('protein gamma', testsets._mini_variant('WAG_g0.5x4', tmp), q_mini, 40, False),
                 ('gamma + fix', testsets._mini_variant('WAG_g0.05x4', tmp), q_mini, 40, True),
                 ('DNA gamma', os.path.join(models, 'like_gtr.refpkg'), q_dna, 15, False),
                 ('mixture', mixture, q_mini, 40, False))
        saved = placer_mod.TENSOR_BYTES, placer_mod.BLOCK_SITES

        def bits(a):
            return np.ascontiguousarray(a).view(np.uint8)

        def py_distal(pl, v):                       # the Python cache build the engine replaced
            t, S, ch, k = pl.tl(), pl.nsites, pl.children[v], pl.iidx[v]
            tA, tB, ten = t.b[0], t.b[1], t.tensor
            a, e = pl.d_pre(ch[0]); pl.evolve_t(tA, a, pl.bl[ch[0]], ten)
            if len(ch) == 1:
                pl.darr[k] = tA; pl.dexp[k] = e
            else:
                b, f = pl.d_pre(ch[1]); pl.evolve_t(tB, b, pl.bl[ch[1]], ten)
                pl.lib.xp_pairwise(core.P(pl.darr[k]), core.P(pl.dexp[k]), core.P(tA), core.P(e), core.P(tB), core.P(f), S)
                for c in ch[2:]:
                    a2, e2 = pl.d_pre(c); pl.evolve_t(tA, a2, pl.bl[c], ten)
                    pl.lib.xp_pairwise(core.P(pl.darr[k]), core.P(pl.dexp[k]), core.P(pl.darr[k]),
                                       core.P(pl.dexp[k]), core.P(tA), core.P(e2), S)
            pl.lib.xp_pull_exponent(core.P(pl.darr[k]), core.P(pl.dexp[k]), S, -50)

        def py_rebuild(pl, top, internal, need_p):
            post, stack = [], [(top, False)]
            while stack:
                v, seen = stack.pop()
                if seen:
                    post.append(v)
                    continue
                stack.append((v, True))
                stack.extend((c, False) for c in reversed(pl.children[v]) if c in internal)
            for v in post:
                py_distal(pl, v)
            if need_p:
                t = pl.tl()
                for v in reversed(post):            # parents first
                    kv, ch = pl.iidx[v], pl.children[v]
                    evp = t.b[2]
                    pl.evolve_t(evp, pl.parr[kv], pl.bl[v], t.tensor)
                    for c in ch:
                        if pl.children[c]:
                            kc = pl.iidx[c]
                            pl._prox_into(pl.parr[kc], pl.pexp[kc], evp, pl.pexp[kv], [r for r in ch if r != c],
                                          t.b[3], t.tensor)

        def snapshot(pl, n):
            return [bits(x[:n].copy()) for x in (pl.darr, pl.dexp, pl.parr, pl.pexp)]

        try:
            for (label, ref, queries, edges, fix), small in [(c, s) for c in cases for s in (False, True)]:
                with self.subTest(model=label, small_chunks=small):
                    placer_mod.BLOCK_SITES = 7                              # several site blocks
                    pl = new_placer(ref, core.read_fasta(queries), 3, fix_gamma_underflow=fix)
                    placer_mod.TENSOR_BYTES = 3 * pl.nrates * pl.ns * pl.ns * 8 if small else saved[0]
                    pl.cache()                                              # untiled: the engine
                    n = len(pl.internal)
                    engine = snapshot(pl, n)
                    py_rebuild(pl, pl.root, set(pl.internal), True)        # the same arrays, in Python
                    for a, b in zip(engine, snapshot(pl, n)):
                        self.assertTrue(np.array_equal(a, b))
                    tp = tiled_placer(pl, edges)
                    self.assertGreater(len(tp.tiles), 2)
                    tp.cache()                                              # boundary pass: the engine
                    with tp._rebuild_threads():
                        for ti in reversed(range(len(tp.tiles))):          # top-down, as the prefilter
                            tile = tp.tiles[ti]
                            n = len(tile['internal']) + len(tile['child_roots'])
                            tp._materialise(ti, need_p=True)
                            engine = snapshot(tp, n)
                            py_rebuild(tp, tile['root'], set(tile['internal']), True)
                            for a, b in zip(engine, snapshot(tp, n)):
                                self.assertTrue(np.array_equal(a, b))
                    tp.release()
        finally:
            placer_mod.TENSOR_BYTES, placer_mod.BLOCK_SITES = saved
            shutil.rmtree(tmp)

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
