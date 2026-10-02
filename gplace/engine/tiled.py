###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""Exact tiled placement, for reference trees whose likelihood vectors do not fit in memory.

The tree is cut into connected tiles of at most `tile_edges` nodes. A tile touches the rest of the
tree only through its root edge (the tree above it) and the roots of the tiles below it. Once the
distal vector of every tile root (bottom-up pass) and the proximal vector of every tile root
(top-down pass) are stored, any tile can be rebuilt on its own with the same operations, operands
and order as the full cache, so every vector is bit-identical to the untiled one while only one
tile is held in memory.

pplacer's algorithm then runs tile by tile, with the work reordered but not changed:
  * prefilter: every query is scored against every edge of a tile (queries in batches) and keeps
    its first `n_cand` edges in pplacer's ranking order (score descending, ties by ascending edge
    number): pplacer's global ranking without a queries x edges score table;
  * baseball in rounds: every first-stage optimisation starts from pplacer's default branch
    lengths, so its result does not depend on which edges were tried before. The next
    `round_size` ranked edges of every query still searching are optimised, grouped by tile, and
    pplacer's stopping rule is replayed on the stored results; a query that would need more than
    its `n_cand` edges is re-ranked with more;
  * the kept edges are refined at the final tolerance in their tiles; like-weight ratios as in
    pplacer.
"""

import gc
import json
import math
import os
import threading
import time
from collections import defaultdict, namedtuple
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from gplace.engine import core, resources
from gplace.engine import placer as placer_mod
from gplace.engine.core import P
from gplace.engine.newick import parse_newick
from gplace.engine.placer import Placer, Baseball, weighted_placements, shared_zeros, _fork_pool
from gplace.exceptions import GPlaceException, GPlaceExit

MIN_TILE_EDGES = 1000                   # smallest tiles --max_memory chooses
PLAN_FRACTION = 0.95                    # --max_memory plans keep 5% headroom below the limit
_NO_EDGE = np.iinfo(np.int64).max       # filler in the candidate lists: ranks after every edge

TilePlan = namedtuple('TilePlan', 'tile_edges n_tiles max_edges n_slots memory untiled_memory')
TilePlan.__doc__ = """How one reference package is placed: tile_edges None = untiled.
memory / untiled_memory: estimated peak bytes of the chosen / the untiled placement."""


def _tile_prefilter_group(g):
    placer_mod._PLACER._tile_prefilter_group(g)


def _tile_eval(job):
    return placer_mod._PLACER._eval_job(job)


# ---------------------------------------------------------------------- tiling plan --
def tile_layout(children, order, root, tile_edges):
    """Cut the tree into connected tiles of at most tile_edges nodes (edges). Walking up in
    post-order, the part of a node's subtree not yet in a tile stays whole while it has at most
    tile_edges nodes; when it would have more, its internal children with the largest parts
    become tile roots first. Leaves never root a tile (a leaf's proximal vector is rebuilt from
    its parent's), so a node with more than tile_edges - 1 leaf children exceeds the bound. The
    root closes the last tile. Returns [(tile root, nodes, slots)] bottom-up (every tile after
    the tiles below it), where slots = the tile's internal nodes plus the roots of the tiles
    below it: the vectors held while the tile is in memory."""
    n = len(children)
    size, slots = [0] * n, [0] * n
    tiles = []
    for v in order:
        ch = children[v]
        s = 1 + sum(size[c] for c in ch)
        k = (1 if ch else 0) + sum(slots[c] for c in ch)
        if s > tile_edges:
            for c in sorted((c for c in ch if children[c]), key=lambda c: -size[c]):
                tiles.append((c, size[c], slots[c]))
                s -= size[c]; k += 1 - slots[c]
                if s <= tile_edges:
                    break
        size[v], slots[v] = s, k
    tiles.append((root, size[root], slots[root]))
    return tiles


def plan_tiling(tree, n_sites, n_queries, n_leaves, aln_len, workers, tile_edges=None, max_memory=None,
                width=resources.NS):
    """Choose how to place n_queries on `tree` (a newick.Tree).

    tile_edges : place on tiles of at most this many edges.
    max_memory : memory limit in bytes: untiled if its estimated peak fits, otherwise the largest
                 tiles whose estimated peak fits (larger tiles are faster: fewer rebuilds and
                 more work per worker pool). Plans use PLAN_FRACTION of the limit as headroom.
    Neither    : untiled.
    width      : doubles per site in a likelihood vector (20 for protein CAT; gamma: rates x states)."""
    n_internal = tree.n_nodes - n_leaves
    untiled = resources.estimate_memory(tree.n_nodes, n_internal, n_sites, n_queries, n_leaves,
                                        aln_len, workers, width=width)
    whole = TilePlan(None, 1, tree.n_nodes, n_internal, untiled, untiled)
    budget = None if max_memory is None else max_memory * PLAN_FRACTION
    if tile_edges is None and (budget is None or untiled <= budget):
        return whole
    order = tree.postorder()

    def plan(n):
        tiles = tile_layout(tree.children, order, tree.root, n)
        mx, sl = max(t[1] for t in tiles), max(t[2] for t in tiles)
        mem = resources.estimate_memory_tiled(sl, mx, len(tiles), n_sites, n_queries, n_leaves,
                                              aln_len, workers, width=width)
        return TilePlan(n, len(tiles), mx, sl, mem, untiled)

    if tile_edges is not None:
        if tile_edges < 2:
            raise GPlaceExit('--tile_edges must be at least 2.')
        return plan(int(tile_edges))
    lo, hi = min(MIN_TILE_EDGES, tree.n_nodes), tree.n_nodes
    best = plan(lo)
    if best.memory > budget:
        GiB = resources.GiB
        need = min(best.memory, untiled) / PLAN_FRACTION
        raise GPlaceExit(f'--max_memory {max_memory / GiB:.1f} GiB is too small: placing {n_queries:,} '
                         f'queries on this reference package needs --max_memory '
                         f'{need / GiB:.1f}G or more (estimated peak: untiled {untiled / GiB:.1f} GiB, '
                         f'smallest tiles of {lo:,} edges {best.memory / GiB:.1f} GiB). Place fewer '
                         f'queries per run, or use fewer --cpus.')
    while lo < hi:                          # largest tile size that fits
        mid = (lo + hi + 1) // 2
        p = plan(mid)
        if p.memory <= budget:
            lo, best = mid, p
        else:
            hi = mid - 1
    return best


def plan_for_placer(placer, tile_edges=None, max_memory=None):
    """plan_tiling for a loaded Placer (its masked sites and queries)."""
    workers = placer.workers or resources.heuristic_cpus(placer.cpus)
    return plan_tiling(placer.tree, placer.nsites, len(placer.query_names), placer.n_leaves,
                       placer.aln_len, workers, tile_edges, max_memory, width=placer.width)


def load_tree(refpkg_dir):
    """The reference tree of a refpkg (for planning without loading the alignments)."""
    files = json.load(open(os.path.join(refpkg_dir, 'CONTENTS.json')))['files']
    return parse_newick(open(os.path.join(refpkg_dir, files['tree'])).read())


# -------------------------------------------------------------------- tiled placer --
class TiledPlacer(Placer):
    """Placer holding one tile of the likelihood cache at a time; same results as Placer."""

    def __init__(self, refpkg_dir, queries, cpus=1, tile_edges=10000, **kw):
        super().__init__(refpkg_dir, queries, cpus)
        self._setup_tiles(tile_edges, **kw)

    @classmethod
    def from_placer(cls, placer, tile_edges, **kw):
        """Tile a loaded Placer (the reference package and queries are not read again)."""
        self = cls.__new__(cls)
        self.__dict__.update(placer.__dict__)
        self._setup_tiles(tile_edges, **kw)
        return self

    def _setup_tiles(self, tile_edges, query_batch=resources.QUERY_BATCH, n_cand=resources.N_CAND,
                     round_size=8):
        if int(tile_edges) < 2:
            raise GPlaceException('tile_edges must be at least 2.')
        self.tile_edges = int(tile_edges)
        self.query_batch = max(1, int(query_batch))
        self.n_cand = max(1, int(n_cand))
        self.round_size = max(1, int(round_size))
        layout = tile_layout(self.children, self.tree.postorder(), self.root, self.tile_edges)
        roots = {r for r, _, _ in layout}
        self.tile_of = [-1] * self.nnodes
        self.tiles = []
        for ti, (r, _, _) in enumerate(layout):            # bottom-up order
            nodes, child_roots, stack = [], [], [r]
            while stack:
                u = stack.pop()
                nodes.append(u)
                self.tile_of[u] = ti
                for c in self.children[u]:
                    (child_roots if c in roots else stack).append(c)
            self.tiles.append({'root': r, 'nodes': nodes, 'child_roots': child_roots,
                               'internal': [u for u in nodes if self.children[u]]})
        self.n_slots = max(len(t['internal']) + len(t['child_roots']) for t in self.tiles)
        self.max_edges = max(len(t['nodes']) for t in self.tiles)
        self.stats = defaultdict(float)           # phase timings (seconds) and counts, for benchmark
        self.record_work = False                  # benchmark: keep every placement job, by tile
        self.work_log = defaultdict(list)

    def memory_estimate(self, workers=None):
        """Estimated peak memory (bytes) of this pass (see resources.estimate_memory_tiled)."""
        w = workers or self.workers or resources.heuristic_cpus(self.cpus)
        return resources.estimate_memory_tiled(self.n_slots, self.max_edges, len(self.tiles),
                                               self.nsites, len(self.query_names), self.n_leaves,
                                               self.aln_len, w, self.query_batch, self.n_cand,
                                               width=self.width)

    def _pool_map(self, fn, items, chunksize=1):
        """fn over items in a pool of self.workers forked processes; returns (results, seconds,
        seconds without the pool start-up)."""
        placer_mod._PLACER = self
        t0 = time.time()
        try:
            with _fork_pool(self.workers) as pool:
                t1 = time.time()
                out = list(pool.imap_unordered(fn, items, chunksize=chunksize))
                t2 = time.time()
            return out, time.time() - t0, t2 - t1
        finally:
            placer_mod._PLACER = None
            gc.unfreeze()

    # ------------------------------------------------------------------ tiles --
    def _materialise(self, ti, need_p=True):
        """Rebuild tile ti's distal (and proximal) vectors from the stored boundary vectors, with
        the operations of Placer.cache; stores its root's distal and, with need_p, the proximal
        vectors of the tile roots below it."""
        tile = self.tiles[ti]
        slots = tile['internal'] + tile['child_roots']
        self.iidx = {v: k for k, v in enumerate(slots)}
        self._tile_in_memory = ti
        for c in tile['child_roots']:
            k = self.iidx[c]
            self.darr[k] = self.bD[c]; self.dexp[k] = self.bDe[c]
        internal = set(tile['internal'])
        n = len(tile['internal'])
        lock = threading.Lock(); done = threading.Event(); errs = []; count = [0]
        ex = ThreadPoolExecutor(self.cpus)
        pending = {v: sum(1 for c in self.children[v] if c in internal) for v in tile['internal']}

        def distal_task(v):
            try:
                self._distal(v, self.tl())
                p = self.parent[v]; ready = False
                with lock:
                    count[0] += 1; fin = count[0] == n
                    if p in internal:
                        pending[p] -= 1; ready = pending[p] == 0
                if ready:
                    ex.submit(distal_task, p)
                if fin:
                    done.set()
            except BaseException as e:
                errs.append(e); done.set()

        try:
            for v in tile['internal']:
                if pending[v] == 0:
                    ex.submit(distal_task, v)
            done.wait()
            if errs:
                raise errs[0]
            r = tile['root']
            kr = self.iidx[r]
            self.bD[r] = self.darr[kr].copy(); self.bDe[r] = self.dexp[kr].copy()
            if not need_p:
                return
            if r == self.root:
                self.parr[kr] = 1.0; self.pexp[kr] = 0
            else:
                self.parr[kr] = self.bP[r]; self.pexp[kr] = self.bPe[r]
            count[0] = 0; done.clear()

            def prox_task(v):
                try:
                    kids = self._proximal(v, self.tl())      # also fills the child tile roots
                    with lock:
                        count[0] += 1; fin = count[0] == n
                    for c in kids:
                        if c in internal:
                            ex.submit(prox_task, c)
                    if fin:
                        done.set()
                except BaseException as e:
                    errs.append(e); done.set()

            ex.submit(prox_task, r)
            done.wait()
            if errs:
                raise errs[0]
            for c in tile['child_roots']:
                k = self.iidx[c]
                self.bP[c] = self.parr[k].copy(); self.bPe[c] = self.pexp[k].copy()
        finally:
            ex.shutdown()

    def cache(self):
        """Allocate one tile's vectors; store the distal vector of every tile root (bottom-up)."""
        S = self.nsites
        self.zero_e = np.zeros(S, np.int64)
        self.darr = shared_zeros((self.n_slots, S, self.width)); self.dexp = shared_zeros((self.n_slots, S), np.int64)
        self.parr = shared_zeros((self.n_slots, S, self.width)); self.pexp = shared_zeros((self.n_slots, S), np.int64)
        self.bD, self.bDe, self.bP, self.bPe = {}, {}, {}, {}
        t0 = time.time()
        for ti in range(len(self.tiles)):
            self._materialise(ti, need_p=False)
        self.stats['boundary_s'] = time.time() - t0

    # ---------------------------------------------------------------- prefilter --
    def prefilter(self, probe=0):
        """pplacer's ranking for every query: its first n_cand edges, tile by tile. probe
        (benchmark): also score each tile for the first `probe` queries alone, which separates
        the prefilter's fixed and per-query costs (see job_seconds)."""
        self._query_tables()
        if self.workers is None:                     # --cpus auto without a saved benchmark
            self.workers = resources.heuristic_cpus(self.cpus)
            self.cpus_source = 'heuristic'
        self.rank = self._ranked_candidates(np.arange(len(self.query_names)), self.n_cand,
                                            self.stats, probe)
        self.rank_more = {}

    def _ranked_candidates(self, qidx, K, stats=None, probe=0):
        """The first K edges of pplacer's ranking (h_ranking: score descending, ties by ascending
        edge number) for the queries qidx, as an (n, K) array padded with _NO_EDGE."""
        st = stats if stats is not None else defaultdict(float)
        n = len(qidx)
        cs = np.full((n, K), -np.inf)
        ce = np.full((n, K), _NO_EDGE, np.int64)
        qb = min(self.query_batch, n)
        self.scores = shared_zeros((qb, self.max_edges))
        t0, nt = time.time(), len(self.tiles)
        st['pref_probe_n'] = probe
        try:
            for step, ti in enumerate(reversed(range(nt)), 1):   # top-down: stores boundary proximals
                tm = time.time()
                self._materialise(ti, need_p=True)
                st['pref_rebuild_s'] += time.time() - tm
                edges = np.array([v for v in self.tiles[ti]['nodes'] if v != self.root], np.int64)
                if not len(edges):
                    continue
                self.locs = edges
                E = len(edges)
                groups = [(i, min(E, i + 8)) for i in range(0, E, 8)]
                for b0 in range(0, n, qb):
                    rows = qidx[b0:b0 + qb]
                    if probe and b0 == 0:           # benchmark only: its scores are overwritten
                        self._bcodes = np.ascontiguousarray(self.codes[rows[:probe]])
                        self._qb = (len(rows[:probe]), E)
                        st['pref_probe_s'] += self._pool_map(_tile_prefilter_group, groups, chunksize=2)[1]
                    self._bcodes = np.ascontiguousarray(self.codes[rows])
                    self._qb = (len(rows), E)
                    st['pref_batch_s'] += self._pool_map(_tile_prefilter_group, groups, chunksize=2)[1]
                    tm = time.time()
                    self._merge(cs, ce, b0, edges, self.scores[:len(rows), :E])
                    st['pref_merge_s'] += time.time() - tm
                self.logger.info(f'Prefilter: tile {step}/{nt} ({E:,} edges) done '
                                 f'({time.time() - t0:.0f} s).')
        finally:
            del self.scores
            self._bcodes = None
        st['pref_s'] += time.time() - t0
        order = np.lexsort((ce, -cs), axis=1)
        return np.take_along_axis(ce, order, 1)

    def _tile_prefilter_group(self, g):
        i0, i1 = g
        nq, _ = self._qb
        t = self._supernodes8(i0, i1)
        st = self.lib.xp_edge_scores8(P(self.scores), self.max_edges, i0, i1 - i0, P(t.sn8), P(t.se8),
                                      P(self.qevtab), P(self.cats), self.nsites, len(self.alpha.codes),
                                      P(self._bcodes), nq)
        if st != 0:
            raise GPlaceException(f'Prefilter kernel failed ({st}) for edges {i0}..{i1}.')

    @staticmethod
    def _merge(cs, ce, b0, edges, sc):
        """Merge one tile's scores sc (rows b0.. of cs/ce) into the candidate lists, keeping each
        query's first K edges in ranking order. The selection is by score; when equal scores
        straddle the cut-off, that row is resolved by (score, edge number) exactly."""
        K, E = cs.shape[1], len(edges)
        step = max(1, resources.MERGE_CELLS // (K + E))
        for a in range(0, sc.shape[0], step):
            b = min(sc.shape[0], a + step)
            s = np.concatenate([cs[b0 + a:b0 + b], sc[a:b]], axis=1)
            e = np.concatenate([ce[b0 + a:b0 + b], np.broadcast_to(edges, (b - a, E))], axis=1)
            part = np.argpartition(-s, K - 1, axis=1)[:, :K]
            kth = np.take_along_axis(s, part, 1).min(axis=1)
            for r in np.nonzero((s >= kth[:, None]).sum(axis=1) > K)[0]:
                part[r] = np.lexsort((e[r], -s[r]))[:K]
            cs[b0 + a:b0 + b] = np.take_along_axis(s, part, 1)
            ce[b0 + a:b0 + b] = np.take_along_axis(e, part, 1)

    def _rank_row(self, qi):
        return self.rank_more.get(qi, self.rank[qi])

    def _extend_ranking(self, qis):
        """Re-rank queries that used up their candidate edges, with 4x as many (up to all)."""
        K = min(self.nnodes - 1, 4 * max(len(self._rank_row(qi)) for qi in qis))
        self.logger.info(f'{len(qis):,} queries need more than their ranked candidate edges '
                         f'(many failed optimisations); ranking {K:,} edges for them.')
        rows = self._ranked_candidates(np.array(qis, np.int64), K)
        for qi, row in zip(qis, rows):
            old = self._rank_row(qi)
            if not np.array_equal(row[:len(old)], old):
                raise GPlaceException(f'Inconsistent edge ranking for {self.query_names[qi]}.')
            self.rank_more[qi] = row

    # ---------------------------------------------------------------- placement --
    def _eval_job(self, job):
        qi, items = job
        evaluate, free = self._query_evaluator(qi)
        try:
            return qi, [(loc,) + tuple(evaluate(loc, pend, dist, tol)) for loc, pend, dist, tol in items]
        finally:
            free()

    def _evaluate(self, work):
        """Run the optimisations work = {qi: [(loc, pend, dist, tol)]}, tile by tile; returns
        [(qi, [(loc, result, status)])]."""
        by_tile = defaultdict(lambda: defaultdict(list))
        for qi, items in work.items():
            for it in items:
                by_tile[self.tile_of[it[0]]][qi].append(it)
        out = []
        for ti in sorted(by_tile, reverse=True):
            self._materialise(ti, need_p=True)
            jobs = list(by_tile[ti].items())
            if self.record_work:
                self.work_log[ti].extend(jobs)
            res, _, inner = self._pool_map(_tile_eval, jobs,
                                           chunksize=max(1, min(8, len(jobs) // (16 * self.workers))))
            out.extend(res)
            self.stats['eval_s'] += inner
            self.stats['evals'] += sum(len(items) for _, items in jobs)
        return out

    def _replay(self, qi, stage1):
        """pplacer's baseball on the stored first-stage results in ranked order: the Baseball when
        the search has ended, 'next' when the next ranked edge is not evaluated yet, 'more' when
        the query needs more ranked edges than it has."""
        bb = Baseball()
        row = self._rank_row(qi)
        for loc in row:
            if loc == _NO_EDGE:                      # every edge of the tree tried
                return bb
            loc = int(loc)
            if loc not in stage1:
                return 'next'
            r, st = stage1[loc]
            if not bb.add(loc, r, st):
                return bb
        return bb if len(row) >= self.nnodes - 1 else 'more'

    def place_all(self, workers=None):
        """Place every query; returns [(name, placements, warnings)] in input order."""
        if workers:
            self.workers = workers
        t_start = time.time()
        Q = len(self.query_names)
        stage1 = [dict() for _ in range(Q)]
        nxt = [0] * Q
        active, done, rnd = list(range(Q)), {}, 0
        while active:
            rnd += 1
            work = {}
            for qi in active:
                row = self._rank_row(qi)
                chunk = [int(v) for v in row[nxt[qi]:nxt[qi] + self.round_size] if v != _NO_EDGE]
                nxt[qi] += len(chunk)
                work[qi] = [(loc, core.START_PEND, self.bl[loc] / 2., core.INITIAL_TOL) for loc in chunk]
            for qi, out in self._evaluate(work):
                for loc, r, st in out:
                    stage1[qi][loc] = (r, st)
            still, more = [], []
            for qi in active:
                bb = self._replay(qi, stage1[qi])
                if isinstance(bb, Baseball):
                    done[qi] = bb
                else:
                    still.append(qi)
                    if bb == 'more':
                        more.append(qi)
            if more:
                self._extend_ranking(more)
            active = still
            self.logger.info(f'Baseball round {rnd}: {Q - len(active):,}/{Q:,} queries complete.')

        kept, work = {}, {}
        for qi in range(Q):
            if not done[qi].results:
                raise GPlaceException(f'empty results for {self.query_names[qi]}!')
            kept[qi] = done[qi].kept(self.lib)
            work[qi] = [(loc, pend, dist, core.FINAL_TOL) for loc, (lk, pend, dist) in kept[qi]]
        refined = {}
        for qi, out in self._evaluate(work):
            for loc, r, _ in out:
                refined[(qi, loc)] = r
        del work, stage1
        results = []
        for qi in range(Q):
            ref = []
            for loc, (lk, pend, dist) in kept[qi]:
                r = refined[(qi, loc)]
                ref.append((loc, r if r is not None else (lk, pend, dist)))
            results.append((self.query_names[qi], weighted_placements(self.lib, ref), done[qi].warns))
        self.stats['rounds'] = rnd
        self.stats['place_s'] += time.time() - t_start
        self.logger.info('Refined the kept placements.')
        return results

    # ---------------------------------------------------------------- benchmark --
    def eval_rates(self, worker_counts, min_evals=4000):
        """Placement speed at each worker count, without tile rebuilds: the tile with the most
        recorded optimisations (record_work) is rebuilt once and its optimisations are replayed,
        repeated to at least min_evals. Returns (tile, evaluations, [(workers, seconds,
        evaluations per second)]); pool start-up is excluded, as in Placer._run_groups."""
        if not self.work_log:
            raise GPlaceException('No recorded placement work: run place_all with record_work set.')
        ti = max(self.work_log, key=lambda t: sum(len(items) for _, items in self.work_log[t]))
        n = sum(len(items) for _, items in self.work_log[ti])
        reps = max(1, math.ceil(min_evals / n))
        jobs = self.work_log[ti] * reps
        self._materialise(ti, need_p=True)
        saved, out = self.workers, []
        try:
            for w in worker_counts:
                self.workers = w
                dt = self._pool_map(_tile_eval, jobs)[2]
                out.append((w, dt, n * reps / max(dt, 1e-9)))
        finally:
            self.workers = saved
        return ti, n * reps, out

    def evals_per_query(self):
        return self.stats['evals'] / len(self.query_names)

    def scoring_seconds(self, batch, small, max_groups):
        """Prefilter scoring cost per (query, edge) at `batch` queries per batch (wall seconds with
        self.workers). The tile in memory is scored for its first max_groups groups of 8 edges,
        with the sample's queries repeated to `batch` rows and to `small` rows; the difference is
        the per-query part (the supernodes and the pool start-up cancel). Measuring at the job's
        batch size matters: the queries' codes then come from main memory, not from cache."""
        ti = self._tile_in_memory
        edges = np.array([v for v in self.tiles[ti]['nodes'] if v != self.root], np.int64)
        groups = [(i, min(len(edges), i + 8)) for i in range(0, len(edges), 8)][:max_groups]
        n_edges = groups[-1][1]
        self.locs = edges
        times = []
        try:
            for b in (small, batch):
                rows = np.resize(np.arange(len(self.query_names)), b)
                self._bcodes = np.ascontiguousarray(self.codes[rows])
                self._qb = (b, len(edges))
                self.scores = shared_zeros((b, self.max_edges))
                times.append(self._pool_map(_tile_prefilter_group, groups, chunksize=2)[1])
                del self.scores
        finally:
            self._bcodes = None
        return max(0.0, times[1] - times[0]) / ((batch - small) * n_edges)

    def job_seconds(self, n_queries, evals_per_second, score_seconds=None):
        """Estimated time of the same tiled job for n_queries, from this run on a sample (cache,
        prefilter(probe=...), place_all). The sample's fixed costs are kept: tile rebuilds, pool
        start-ups, the prefilter's supernodes (once per batch of query_batch queries) and the
        baseball rounds; its per-query costs are scaled: prefilter scoring (score_seconds per
        query and edge, from scoring_seconds, else the sample's own cost; the probe separates it
        from the fixed part), merging, and the optimisations at evals_per_second.
        Returns (total, boundary, prefilter, placement) seconds."""
        s, n = self.stats, len(self.query_names)
        p = int(s['pref_probe_n'])
        if p and n > p:
            per_q = max(0.0, (s['pref_batch_s'] - s['pref_probe_s']) / (n - p))   # scoring
            sweep = max(0.0, s['pref_probe_s'] - p * per_q)         # supernodes + pool start-ups
        else:
            per_q, sweep = s['pref_batch_s'] / n, 0.0
        fixed = max(0.0, s['pref_s'] - s['pref_probe_s'] - s['pref_merge_s'] - n * per_q)
        batches = math.ceil(n_queries / self.query_batch)
        score_q = per_q if score_seconds is None else score_seconds * (self.nnodes - 1)
        t_pref = fixed + (batches - 1) * sweep + n_queries * (score_q + s['pref_merge_s'] / n)
        t_place = (max(0.0, s['place_s'] - s['eval_s'])
                   + n_queries * self.evals_per_query() / evals_per_second)
        return s['boundary_s'] + t_pref + t_place, s['boundary_s'], t_pref, t_place

    def release(self):
        for attr in ('bD', 'bDe', 'bP', 'bPe', 'rank', 'rank_more', '_bcodes'):
            if hasattr(self, attr):
                delattr(self, attr)
        super().release()

    def run(self):
        """Boundary pass -> prefilter -> baseball rounds -> refinement, tile by tile."""
        t0 = time.time()
        self.cache()
        self.logger.info(f'Boundary vectors of {len(self.tiles):,} tiles computed '
                         f'({(self.darr.nbytes + self.parr.nbytes) / 1e9:.1f} GB of vectors per tile; '
                         f'{time.time() - t0:.0f} s).')
        t0 = time.time()
        self.prefilter()
        self.logger.info(f'Ranked {self.nnodes - 1:,} edges for {len(self.query_names):,} queries '
                         f'({time.time() - t0:.0f} s).')
        results = self.place_all()
        self.release()
        return results
