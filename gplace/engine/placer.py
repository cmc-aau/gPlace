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

"""Exact, memory-bounded re-implementation of `pplacer -c <refpkg>` (v1.1.alpha19), for every model
pplacer supports (WAG/LG/JTT/GTR, CAT or gamma rates; see models.py; gamma in gamma.py).

Memory: pplacer caches distal, proximal and supernode vectors for every node
(3 x nodes x sites x 20 doubles, ~104 GB for the largest GTDB-Tk R226 class-level tree).
gPlace stores distal/proximal vectors for internal nodes only, in their state before pplacer's
final "pull exponents (-10)" pass, and re-applies that pass to a copy at use time; power-of-two
rescaling is exact, so every value is bit-identical. Leaf vectors are rebuilt on demand in
pplacer's operation order and supernodes are computed during the prefilter sweep.

Speed: the likelihood-vector cache is built by a dependency-driven thread pool; the prefilter
and the placements run in forked worker processes, with the large arrays in shared memory.
"""

import gc
import json
import logging
import math
import mmap
import multiprocessing
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from gplace.engine import core, models
from gplace.engine.core import P
from gplace.engine import resources
from gplace.engine.newick import parse_newick, to_newick
from gplace.exceptions import RefpkgMalformed, AlignmentMalformed, GPlaceException


def shared_zeros(shape, dtype=np.float64):
    """Zero-filled array in shared anonymous memory: visible to forked workers, never copied
    on fork, charged once to a memory cgroup."""
    count = int(np.prod(shape))
    buf = mmap.mmap(-1, max(count * np.dtype(dtype).itemsize, 1),
                    flags=mmap.MAP_SHARED | mmap.MAP_ANONYMOUS)
    return np.frombuffer(buf, dtype=dtype, count=count).reshape(shape)


_PLACER = None                                  # the Placer seen by forked workers


def _fork_pool(n):
    gc.collect()
    gc.freeze()                                 # keep object pages shared after fork
    return multiprocessing.get_context('fork').Pool(n)


def _prefilter_group(g):
    _PLACER._prefilter_group(g)


def _place_chunk(qis):
    return [(qi,) + tuple(_PLACER.place_query(qi)) for qi in qis]


class Baseball(object):
    """pplacer's stopping rule for the first-stage optimisations ("baseball"), fed one location at
    a time in ranked order: stop after MAX_STRIKES locations more than STRIKE_BOX log-units below
    the best so far, or after MAX_PITCHES locations."""

    def __init__(self):
        self.results, self.warns = [], []
        self.like_record, self.n_strikes = -math.inf, 0

    def add(self, loc, r, st):
        """Record one location's result (r = (like, pend, dist), or None with GSL status st);
        returns False when pplacer tries no further locations."""
        if r is None:
            self.warns.append((loc, st))
            return True
        like = r[0]; prev = len(self.results); self.results.append((loc, r))
        if prev >= core.MAX_PITCHES:
            return False
        elif like > self.like_record:
            self.like_record = like
        elif like < self.like_record - core.STRIKE_BOX:
            if self.n_strikes + 1 >= core.MAX_STRIKES:
                return False
            self.n_strikes += 1
        return True

    def kept(self, lib):
        """Locations refined at the final tolerance: up to KEEP_AT_MOST within log(KEEP_FACTOR)
        of the best, best first."""
        srt = sorted(self.results, key=lambda x: -x[1][0])
        best = srt[0][1][0]
        lkf = lib.xp_log(core.KEEP_FACTOR)
        return [x for i, x in enumerate(srt) if i < core.KEEP_AT_MOST and x[1][0] >= lkf + best]


def weighted_placements(lib, refined):
    """Like-weight ratios and pplacer's placement order for [(loc, (like, pend, dist))];
    returns [[distal_length, edge_num, like_weight_ratio, likelihood, pendant_length]]."""
    lls = [r[0] for _, r in refined]
    lwr = [1. / sum(lib.xp_exp(o - ll) for o in lls) for ll in lls]
    pl = [[r[2], loc, w, r[0], r[1]] for (loc, r), w in zip(refined, lwr)]
    pl.sort(key=lambda p: -p[2])            # make_ml_sorted (stable)
    return pl


def load_refpkg(refpkg_dir):
    """Read CONTENTS.json and return (tree text, reference alignment, phylo_model dict). The model is
    checked when it is built (models.model_from_refpkg)."""
    cpath = os.path.join(refpkg_dir, 'CONTENTS.json')
    if not os.path.exists(cpath):
        raise RefpkgMalformed(f'No CONTENTS.json in reference package: {refpkg_dir}')
    files = json.load(open(cpath)).get('files', {})
    for key in ('tree', 'aln_fasta', 'phylo_model'):
        if key not in files or not os.path.exists(os.path.join(refpkg_dir, files[key])):
            raise RefpkgMalformed(f'Reference package {refpkg_dir} has no "{key}" file.')
    tree_txt = open(os.path.join(refpkg_dir, files['tree'])).read().strip()
    ref = core.read_fasta(os.path.join(refpkg_dir, files['aln_fasta']))
    pm = json.load(open(os.path.join(refpkg_dir, files['phylo_model'])))
    return tree_txt, ref, pm


class Placer(object):
    """Place aligned queries on one reference package, exactly as pplacer does (CAT models; gamma
    models use gamma.GammaPlacer, see factory.new_placer)."""

    def __init__(self, refpkg_dir, queries, cpus=1):
        self.logger = logging.getLogger('timestamp')
        cpus = resources.parse_cpus(cpus)
        if cpus == 'auto':
            self.cpus = resources.available_cpus()          # threads for the cache
            tuned = resources.load_tuned_cpus()
            self.workers = min(tuned, self.cpus) if tuned else None   # None: calibrate in the prefilter
            self.cpus_source = 'saved benchmark' if tuned else 'calibration'
        else:
            self.cpus = self.workers = max(1, int(cpus))
            self.cpus_source = 'user'
        tree_txt, ref, pm = load_refpkg(refpkg_dir)
        self.alpha = models.alphabet_of(pm)                  # protein or DNA
        self.ns = self.alpha.ns
        self.lib = core.lib(self.ns)
        self.tree = parse_newick(tree_txt)
        self.children = self.tree.children
        self.bl = self.tree.branch_lengths_for_placement()
        self.label = self.tree.label
        self.root = self.tree.root
        self.nnodes = self.tree.n_nodes

        # pplacer_run.partition_queries: reference sequences inside the query file
        leaf_names = [self.label[v] for v in self.tree.leaves()]
        ref_names = set(leaf_names)
        if len(ref_names) != len(leaf_names):
            raise RefpkgMalformed('Duplicate leaf names in the reference tree.')
        in_ref = [q for q in queries if q[0] in ref_names]
        queries = [q for q in queries if q[0] not in ref_names]
        if in_ref:
            if len(in_ref) != len(ref_names):
                raise AlignmentMalformed('Some, but not all reference sequences found in provided '
                                         'alignment file. Are there overlaps between query and '
                                         'reference sequence IDs?')
            self.logger.info('Found reference sequences in the query alignment; using those as '
                             'the reference alignment (as pplacer does).')
            ref = in_ref
        if not queries:
            raise AlignmentMalformed('No query sequences to place.')
        names = [n for n, _ in queries]
        if len(set(names)) != len(names):
            raise AlignmentMalformed('Duplicate query names; every query needs a unique name.')

        # model (Refpkg.get_aln_fasta uppercases the alignment; the model is built before pre-masking)
        ref = [(n, s.upper()) for n, s in ref]
        self.model = models.model_from_refpkg(pm, ref, self.lib)
        self.statd = np.array(self.model.statd)
        self.x, self.l, self.xit, _ = core.build_diagd(self.model.R, list(self.statd), self.lib)
        self.rates = np.array(self.model.rates)
        self.nrates = len(self.rates)
        self.width = self.ns                                 # doubles per site (gamma: rates x states)
        if self.model.ras == 'cat':
            cats_full = self.model.site_categories
        else:                                                # gamma: every site under every rate
            cats_full = np.zeros(len(ref[0][1]), dtype=np.int32)

        # alignments and pplacer's pre-mask (informative in >=1 reference AND >=1 query)
        queries = [(n, s.upper()) for n, s in queries]
        L = len(ref[0][1])
        for nm, sq in ref + queries:
            if len(sq) != L:
                raise AlignmentMalformed(f'Sequence {nm} has length {len(sq)}, expected {L} '
                                         f'(the reference alignment length).')
        if len(cats_full) != L:
            raise RefpkgMalformed(f'The CAT site categories ({len(cats_full)}) do not match the '
                                  f'alignment length ({L}).')
        RA = np.frombuffer(''.join(s for _, s in ref).encode(), np.uint8).reshape(len(ref), L)
        QA = np.frombuffer(''.join(s for _, s in queries).encode(), np.uint8).reshape(len(queries), L)
        for A, seqs in ((RA, ref), (QA, queries)):
            ok = self.alpha.valid[A].all(axis=1)
            if not ok.all():
                i = int(np.nonzero(~ok)[0][0])
                raise AlignmentMalformed(f'Unknown characters {set(seqs[i][1].upper()) - set(self.alpha.char_map)} '
                                         f'in {seqs[i][0]}.')
        mask = core.informative(RA).any(axis=0) & core.informative(QA).any(axis=0)
        keep = np.nonzero(mask)[0]
        self.nsites = len(keep)
        self.cats = np.ascontiguousarray(cats_full[keep])
        self.occupied = np.zeros(self.nrates, np.uint8)
        self.occupied[np.unique(self.cats)] = 1
        self.refA = np.ascontiguousarray(RA[:, keep])
        self.QM = np.ascontiguousarray(QA[:, keep])
        self.query_names = names
        rowof = {n: i for i, (n, _) in enumerate(ref)}
        missing = [nm for nm in leaf_names if nm not in rowof]
        if missing:
            raise RefpkgMalformed(f'{len(missing)} tree leaves have no sequence in the reference '
                                  f'alignment, e.g. {missing[0]}.')
        self.leaf_row = {v: rowof[self.label[v]] for v in range(self.nnodes) if not self.children[v]}

        self.parent = [-1] * self.nnodes
        for v, ch in enumerate(self.children):
            for c in ch:
                self.parent[c] = v
        self.internal = [v for v in range(self.nnodes) if self.children[v]]
        self.iidx = {v: k for k, v in enumerate(self.internal)}
        self.local = threading.local()
        self.aln_len = L
        self.n_leaves = len(leaf_names)

    def memory_estimate(self, workers=None):
        """Estimated peak memory (bytes) of this pass (see resources.estimate_memory)."""
        w = workers or self.workers or resources.heuristic_cpus(self.cpus)
        return resources.estimate_memory(self.nnodes, len(self.internal), self.nsites,
                                         len(self.query_names), self.n_leaves, self.aln_len, w,
                                         width=self.width)

    # ------------------------------------------------------------------ helpers --
    def tl(self):
        """Per-thread (and per-process) scratch buffers."""
        t = self.local
        if not hasattr(t, 'tensor'):
            S = self.nsites
            t.tensor = np.zeros((self.nrates, self.ns, self.ns))
            t.b = [np.empty((S, self.width)) for _ in range(6)]
            t.e = [np.zeros(S, np.int64) for _ in range(4)]
        return t

    def evolve_t(self, dst, src, bl, tensor):
        self.lib.xp_multi_exp(P(tensor), P(self.x), P(self.l), P(self.xit), P(self.rates),
                              self.nrates, P(self.occupied), bl)
        self.lib.xp_evolve(P(dst), P(src), self.nsites, P(self.cats), P(tensor))

    def leaf_lv(self, v):
        return np.ascontiguousarray(self.alpha.lut[self.refA[self.leaf_row[v]]])

    def d_pre(self, v):
        """Distal vector before the final (-10) pull: (array, exponents)."""
        if self.children[v]:
            k = self.iidx[v]
            return self.darr[k], self.dexp[k]
        return self.leaf_lv(v), self.zero_e

    def _prox_into(self, dst, de, evp, evpe, rest, tmp, ten):
        """Like_stree.calc_proximal for one child: evp .* prod(evolved sibling distals), pull(-50)."""
        S = self.nsites
        if not rest:
            dst[:] = evp; de[:] = evpe
        else:
            a, e = self.d_pre(rest[0]); self.evolve_t(tmp, a, self.bl[rest[0]], ten)
            self.lib.xp_pairwise(P(dst), P(de), P(evp), P(evpe), P(tmp), P(e), S)
            for r in rest[1:]:
                a, e = self.d_pre(r); self.evolve_t(tmp, a, self.bl[r], ten)
                self.lib.xp_pairwise(P(dst), P(de), P(dst), P(de), P(tmp), P(e), S)
        self.lib.xp_pull_exponent(P(dst), P(de), S, -50)

    def get_d(self, v, out, oute):
        """Distal vector of node v as pplacer holds it (after the -10 pull)."""
        a, e = self.d_pre(v)
        out[:] = a; oute[:] = e
        self.lib.xp_pull_exponent(P(out), P(oute), self.nsites, -10)

    def get_p(self, v, out, oute):
        """Proximal vector of node v as pplacer holds it; leaf proximals are rebuilt exactly."""
        t = self.tl()
        if self.children[v]:
            k = self.iidx[v]
            out[:] = self.parr[k]; oute[:] = self.pexp[k]
        else:
            pv = self.parent[v]; kp = self.iidx[pv]
            self.evolve_t(t.b[4], self.parr[kp], self.bl[pv], t.tensor)
            self._prox_into(out, oute, t.b[4], self.pexp[kp],
                            [r for r in self.children[pv] if r != v], t.b[5], t.tensor)
        self.lib.xp_pull_exponent(P(out), P(oute), self.nsites, -10)

    # -------------------------------------------------------------------- cache --
    def _distal(self, v, t):
        S, ch, k = self.nsites, self.children[v], self.iidx[v]
        tA, tB, ten = t.b[0], t.b[1], t.tensor
        a, e = self.d_pre(ch[0]); self.evolve_t(tA, a, self.bl[ch[0]], ten)
        if len(ch) == 1:
            self.darr[k] = tA; self.dexp[k] = e
        else:
            b, f = self.d_pre(ch[1]); self.evolve_t(tB, b, self.bl[ch[1]], ten)
            self.lib.xp_pairwise(P(self.darr[k]), P(self.dexp[k]), P(tA), P(e), P(tB), P(f), S)
            for c in ch[2:]:
                a2, e2 = self.d_pre(c); self.evolve_t(tA, a2, self.bl[c], ten)
                self.lib.xp_pairwise(P(self.darr[k]), P(self.dexp[k]), P(self.darr[k]),
                                     P(self.dexp[k]), P(tA), P(e2), S)
        self.lib.xp_pull_exponent(P(self.darr[k]), P(self.dexp[k]), S, -50)

    def _proximal(self, v, t):
        kv, ch = self.iidx[v], self.children[v]
        evp = t.b[2]
        self.evolve_t(evp, self.parr[kv], self.bl[v], t.tensor)
        kids = []
        for c in ch:
            if not self.children[c]:
                continue                        # leaf proximal: rebuilt on demand
            kc = self.iidx[c]
            self._prox_into(self.parr[kc], self.pexp[kc], evp, self.pexp[kv],
                            [r for r in ch if r != c], t.b[3], t.tensor)
            kids.append(c)
        return kids

    def cache(self):
        """Like_stree.calc_distal_and_proximal for internal nodes (pre -10 pull)."""
        S, n = self.nsites, len(self.internal)
        self.zero_e = np.zeros(S, np.int64)
        self.darr = shared_zeros((n, S, self.width)); self.dexp = shared_zeros((n, S), np.int64)
        self.parr = shared_zeros((n, S, self.width)); self.pexp = shared_zeros((n, S), np.int64)
        lock = threading.Lock(); done = threading.Event(); errs = []; count = [0]
        ex = ThreadPoolExecutor(self.cpus)
        pending = [0] * self.nnodes
        for v in self.internal:
            pending[v] = sum(1 for c in self.children[v] if self.children[c])

        def distal_task(v):
            try:
                self._distal(v, self.tl())
                p = self.parent[v]; ready = False
                with lock:
                    count[0] += 1; fin = count[0] == n
                    if p >= 0:
                        pending[p] -= 1; ready = pending[p] == 0
                if ready:
                    ex.submit(distal_task, p)
                if fin:
                    done.set()
            except BaseException as e:
                errs.append(e); done.set()

        for v in self.internal:
            if pending[v] == 0:
                ex.submit(distal_task, v)
        done.wait()
        if errs:
            ex.shutdown(); raise errs[0]

        kr = self.iidx[self.root]
        self.parr[kr] = 1.0; self.pexp[kr] = 0
        count[0] = 0; done.clear()

        def prox_task(v):
            try:
                kids = self._proximal(v, self.tl())
                with lock:
                    count[0] += 1; fin = count[0] == n
                for c in kids:
                    ex.submit(prox_task, c)
                if fin:
                    done.set()
            except BaseException as e:
                errs.append(e); done.set()

        ex.submit(prox_task, self.root)
        done.wait()
        ex.shutdown()
        if errs:
            raise errs[0]

    # ---------------------------------------------------------------- prefilter --
    def _query_tables(self):
        """Prefilter inputs: every query character evolved by START_PEND under every rate
        (qevtab), and the queries as character codes."""
        ns, alpha = self.ns, self.alpha
        ten = np.zeros((self.nrates, ns, ns))
        self.lib.xp_multi_exp(P(ten), P(self.x), P(self.l), P(self.xit), P(self.rates),
                              self.nrates, P(self.occupied), core.START_PEND)
        nc = len(alpha.codes)
        self.qevtab = np.zeros((self.nrates, nc, ns))
        one = np.zeros((1, ns)); outv = np.empty((1, ns))
        for c in range(self.nrates):
            if not self.occupied[c]:
                continue
            cc = np.array([c], np.int32)
            for k, ch in enumerate(alpha.codes):
                one[0] = alpha.char_map[ch]
                self.lib.xp_evolve(P(outv), P(one), 1, P(cc), P(ten))
                self.qevtab[c, k] = outv[0]
        self.codes = np.ascontiguousarray(alpha.code_of[self.QM])

    def prefilter(self):
        """Score every query against every edge (pplacer's h_ranking scores), 8 edges per call."""
        self._query_tables()
        Qn = len(self.query_names)
        self.locs = np.array([v for v in range(self.nnodes) if v != self.root])
        E = len(self.locs)
        self.scores = shared_zeros((Qn, E))
        groups = [(i, min(E, i + 8)) for i in range(0, E, 8)]
        start = 0
        if self.workers is None:                    # --cpus auto without a saved benchmark
            start = self._calibrate(groups)
        self._run_groups(groups[start:], self.workers)

    def prefilter_seconds(self, n_queries, workers):
        """Benchmark helper: time a full prefilter sweep restricted to the first n_queries
        queries (the per-edge work is fixed, the kernel work grows with the queries)."""
        names, QM = self.query_names, self.QM
        self.query_names, self.QM = names[:n_queries], np.ascontiguousarray(QM[:n_queries])
        saved = self.workers
        try:
            self.workers = workers
            t0 = time.time()
            self.prefilter()
            return time.time() - t0
        finally:
            self.query_names, self.QM, self.workers = names, QM, saved

    def _run_groups(self, groups, workers, chunksize=2):
        """Run prefilter groups in a pool of `workers` processes; returns the elapsed time
        (pool start-up excluded)."""
        global _PLACER
        if not groups:
            return 0.0
        _PLACER = self
        try:
            with _fork_pool(workers) as pool:
                t0 = time.time()
                for _ in pool.imap_unordered(_prefilter_group, groups, chunksize=chunksize):
                    pass
                return time.time() - t0
        finally:
            _PLACER = None
            gc.unfreeze()

    def _calibrate(self, groups, per_worker=10):
        """--cpus auto: time consecutive slices of the prefilter at several worker counts (the
        work is real, not repeated) and keep the fastest. Small jobs use the heuristic."""
        n = self.cpus
        cands = [c for c in resources.default_candidates(n) if c >= max(1, n // 3)]
        need = per_worker * sum(cands)
        if len(groups) < 2 * need or len(cands) < 2:
            self.workers = resources.heuristic_cpus(n)
            self.cpus_source = 'heuristic'
            self.logger.info(f'--cpus auto: job too small to calibrate; using {self.workers} of {n} CPUs.')
            return 0
        start, rates = 0, {}
        for c in cands:
            chunk = groups[start:start + per_worker * c]
            dt = self._run_groups(chunk, c, chunksize=1)
            rates[c] = len(chunk) / max(dt, 1e-9)
            start += len(chunk)
        top = max(rates.values())                   # within 3% of the fastest: take fewer CPUs
        self.workers = min(c for c in cands if rates[c] >= 0.97 * top)
        self.cpus_source = 'calibration'
        best = rates[self.workers]
        self.logger.info('--cpus auto calibration (relative speed): ' +
                         ', '.join(f'{c}: {rates[c] / best:.2f}' for c in cands) +
                         f' -> using {self.workers} of {n} CPUs.')
        return start

    def _supernodes8(self, i0, i1):
        """Supernode vectors of the edges locs[i0:i1] (at most 8) as pplacer's prefilter builds
        them (distal and proximal evolved by half the branch length, times the frequencies);
        returns the thread's scratch buffers holding them."""
        S = self.nsites
        t = self.tl()
        if not hasattr(t, 'sn8'):
            t.sn8 = np.empty((8, S, self.width)); t.se8 = np.zeros((8, S), np.int64)
        for lane, k in enumerate(range(i0, i1)):
            v = int(self.locs[k])
            self.get_d(v, t.b[0], t.e[0]); self.get_p(v, t.b[1], t.e[1])
            self.evolve_t(t.b[2], t.b[0], self.bl[v] / 2., t.tensor)
            self.evolve_t(t.b[3], t.b[1], self.bl[v] / 2., t.tensor)
            self.lib.xp_statd_pairwise(P(t.sn8[lane]), P(t.se8[lane]), P(self.statd), P(t.b[2]),
                                       P(t.e[0]), P(t.b[3]), P(t.e[1]), S)
        return t

    def _prefilter_group(self, g):
        i0, i1 = g
        S, Qn, E, nc = self.nsites, len(self.query_names), self.scores.shape[1], len(self.alpha.codes)
        t = self._supernodes8(i0, i1)
        st = self.lib.xp_edge_scores8(P(self.scores), E, i0, i1 - i0, P(t.sn8), P(t.se8),
                                      P(self.qevtab), P(self.cats), S, nc, P(self.codes), Qn)
        if st != 0:
            raise GPlaceException(f'Prefilter kernel failed ({st}) for edges {i0}..{i1}.')

    def ranking(self, qi, K=96):
        """h_ranking: edges by prefilter score (descending), ties by ascending edge id."""
        row = self.scores[qi]
        E = len(row)
        if K < E:
            kth = np.partition(row, E - K)[E - K]
            cand = np.nonzero(row >= kth)[0]
        else:
            cand = np.arange(E)
        for k in sorted(cand, key=lambda k: (-row[k], self.locs[k])):
            yield int(self.locs[k])
        if len(cand) < E:                                  # rare: exhausted candidates
            rest = np.setdiff1d(np.arange(E), cand)
            for k in sorted(rest, key=lambda k: (-row[k], self.locs[k])):
                yield int(self.locs[k])

    # ---------------------------------------------------------------- placement --
    def _query_evaluator(self, qi):
        """(evaluate, free) for query qi. evaluate(loc, pend, dist, tol) sets the three-taxon tree
        on edge loc with starting branch lengths pend/dist and optimises it as pplacer's
        safe_ml_optimize_location does (after a GSL error: retry from the default lengths);
        it returns ((like, pend, dist), None), or (None, status) if the location failed."""
        L = self.lib
        t = self.tl()
        name = self.query_names[qi]
        qrow = self.QM[qi]
        reind = np.nonzero(core.informative(qrow))[0].astype(np.int32)
        m = len(reind)
        qm = np.ascontiguousarray(self.alpha.lut[qrow][reind])
        catm = np.ascontiguousarray(self.cats[reind])
        tt = L.xp_tt_new(m, P(catm), P(self.statd), P(self.x), P(self.l), P(self.xit),
                         P(self.rates), self.nrates, P(qm))
        if not tt:
            raise GPlaceException(f'Could not set up the likelihood for {name}.')
        res = np.empty(4)
        D, De, Pp, Pe = t.b[0], t.e[0], t.b[1], t.e[1]
        cur = [None]

        def set_edges(loc, pend, dist):
            if cur[0] != loc:
                self.get_d(loc, D, De); self.get_p(loc, Pp, Pe); cur[0] = loc
            L.xp_tt_set(tt, P(D), P(De), P(Pp), P(Pe), P(reind), self.bl[loc], pend, dist)

        def optimize(tol):
            st = L.xp_tt_optimize(tt, tol, core.MAX_PEND, core.MAX_ITER)
            if st in (core.XP_OK, core.XP_EXCEEDED):
                L.xp_tt_results(tt, P(res))
                return (float(res[0]), float(res[1]), float(res[2])), None
            return None, st

        def evaluate(loc, pend, dist, tol):     # safe_ml_optimize_location
            set_edges(loc, pend, dist)
            r, st = optimize(tol)
            if st == core.XP_GSL:
                set_edges(loc, core.START_PEND, self.bl[loc] / 2.)
                r, st = optimize(tol)
            return r, st

        arrays = [qm, catm]                     # the engine holds pointers into these until freed

        def free():
            L.xp_tt_free(tt)
            arrays.clear()

        return evaluate, free

    def place_query(self, qi):
        """Baseball + three-taxon optimisation + keep/refine + LWR (core.ml)."""
        name = self.query_names[qi]
        evaluate, free = self._query_evaluator(qi)
        try:
            bb = Baseball()
            for loc in self.ranking(qi):
                r, st = evaluate(loc, core.START_PEND, self.bl[loc] / 2., core.INITIAL_TOL)
                if not bb.add(loc, r, st):
                    break
            if not bb.results:
                raise GPlaceException(f'empty results for {name}!')
            refined = []
            for loc, (lk, pend, dist) in bb.kept(self.lib):
                r, _ = evaluate(loc, pend, dist, core.FINAL_TOL)
                refined.append((loc, r if r is not None else (lk, pend, dist)))
            pl = weighted_placements(self.lib, refined)
        finally:
            free()
        return name, pl, bb.warns

    def place_all(self, workers=None):
        """Place every query; returns [(name, placements, warnings)] in input order.
        placements: [distal_length, edge_num, like_weight_ratio, likelihood, pendant_length]."""
        global _PLACER
        workers = workers or self.workers or resources.heuristic_cpus(self.cpus)
        n = len(self.query_names)
        results = [None] * n
        _PLACER = self
        chunks = [range(i, min(n, i + 4)) for i in range(0, n, 4)]
        try:
            with _fork_pool(workers) as pool:
                for out in pool.imap_unordered(_place_chunk, chunks):
                    for qi, nm, pl, w in out:
                        results[qi] = (nm, pl, w)
        finally:
            _PLACER = None
            gc.unfreeze()
        return results

    def release(self):
        """Free the likelihood cache and prefilter scores (the shared-memory arrays)."""
        for attr in ('darr', 'dexp', 'parr', 'pexp', 'scores', 'codes', 'qevtab'):
            if hasattr(self, attr):
                delattr(self, attr)
        self.local = threading.local()
        gc.collect()

    def jplace_tree(self):
        """The reference tree with edge numbers, formatted as pplacer writes it."""
        return to_newick(self.tree, with_node_numbers=True)

    def run(self):
        """cache -> prefilter -> placement."""
        self.cache()
        self.logger.info(f'Cached likelihood vectors for {len(self.internal):,} internal nodes '
                         f'({(self.darr.nbytes + self.parr.nbytes) / 1e9:.1f} GB).')
        self.prefilter()
        self.logger.info(f'Scored {len(self.query_names):,} queries against {len(self.locs):,} edges.')
        results = self.place_all()
        self.release()
        return results
