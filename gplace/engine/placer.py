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
(3 x nodes x sites x 20 doubles, ~104 GB for the largest GTDB R226 class-level tree).
gPlace stores distal/proximal vectors for internal nodes only, in their state before pplacer's
final "pull exponents (-10)" pass, and re-applies that pass to a copy at use time; power-of-two
rescaling is exact, so every value is bit-identical. Leaf vectors are rebuilt on demand in
pplacer's operation order and supernodes are computed during the prefilter sweep.

Speed: the likelihood-vector cache is built by a dependency-driven thread pool; the prefilter
and the placements run in forked worker processes, with the large arrays in shared memory.
"""

import ctypes
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
from gplace.engine.reference import Reference
from gplace.exceptions import RefpkgMalformed, AlignmentMalformed, GPlaceException, GPlaceExit


def shared_zeros(shape, dtype=np.float64):
    """Zero-filled array in shared anonymous memory: visible to forked workers, never copied
    on fork, charged once to a memory cgroup."""
    count = int(np.prod(shape))
    buf = mmap.mmap(-1, max(count * np.dtype(dtype).itemsize, 1),
                    flags=mmap.MAP_SHARED | mmap.MAP_ANONYMOUS)
    return np.frombuffer(buf, dtype=dtype, count=count).reshape(shape)


_PLACER = None                                  # the Placer seen by forked workers
TENSOR_BYTES = 128e6                            # rebuild: transition matrices held at a time
BLOCK_SITES = 160                               # rebuild: fewest sites per thread


def _free_memory_call():
    """The C library's call that gives freed malloc memory back to the system, or None."""
    try:
        libc = ctypes.CDLL(None)
    except OSError:
        return None
    if hasattr(libc, 'malloc_trim'):                        # glibc
        return lambda: libc.malloc_trim(0)
    if hasattr(libc, 'malloc_zone_pressure_relief'):        # macOS: every zone (NULL), as much as possible (0)
        f = libc.malloc_zone_pressure_relief
        f.argtypes, f.restype = [ctypes.c_void_p, ctypes.c_size_t], ctypes.c_size_t
        return lambda: f(None, 0)
    return None


_free_memory = _free_memory_call()


def _fork_pool(n):
    gc.collect()
    gc.freeze()                                 # keep object pages shared after fork
    if _free_memory is not None:
        _free_memory()                          # freed memory back to the system: not shared, so not copied
    return multiprocessing.get_context('fork').Pool(n)


def _batch_prefilter_group(g):
    _PLACER._batch_prefilter_group(g)


def _place_chunk(qis):
    return [(qi,) + tuple(_PLACER.place_query(qi)) for qi in qis]


class Baseball(object):
    """pplacer's stopping rule for the first-stage optimisations ("baseball"), fed one location at
    a time in ranked order: stop after max_strikes locations more than strike_box log-units below
    the best so far, or after max_pitches locations (max_strikes 0: never; see core.Search)."""

    def __init__(self, search):
        self.search = search
        self.results, self.warns = [], []
        self.like_record, self.n_strikes = -math.inf, 0

    def add(self, loc, r, st):
        """Record one location's result (r = (like, pend, dist), or None with GSL status st);
        returns False when pplacer tries no further locations."""
        if r is None:
            self.warns.append((loc, st))
            return True
        like = r[0]; prev = len(self.results); self.results.append((loc, r))
        if prev >= self.search.t_max_pitches:
            return False
        elif like > self.like_record:
            self.like_record = like
        elif like < self.like_record - self.search.strike_box:
            if self.n_strikes + 1 >= self.search.t_max_strikes:
                return False
            self.n_strikes += 1
        return True

    def kept(self, lib):
        """Locations refined at the final tolerance: up to keep_at_most within log(keep_factor)
        of the best, best first."""
        srt = sorted(self.results, key=lambda x: -x[1][0])
        best = srt[0][1][0]
        lkf = lib.xp_log(self.search.keep_factor)
        return [x for i, x in enumerate(srt) if i < self.search.keep_at_most and x[1][0] >= lkf + best]


_X86_NAN = math.copysign(math.nan, -1.0)     # x86-64's default NaN: its sign bit is set


def _x86_sub(a, b):
    """a - b as on x86-64, where pplacer's release binary runs: inf - inf gives the default NaN, which has
    its sign bit set, so pplacer writes -nan. ARM64's default NaN is positive (nan)."""
    return _X86_NAN if a == b and math.isinf(a) else a - b


def ll_normalized_prob(lib, lls):
    """pplacer's ll_normalized_prob: 1 / sum(exp(other - this)) for each log likelihood."""
    return [1. / sum(lib.xp_exp(_x86_sub(o, ll)) for o in lls) for ll in lls]


def weighted_placements(lib, refined, marginals=None):
    """Like-weight ratios and pplacer's placement order for [(loc, (like, pend, dist))];
    returns [[distal_length, edge_num, like_weight_ratio, likelihood, pendant_length]].
    With marginals (the marginal likelihood of each placement, -p): rows
    [distal_length, edge_num, like_weight_ratio, likelihood, marginal_like, pendant_length, post_prob].
    A NaN keeps its sign through exp, the sum and the division, on x86-64 and ARM64 alike."""
    lwr = ll_normalized_prob(lib, [r[0] for _, r in refined])
    if marginals is None:
        pl = [[r[2], loc, w, r[0], r[1]] for (loc, r), w in zip(refined, lwr)]
    else:
        pp = ll_normalized_prob(lib, marginals)
        pl = [[r[2], loc, w, r[0], m, r[1], q] for (loc, r), w, m, q in zip(refined, lwr, marginals, pp)]
    pl.sort(key=lambda p: -p[2])            # make_ml_sorted (stable)
    return pl


def tree_length(children, bl, root):
    """Gtree.tree_length: every branch length but the root's, each node's own length first and
    then its subtrees' totals, left to right."""
    tot = [0.] * len(children)
    order, stack = [], [root]
    while stack:
        v = stack.pop(); order.append(v); stack.extend(children[v])
    for v in reversed(order):                   # children before parents
        acc = 0. if v == root else bl[v]
        for c in children[v]:
            acc = acc + tot[c]
        tot[v] = acc
    return tot[root]


def prior_mean_map(children, bl, root, prior_lower):
    """Core.prior_mean_map ((+.) prior_lower): for each edge, prior_lower + its length plus the
    average distance from its lower end to the leaves below (pplacer's --inform-prior means),
    summed in pplacer's order (leaf distances in subtree order, a left fold)."""
    means, dist = {}, {}
    order, stack = [], [root]
    while stack:
        v = stack.pop(); order.append(v); stack.extend(children[v])
    for v in reversed(order):
        top = prior_lower + bl[v]
        if not children[v]:
            means[v] = top
            dist[v] = np.array([bl[v]])
            continue
        d = np.concatenate([dist.pop(c) for c in children[v]])
        means[v] = float(np.add.accumulate(top + d)[-1]) / float(len(d))
        dist[v] = bl[v] + d
    del means[root]
    return means


def pp_prior(tree, bl, search):
    """The prior on the pendant length for -p (pplacer_run.ml): (uniform, mean) where mean is
    a float (Flat_exp_prior: tree length / number of edges) or {edge: mean} (--inform-prior)."""
    if search.uniform_prior:
        return True, 0.
    if search.informative_prior:
        return False, prior_mean_map(tree.children, bl, tree.root, search.prior_lower)
    return False, tree_length(tree.children, bl, tree.root) / float(tree.n_nodes - 1)


_PP_ERRORS = {core.XP_EXCEEDED: 'Minimization.ExceededMaxIter', core.XP_FINDSTART: 'Minimization.FindStartFailure',
              core.XP_INVALIDSTART: 'Minimization.InvalidStartValues'}


_NO_EDGE = np.iinfo(np.int64).max       # filler in the candidate lists: ranks after every edge


def candidates_needed(search, n_edges, n_cand=resources.N_CAND):
    """Ranked edges kept per query: 4 per edge pplacer's search can use (max_pitches + 1: 164 with
    the default 40), at least n_cand; every edge when --max_strikes 0 lifts the limit."""
    tp = search.t_max_pitches
    if tp == math.inf:
        return n_edges
    return min(max(int(n_cand), 4 * (tp + 1)), n_edges)


def merge_candidates(cs, ce, b0, edges, sc):
    """Merge the scores sc of the edges `edges` (rows b0.. of cs/ce) into the candidate lists,
    keeping each query's first K edges in ranking order. The selection is by score; when equal
    scores straddle the cut-off, that row is resolved by (score, edge number) exactly."""
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


def _distinct_rows(A):
    """(row -> index of its distinct row, first row of every distinct row), distinct rows numbered
    in order of first occurrence. Rows are hashed (BLAKE2b), and equal hashes are verified."""
    import hashlib
    seen, first = {}, []
    rep_of = np.empty(len(A), np.int64)
    for i in range(len(A)):
        key = hashlib.blake2b(A[i].tobytes(), digest_size=16).digest()
        for r in seen.get(key, ()):
            if np.array_equal(A[first[r]], A[i]):
                break
        else:
            r = len(first)
            first.append(i)
            seen.setdefault(key, []).append(r)
        rep_of[i] = r
    return rep_of, first


class Placer(object):
    """Place aligned queries on one reference (a reference package directory or a
    reference.Reference), exactly as pplacer does (CAT models; gamma models use gamma.GammaPlacer,
    see factory.new_placer)."""

    def __init__(self, refpkg_dir, queries, cpus=1, discard_nonoverlapped=False, search=None,
                 always_refine=False, ranked_edges=None, max_memory=None):
        self.logger = logging.getLogger('timestamp')
        self.search = search = search or core.Search()
        if ranked_edges is not None and int(ranked_edges) < 1:
            raise GPlaceException('--ranked_edges must be at least 1.')
        self.ranked_edges = None if ranked_edges is None else int(ranked_edges)
        self.rank = None
        cpus = resources.parse_cpus(cpus)
        if cpus == 'auto':
            self.cpus = resources.available_cpus()          # threads for the cache
            tuned = resources.load_tuned_cpus()
            self.workers = min(tuned, self.cpus) if tuned else None   # None: calibrate in the prefilter
            self.cpus_source = 'saved benchmark' if tuned else 'calibration'
        else:
            self.cpus = self.workers = max(1, int(cpus))
            self.cpus_source = 'user'
        self.reference = Reference.of(refpkg_dir)
        tree_txt, pm = self.reference.tree_text(), self.reference.phylo_model()
        self.alpha = models.alphabet_of(pm)                  # protein or DNA
        self.ns = self.alpha.ns
        self.lib = core.lib(self.ns)
        self.tree = parse_newick(tree_txt)
        del tree_txt
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
        # the alignments are core.Alignment matrices (uppercased), never lists of strings
        queries = core.Alignment.of(queries)
        in_ref = [i for i, n in enumerate(queries.names) if n in ref_names]
        self.from_input_alignment = bool(in_ref)
        if in_ref:
            if len(in_ref) != len(ref_names):
                raise AlignmentMalformed('Some, but not all reference sequences found in provided '
                                         'alignment file. Are there overlaps between query and '
                                         'reference sequence IDs?')
            seen = set()
            for i in in_ref:
                nm = queries.names[i]
                if nm in seen:
                    raise AlignmentMalformed(f"duplicate reference sequence '{nm}' in the query file")
                seen.add(nm)
            self.logger.info('Found reference sequences in the query alignment; using those as '
                             'the reference alignment (as pplacer does).')
            ref = queries.take(in_ref)
            queries = queries.take([i for i, n in enumerate(queries.names) if n not in ref_names])
        else:
            if max_memory is not None:
                self._check_loading_memory(max_memory, len(leaf_names), queries, pm)
            ref = self.reference.alignment_matrix()
            if ref is None:
                raise RefpkgMalformed('Please specify a reference alignment with --ref_aln or --refpkg, or '
                                      'include all reference sequences in the query alignment.')
        if not len(queries):
            raise AlignmentMalformed('No query sequences to place.')
        names = list(queries.names)
        if len(set(names)) != len(names):
            raise AlignmentMalformed('Duplicate query names; every query needs a unique name.')

        # model (Refpkg.get_aln_fasta uppercases the alignment; the model is built before pre-masking)
        self.model = models.model_from_refpkg(pm, ref, self.lib)
        self.statd = np.array(self.model.statd)
        self.x, self.l, self.xit, _ = core.build_diagd(self.model.R, list(self.statd), self.lib)
        self.rates = np.array(self.model.rates)
        self.pp_prior = pp_prior(self.tree, self.bl, search) if search.calc_pp else None
        self.nrates = len(self.rates)
        self.width = self.ns                                 # doubles per site (gamma: rates x states)
        # pplacer refines the model when the reference sequences come from the query file or with
        # --always-refine; only the CAT refinement does something (Gmix_model.refine is a no-op)
        self.refined = self.model.ras == 'cat' and (self.from_input_alignment or always_refine)
        if self.model.ras == 'cat' and not self.refined:
            cats_full = self.model.site_categories
        else:                                                # gamma: every site under every rate
            cats_full = np.zeros(ref.length, dtype=np.int32)    # refined: set below

        # alignments and pplacer's pre-mask (informative in >=1 reference AND >=1 query)
        L = ref.length                                       # each alignment has one length (Alignment.of)
        if queries.length != L:
            raise AlignmentMalformed(f'Sequence {names[0]} has length {queries.length}, expected {L} '
                                     f'(the reference alignment length).')
        if len(cats_full) != L:
            raise RefpkgMalformed(f'The CAT site categories ({len(cats_full)}) do not match the '
                                  f'alignment length ({L}).')
        RA, QA = ref.matrix, queries.matrix
        for A, seq_names in ((RA, ref.names), (QA, names)):
            for i, B in core.row_blocks(A):
                ok = self.alpha.valid[B].all(axis=1)
                if not ok.all():
                    j = i + int(np.nonzero(~ok)[0][0])
                    bad = {chr(c) for c in np.unique(A[j])} - set(self.alpha.char_map)
                    raise AlignmentMalformed(f'Unknown characters {bad} in {seq_names[j]}.')
        # pplacer_run.premask: every query must have a residue in a column where some reference
        # has one ("doesn't overlap any reference sequences"), unless --discard-nonoverlapped.
        # --no-pre-mask: every column is kept and nothing is checked or discarded; a query without
        # a residue stops pplacer (core.ml)
        if search.no_pre_mask:
            for i, B in core.row_blocks(QA):
                empty = np.nonzero(~core.informative(B).any(axis=1))[0]
                if len(empty):
                    raise AlignmentMalformed(f"sequence '{names[i + empty[0]]}' has no informative sites.")
            ref_inf = overlaps = np.ones(L, bool)
        else:
            ref_inf = core.any_informative(RA)
            overlaps = np.concatenate([(core.informative(B) & ref_inf).any(axis=1) for _, B in core.row_blocks(QA)])
        self.discarded = [] if search.no_pre_mask else [names[i] for i in np.nonzero(~overlaps)[0]]
        if self.discarded:
            if not discard_nonoverlapped:
                raise AlignmentMalformed(f'Sequence {self.discarded[0]} doesn\'t overlap any reference '
                                         f'sequences.' + (f' ({len(self.discarded)} such sequences; '
                                                          f'--discard_nonoverlapped drops them.)'
                                                          if len(self.discarded) > 1 else
                                                          ' (--discard_nonoverlapped drops it.)'))
            keep_q = np.nonzero(overlaps)[0]
            names = [names[i] for i in keep_q]
            QA = np.ascontiguousarray(QA[keep_q])
            if not names:
                raise AlignmentMalformed('No query sequences to place: none overlaps the reference sequences.')
        mask = ref_inf if search.no_pre_mask else ref_inf & core.any_informative(QA)
        keep = np.nonzero(mask)[0]
        self.nsites = len(keep)
        self.cats = np.ascontiguousarray(cats_full[keep])
        self.occupied = np.zeros(self.nrates, np.uint8)
        self.occupied[np.unique(self.cats)] = 1
        self.refA = core.drop_columns(RA, keep)              # RA is ours (read here): masked in place
        ref.matrix = RA = None
        # pplacer places each distinct (pre-masked) sequence once: query_names/QM hold the first
        # query of every distinct sequence; all_names and rep_of map every query to its row
        QM = QA if len(keep) == L else np.take(QA, keep, axis=1)
        rep_of, first = _distinct_rows(QM)
        self.all_names, self.rep_of = names, rep_of
        self.QM = QM if len(first) == len(names) else np.ascontiguousarray(QM[first])
        self.query_names = [names[i] for i in first]
        self._name_index = None
        rowof = {n: i for i, n in enumerate(ref.names)}
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
        self._tensors = None                         # _rebuild's transition matrices (shared memory)
        self.aln_len = L
        self.n_leaves = len(leaf_names)
        if self.refined:
            self.logger.info('Refining the CAT site categories on the reference tree (pplacer refines the '
                             'model when ' + ('the reference sequences are in the query alignment).'
                                              if self.from_input_alignment else '--always_refine is given).'))
            self.refine_site_categories()

    def _check_loading_memory(self, max_memory, n_leaves, queries, pm):
        """Stop before reading the reference alignment if reading it alone would take more than
        max_memory: under a hard memory limit the process would be killed before the plan's
        message. The message is the plan's (tiled.plan_tiling) where the model allows planning
        here, with the sites informative in a query (at least the pre-masked sites)."""
        load = resources.current_memory() + n_leaves * queries.length + core.Alignment.BLOCK_BYTES
        if load <= max_memory:
            return
        width = models.width_of(pm)
        if width is not None:
            from gplace.engine.tiled import plan_tiling
            sites = queries.length if self.search.no_pre_mask else int(core.any_informative(queries.matrix).sum())
            plan_tiling(self.tree, sites, len(queries), n_leaves, queries.length,
                        self.workers or resources.heuristic_cpus(self.cpus), None, max_memory, width=width,
                        n_cand=self.ranked_edges)                    # raises when no plan fits
        GiB = resources.GiB
        raise GPlaceExit(f'--max_memory {max_memory / GiB:.1f} GiB is too small: reading the reference '
                         f'alignment ({n_leaves:,} sequences of {queries.length:,} columns) needs '
                         f'--max_memory {load / GiB:.1f}G or more.')

    # --------------------------------------------------------------- refinement --
    def _heavy_first_postorder(self):
        """Postorder with the largest subtree of every node first, so that a pass keeping only the
        vectors of finished subtrees holds O(log n) of them."""
        size = [1] * self.nnodes
        for v in self.tree.postorder():
            for c in self.children[v]:
                size[v] += size[c]
        order, stack = [], [(self.root, False)]
        while stack:
            v, seen = stack.pop()
            if seen or not self.children[v]:
                order.append(v)
                continue
            stack.append((v, True))
            for c in sorted(self.children[v], key=lambda c: size[c]):   # largest popped first
                stack.append((c, False))
        return order

    def _site_log_likes(self, cat, order):
        """Like_stree.site_log_like_arr with every site in category cat: the distal vectors of
        calc_distal_and_evolv_dist (the same operations as the cache), then
        site_log_like_arr3 at the root; returns the per-site log-likelihoods."""
        S, ns, R, L = self.nsites, self.ns, self.nrates, self.lib
        cats = np.full(S, cat, np.int32)
        occ = np.zeros(R, np.uint8); occ[cat] = 1
        ten = np.zeros((R, ns, ns))
        zero_e = np.zeros(S, np.int64)

        def evolve(src, bl):
            dst = np.empty((S, ns))
            L.xp_multi_exp(P(ten), P(self.x), P(self.l), P(self.xit), P(self.rates), R, P(occ), bl)
            L.xp_evolve(P(dst), P(src), S, P(cats), P(ten))
            return dst

        ed = {}                                  # finished node -> (evolved distal, exponents)
        for v in order:
            if v == self.root:                   # the root's own distal vector is not used
                break
            ch = self.children[v]
            if not ch:
                a, e = self.leaf_lv(v), zero_e
            else:
                a, e = ed.pop(ch[0])
                if len(ch) > 1:                  # Glv.listwise_prod
                    b, f = ed.pop(ch[1])
                    d, de = np.empty((S, ns)), np.empty(S, np.int64)
                    L.xp_pairwise(P(d), P(de), P(a), P(e), P(b), P(f), S)
                    for c in ch[2:]:
                        b, f = ed.pop(c)
                        L.xp_pairwise(P(d), P(de), P(d), P(de), P(b), P(f), S)
                    a, e = d, de
                else:
                    a, e = a.copy(), e.copy()
                L.xp_pull_exponent(P(a), P(e), S, -50)
            ed[v] = (evolve(a, self.bl[v]), e)
        top = self.children[self.root]
        if len(top) < 2:
            raise GPlaceException('site_log_like_arr: degree two node')
        (x, xe), (y, ye) = ed[top[0]], ed[top[1]]
        if len(top) == 2:
            z, ze = np.ones((S, ns)), zero_e     # Glv.set_unit
        elif len(top) == 3:
            z, ze = ed[top[2]]
        else:
            (a, e), (b, f) = ed[top[2]], ed[top[3]]
            z, ze = np.empty((S, ns)), np.empty(S, np.int64)
            L.xp_pairwise(P(z), P(ze), P(a), P(e), P(b), P(f), S)
            for c in top[4:]:
                b, f = ed[c]
                L.xp_pairwise(P(z), P(ze), P(z), P(ze), P(b), P(f), S)
        like = np.zeros(S)                       # site_log_like_arr3, state by state
        for st in range(ns):
            like = like + ((self.statd[st] * x[:, st]) * y[:, st]) * z[:, st]
        bad = np.nonzero(like <= 0.)[0]
        if len(bad):
            raise GPlaceException(f'Site {int(bad[0])} has zero likelihood.')
        log2 = L.xp_log(2.)
        twoexp = xe + ye + ze
        return np.array([L.xp_log(float(like[i])) + log2 * float(twoexp[i]) for i in range(S)])

    def refine_site_categories(self):
        """Gcat_model.refine: each site gets the category c with the largest
        site log-likelihood (all sites in c) + 2 log(rate c) - 3 rate c, ties to the lower c."""
        S, R = self.nsites, self.nrates
        order = self._heavy_first_postorder()
        with ThreadPoolExecutor(max(1, min(self.cpus, R))) as ex:
            lls = list(ex.map(lambda c: self._site_log_likes(c, order), range(R)))
        best, best_cat = np.full(S, -np.inf), np.full(S, -1, np.int32)
        for c in range(R):                       # record_best_log_lks in category order
            x = (lls[c] + 2. * self.lib.xp_log(float(self.rates[c]))) - 3. * float(self.rates[c])
            up = x > best
            best[up], best_cat[up] = x[up], c
        if (best_cat < 0).any():
            raise GPlaceException(f'No rate category for site {int(np.nonzero(best_cat < 0)[0][0])} '
                                  f'(its log-likelihoods are not numbers).')
        self.cats = np.ascontiguousarray(best_cat)
        self.occupied = np.zeros(R, np.uint8)
        self.occupied[np.unique(self.cats)] = 1

    def expand(self, results):
        """Results of the distinct sequences (run/place_all) -> one per query, in input order."""
        return [(nm,) + tuple(results[r][1:]) for nm, r in zip(self.all_names, self.rep_of)]

    def records(self, names):
        """jplace records for the queries `names` (one output file), as pplacer groups them:
        queries with the same pre-masked sequence share a record, names in input order. The mask is
        the file's own (columns with a residue in some reference and some query of the file), as
        when pplacer places that file alone (--no-pre-mask: the whole sequences). Returns
        [[name, ...], ...] in first-occurrence order."""
        if self._name_index is None:
            self._name_index = {n: i for i, n in enumerate(self.all_names)}
        ii = np.array([self._name_index[n] for n in names], np.int64)
        rows = self.QM[self.rep_of[ii]]
        if not self.search.no_pre_mask:
            rows = np.ascontiguousarray(rows[:, core.informative(rows).any(axis=0)])
        rep_of, first = _distinct_rows(rows)
        out = [[] for _ in first]
        for n, r in zip(names, rep_of):
            out[r].append(n)
        return out

    def memory_estimate(self, workers=None):
        """Estimated peak memory (bytes) of this pass (see resources.estimate_memory)."""
        w = workers or self.workers or resources.heuristic_cpus(self.cpus)
        return resources.estimate_memory(self.nnodes, len(self.internal), self.nsites,
                                         len(self.query_names), self.n_leaves, self.aln_len, w,
                                         width=self.width, n_cand=self.ranked_edges)

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
    def cache(self):
        """Like_stree.calc_distal_and_proximal for internal nodes (pre -10 pull), by the engine (_rebuild)."""
        S, n = self.nsites, len(self.internal)
        self.zero_e = np.zeros(S, np.int64)
        self.darr = shared_zeros((n, S, self.width)); self.dexp = shared_zeros((n, S), np.int64)
        self.parr = shared_zeros((n, S, self.width)); self.pexp = shared_zeros((n, S), np.int64)
        kr = self.iidx[self.root]
        self.parr[kr] = self.root_proximal(); self.pexp[kr] = 0
        with ThreadPoolExecutor(self.cpus) as ex:
            self._rebuild(self._schedule(self.root, set(self.internal)), True, ex)
        self._tensors = None

    def _rebuild_model(self):
        """What xp_rebuild needs of the model: rate categories per site (1), gamma (0), every site's
        matrix (its CAT category) and the eigen-system for xp_multi_exp_many (one, not per category)."""
        return 1, 0, self.cats, (self.x, self.l, self.xit, 0)

    def _schedule(self, top, internal):
        """The internal nodes below top (top included) for xp_rebuild: the nodes children first, their
        rows (slot, from self.iidx), their children in tree order (ptr, ch: a row, or -1 - the leaf's row
        in refA; a child with children outside `internal` is the root of a tile below, whose distal row
        is set by the caller) and chunks (first node, end, branch lengths of the chunk's child edges then
        of its nodes' own edges, number of child edges) of at most TENSOR_BYTES of matrices."""
        post, stack = [], [(top, False)]
        while stack:
            v, seen = stack.pop()
            if seen:
                post.append(v)
                continue
            stack.append((v, True))
            stack.extend((c, False) for c in reversed(self.children[v]) if c in internal)
        ptr, ch, ch_bl = [0], [], []
        for v in post:
            for c in self.children[v]:
                ch.append(self.iidx[c] if self.children[c] else -1 - self.leaf_row[c])
                ch_bl.append(self.bl[c])
            ptr.append(len(ch))
        cap = max(1, int(TENSOR_BYTES // (self.nrates * self.ns * self.ns * 8)))
        chunks, i0, n = [], 0, len(post)
        while i0 < n:
            i1 = i0 + 1
            while i1 < n and ptr[i1 + 1] - ptr[i0] + i1 + 1 - i0 <= cap:
                i1 += 1
            bls = np.array(ch_bl[ptr[i0]:ptr[i1]] + [self.bl[v] for v in post[i0:i1]], np.float64)
            chunks.append((i0, i1, bls, ptr[i1] - ptr[i0]))
            i0 = i1
        return {'slot': np.array([self.iidx[v] for v in post], np.int32), 'ptr': np.array(ptr, np.int64),
                'ch': np.array(ch, np.int32), 'chunks': chunks}

    def _rebuild(self, sch, need_p, ex):
        """The distal vectors of a schedule's nodes (_schedule) and, with need_p, the proximal vectors
        of their internal children (the top node's proximal set by the caller), by the engine: every
        thread of `ex` computes every node for a block of sites (xp_rebuild: each site on its own, so the
        same values as one thread), chunk by chunk, after the chunk's transition matrices (computed on
        the threads too). One call per thread and chunk: the threads never wait for each other's Python."""
        R, gamma, cat, (x, l, xit, per_cat) = self._rebuild_model()
        tsize = self.nrates * self.ns * self.ns
        need = max(len(c[2]) for c in sch['chunks'])
        if self._tensors is None or len(self._tensors) < need:
            self._tensors = shared_zeros((need, self.nrates, self.ns, self.ns))
        ten, L, S = self._tensors, self.lib, self.nsites
        nb = max(1, min(self.cpus, S // BLOCK_SITES))   # more threads only add memory traffic
        step = -(-S // nb)
        blocks = [(s, min(S, s + step)) for s in range(0, S, step)]
        head = (P(self.darr), P(self.dexp), P(self.parr), P(self.pexp), S, R, gamma)
        mid = (P(cat), P(ten), tsize)
        tail = (P(sch['slot']), P(sch['ptr']), P(sch['ch']), P(self.refA), P(self.alpha.lut))
        model = (P(x), P(l), P(xit), P(self.rates), self.nrates, P(self.occupied), per_cat)
        t0 = ten.ctypes.data

        def run(chunk, proximal, n):
            i0, i1, bls, n_ch = chunk
            if n:                                   # the chunk's matrices: n_ch child edges (+ own edges)
                cut = [n * k // min(n, len(blocks)) for k in range(min(n, len(blocks)) + 1)]
                b0 = bls.ctypes.data
                list(ex.map(lambda a, b: L.xp_multi_exp_many(ctypes.c_void_p(t0 + a * tsize * 8),
                                                             ctypes.c_void_p(b0 + a * 8), b - a, *model),
                            cut[:-1], cut[1:]))
            rc = list(ex.map(lambda b: L.xp_rebuild(*head, b[0], b[1], *mid, proximal, i0, i1, *tail), blocks))
            if any(rc):
                raise MemoryError('Not enough memory to rebuild the likelihood vectors.')

        chunks = sch['chunks']
        for k, chunk in enumerate(chunks):          # distal: children first; the last chunk's matrices
            last = need_p and k == len(chunks) - 1  # also serve the first proximal chunk
            run(chunk, 0, len(chunk[2]) if last else chunk[3])
        if need_p:
            for k, chunk in enumerate(reversed(chunks)):    # proximal: parents first
                run(chunk, 1, len(chunk[2]) if k else 0)

    def root_proximal(self):
        """The root's proximal vector (Glv.set_unit): ones."""
        return 1.0

    def _multi_exp(self, tensor, bl):
        """Diagd.multi_exp: the transition matrices of the occupied rates at branch length bl."""
        self.lib.xp_multi_exp(P(tensor), P(self.x), P(self.l), P(self.xit), P(self.rates),
                              self.nrates, P(self.occupied), bl)

    # ---------------------------------------------------------------- prefilter --
    def _query_tables(self):
        """Prefilter inputs: every query character evolved by start_pend under every rate
        (qevtab), and the queries as character codes."""
        ns, alpha = self.ns, self.alpha
        ten = np.zeros((self.nrates, ns, ns))
        self._multi_exp(ten, self.search.start_pend)
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

    def prefilter(self, n_cand=None, query_batch=resources.QUERY_BATCH):
        """pplacer's ranking (h_ranking) for every query, scored 8 edges per call.

        Default: every query's score against every edge is kept (a queries x edges table) and
        ranking orders a query's edges when it is placed. With n_cand (or the placer's
        ranked_edges, --ranked_edges): the queries are scored in batches of query_batch and each
        keeps only its first n_cand edges in ranking order (score descending, ties by ascending
        edge number); memory then barely grows with the number of queries. A query whose search
        goes past them is ranked in full when it gets there (ranking); results are identical."""
        self._query_tables()
        Qn = len(self.query_names)
        self.locs = np.array([v for v in range(self.nnodes) if v != self.root], np.int64)
        E = len(self.locs)
        self.max_edges = E
        K = n_cand if n_cand is not None else self.ranked_edges
        groups = [(i, min(E, i + 8)) for i in range(0, E, 8)]
        if K is None:                                  # every score, one batch of all queries
            self.rank = None
            self.scores = shared_zeros((Qn, E))
            self._bcodes, self._qb = self.codes, (Qn, E)
            start = 0
            if self.workers is None:                   # --cpus auto without a saved benchmark
                start = self._calibrate(groups)
            self._run_groups(groups[start:], self.workers)
            self._bcodes = None
            return
        K = min(max(1, int(K)), E)
        qb = max(1, min(int(query_batch), Qn))
        full = K >= E                                  # every edge, in order
        cs = None if full else np.full((Qn, K), -np.inf)
        ce = np.full((Qn, K), _NO_EDGE, np.int64)
        self.scores = shared_zeros((qb, E))
        try:
            for b0 in range(0, Qn, qb):
                n = min(qb, Qn - b0)
                self._bcodes = np.ascontiguousarray(self.codes[b0:b0 + n])
                self._qb = (n, E)
                start = 0
                if self.workers is None:
                    start = self._calibrate(groups)
                self._run_groups(groups[start:], self.workers)
                sc = self.scores[:n]
                if full:
                    ce[b0:b0 + n] = self.locs[np.lexsort((np.broadcast_to(self.locs, sc.shape), -sc), axis=1)]
                else:
                    merge_candidates(cs, ce, b0, self.locs, sc)
        finally:
            del self.scores
            self._bcodes = None
        self.rank = ce if full else np.take_along_axis(ce, np.lexsort((ce, -cs), axis=1), 1)

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
                for _ in pool.imap_unordered(_batch_prefilter_group, groups, chunksize=chunksize):
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

    def _batch_prefilter_group(self, g):
        """Prefilter scores of the current query batch (_bcodes, _qb) against edges locs[i0:i1]."""
        i0, i1 = g
        nq, _ = self._qb
        t = self._supernodes8(i0, i1)
        st = self.lib.xp_edge_scores8(P(self.scores), self.max_edges, i0, i1 - i0, P(t.sn8), P(t.se8),
                                      P(self.qevtab), P(self.cats), self.nsites, len(self.alpha.codes),
                                      P(self._bcodes), nq)
        if st != 0:
            raise GPlaceException(f'Prefilter kernel failed ({st}) for edges {i0}..{i1}.')

    def ranking(self, qi, K=96):
        """h_ranking for query qi: edges by prefilter score (descending), ties by ascending edge
        number. From the score table (default), or from the query's stored first edges
        (--ranked_edges): a search that goes past them (rare: many failed optimisations) scores qi
        against every edge again and continues in the same order."""
        if self.rank is None:
            row = self.scores[qi]
            E = len(row)
            if K < E:
                kth = np.partition(row, E - K)[E - K]
                cand = np.nonzero(row >= kth)[0]
            else:
                cand = np.arange(E)
            for k in sorted(cand, key=lambda k: (-row[k], self.locs[k])):
                yield int(self.locs[k])
            if len(cand) < E:                              # rare: exhausted candidates
                rest = np.setdiff1d(np.arange(E), cand)
                for k in sorted(rest, key=lambda k: (-row[k], self.locs[k])):
                    yield int(self.locs[k])
            return
        row = self.rank[qi]
        for loc in row:
            yield int(loc)
        if len(row) < len(self.locs):
            full = self._full_ranking(qi)
            if not np.array_equal(full[:len(row)], row):
                raise GPlaceException(f'Inconsistent edge ranking for {self.query_names[qi]}.')
            for loc in full[len(row):]:
                yield int(loc)

    def _full_ranking(self, qi):
        """Every edge in ranking order for query qi, scored in this process."""
        E = len(self.locs)
        saved = {k: getattr(self, k, None) for k in ('scores', '_bcodes', '_qb', 'max_edges')}
        try:
            self.scores = np.zeros((1, E))
            self._bcodes = np.ascontiguousarray(self.codes[qi:qi + 1])
            self._qb, self.max_edges = (1, E), E
            for i in range(0, E, 8):
                self._batch_prefilter_group((i, min(E, i + 8)))
            row = self.scores[0].copy()
        finally:
            for k, v in saved.items():
                setattr(self, k, v)
        return self.locs[np.lexsort((self.locs, -row))]

    # ---------------------------------------------------------------- placement --
    def _marg_prob(self, tt, loc):
        """Three_tax.calc_marg_prob for the placement tt is set to (-p): (marginal log likelihood,
        True if pplacer warned that the integration did not converge after changing bounds)."""
        unif, mu = self.pp_prior
        out = np.empty(2)
        st = self.lib.xp_tt_marg_prob(tt, int(unif), mu if unif or isinstance(mu, float) else mu[loc],
                                      self.search.pp_rel_err, self.search.max_pend, P(out))
        if st == core.XP_GSL:
            raise GPlaceException(f'GSL error {int(out[1])} in the posterior probability calculation '
                                  f'(pplacer: Fatal error: exception Gsl_error.Gsl_exn).')
        if st != core.XP_OK:
            raise GPlaceException(f'pplacer stops here: Fatal error: exception {_PP_ERRORS[st]} '
                                  f'(posterior probability calculation).')
        return float(out[0]), out[1] == 1.

    def _query_evaluator(self, qi):
        """(evaluate, marginal, free) for query qi. evaluate(loc, pend, dist, tol) sets the
        three-taxon tree on edge loc with starting branch lengths pend/dist and optimises it as
        pplacer's safe_ml_optimize_location does (after a GSL error: retry from the default
        lengths); it returns ((like, pend, dist), None), or (None, status) if the location failed.
        marginal(loc, pend, dist) sets it to a placement and returns _marg_prob's result (-p)."""
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
            st = L.xp_tt_optimize(tt, tol, self.search.max_pend, core.MAX_ITER)
            if st in (core.XP_OK, core.XP_EXCEEDED):
                L.xp_tt_results(tt, P(res))
                return (float(res[0]), float(res[1]), float(res[2])), None
            return None, st

        def evaluate(loc, pend, dist, tol):     # safe_ml_optimize_location
            set_edges(loc, pend, dist)
            r, st = optimize(tol)
            if st == core.XP_GSL:
                set_edges(loc, self.search.start_pend, self.bl[loc] / 2.)
                r, st = optimize(tol)
            return r, st

        def marginal(loc, pend, dist):          # tt_edges_from_placement; calc_marg_prob
            set_edges(loc, pend, dist)
            return self._marg_prob(tt, loc)

        arrays = [qm, catm]                     # the engine holds pointers into these until freed

        def free():
            L.xp_tt_free(tt)
            arrays.clear()

        return evaluate, marginal, free

    def place_query(self, qi):
        """Baseball + three-taxon optimisation + keep/refine + LWR (core.ml)."""
        name = self.query_names[qi]
        evaluate, marginal, free = self._query_evaluator(qi)
        try:
            bb = Baseball(self.search)
            for loc in self.ranking(qi):
                r, st = evaluate(loc, self.search.start_pend, self.bl[loc] / 2., self.search.ml_tolerance)
                if not bb.add(loc, r, st):
                    break
            if not bb.results:
                raise GPlaceException(f'empty results for {name}!')
            refined = []
            for loc, (lk, pend, dist) in bb.kept(self.lib):
                r, _ = evaluate(loc, pend, dist, core.FINAL_TOL)
                refined.append((loc, r if r is not None else (lk, pend, dist)))
            margs = None
            if self.search.calc_pp:
                margs = []
                for loc, r in refined:
                    m, warn = marginal(loc, r[1], r[2])
                    margs.append(m)
                    if warn:
                        bb.warns.append((loc, 'pp'))
            pl = weighted_placements(self.lib, refined, margs)
        finally:
            free()
        return name, pl, bb.warns

    def place_all(self, workers=None, progress=True):
        """Place every query; returns [(name, placements, warnings)] in input order.
        placements: [distal_length, edge_num, like_weight_ratio, likelihood, pendant_length].
        progress: log 'Placed x/n queries.' about every 5% (off for the benchmark's timings)."""
        global _PLACER
        workers = workers or self.workers or resources.heuristic_cpus(self.cpus)
        n = len(self.query_names)
        results = [None] * n
        _PLACER = self
        chunks = [range(i, min(n, i + 4)) for i in range(0, n, 4)]
        try:
            with _fork_pool(workers) as pool:
                done_q, step = 0, max(1, n // 20)
                for out in pool.imap_unordered(_place_chunk, chunks):
                    for qi, nm, pl, w in out:
                        results[qi] = (nm, pl, w)
                    if progress and ((done_q + len(out)) // step != done_q // step or done_q + len(out) == n):
                        self.logger.info(f'Placed {done_q + len(out):,}/{n:,} queries.')
                    done_q += len(out)
        finally:
            _PLACER = None
            gc.unfreeze()
        return results

    def release(self):
        """Free the likelihood cache and prefilter scores (the shared-memory arrays)."""
        for attr in ('darr', 'dexp', 'parr', 'pexp', 'scores', 'codes', 'qevtab', 'rank'):
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
