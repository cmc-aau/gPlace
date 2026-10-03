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

"""pplacer's gamma mixture model (Gmix_model): every site is evaluated under all R discrete-gamma
rates and the site likelihood is their average.

Likelihood vectors are nsites x R x NS (site-major; pplacer stores them rate-major, which changes
no value: only the per-site accumulation order matters, and it is pplacer's), with ONE exponent
per site shared by the R rates. The placement algorithm (cache, prefilter, baseball, refinement;
untiled or tiled) is Placer's / TiledPlacer's, unchanged: GammaMixin only swaps the model-specific
operations.

pplacer's exponent rescaling for gamma has a flaw: frexp(0.0) returns exponent 0, so a rate whose
values have all underflowed to 0.0 blocks the rescaling of its site, and on large trees the other
rates then underflow too (likelihood -inf). By default gPlace reproduces it; set_gamma_fix(True)
skips all-zero rates when the exponent is chosen (not pplacer's behaviour).
"""

import numpy as np

from gplace.engine import core
from gplace.engine.core import P
from gplace.engine.placer import Placer
from gplace.engine.tiled import TiledPlacer
from gplace.exceptions import GPlaceException


class _GammaLib(object):
    """The engine seen through Gmix: the generic pairwise product, statd product and exponent pull
    take the site count and act on all R rates of each site (one exponent per site)."""

    def __init__(self, lib, R):
        self._l, self.R = lib, R

    def __getattr__(self, name):
        return getattr(self._l, name)

    def xp_pairwise(self, dst, de, x, xe, y, ye, S):
        self._l.xp_g_pairwise(dst, de, x, xe, y, ye, S, self.R)

    def xp_statd_pairwise(self, dst, de, statd, a, ae, b, be, S):
        self._l.xp_g_statd_pairwise(dst, de, statd, a, ae, b, be, S, self.R)

    def xp_pull_exponent(self, a, e, S, min_twoexp):
        self._l.xp_g_pull_exponent(a, e, S, self.R, min_twoexp)


class GammaMixin(object):
    """Gmix_model operations for Placer and TiledPlacer (GammaPlacer, GammaTiledPlacer)."""

    def __init__(self, refpkg_dir, queries, cpus=1, **kw):
        super().__init__(refpkg_dir, queries, cpus, **kw)
        self._gamma_setup()

    def _gamma_setup(self):
        if self.model.ras != 'gamma':
            raise GPlaceException('GammaPlacer needs a gamma model.')
        self.R = self.nrates
        self.occupied = np.ones(self.R, np.uint8)           # Diagd.multi_exp without a mask
        self.width = self.R * self.ns                        # doubles per site
        self.cat_v = np.ascontiguousarray(np.tile(np.arange(self.R, dtype=np.int32), self.nsites))
        self.lib = _GammaLib(self.lib, self.R)
        self.set_gamma_fix(False)

    def set_gamma_fix(self, on):
        """pplacer's rescaling (False, default) or the underflow fix (True); per process, inherited
        by the forked workers."""
        self.gamma_fix = bool(on)
        self.lib._l.xp_set_gamma_fix(1 if on else 0)

    def evolve_t(self, dst, src, bl, tensor):
        """Gmix evolve_into: tensor_mul (gemmish, rate by rate) performs the same products and
        additions as the CAT evolve on nsites x R virtual sites where virtual site (s, r) uses
        rate r."""
        self.lib.xp_multi_exp(P(tensor), P(self.x), P(self.l), P(self.xit), P(self.rates),
                              self.R, P(self.occupied), bl)
        self.lib.xp_evolve(P(dst), P(src), self.nsites * self.R, P(self.cat_v), P(tensor))

    def leaf_lv(self, v):
        """make_constant_rate_glv_from_lv_arr: the leaf's likelihood vectors under every rate."""
        return np.ascontiguousarray(np.repeat(self.alpha.lut[self.refA[self.leaf_row[v]]], self.R, axis=0)
                                    .reshape(self.nsites, self.width))

    def _query_tables(self):
        super()._query_tables()                              # qevtab for every rate (all occupied)
        self.nused = self._count_used(self.codes)

    def _count_used(self, codes):
        """Informative sites of each query (n_used of ten_masked_logdot_c)."""
        return np.ascontiguousarray((codes != self.alpha.gap_code).sum(axis=1).astype(np.int32))

    def _prefilter_group(self, g):
        i0, i1 = g
        S, Qn, E, nc = self.nsites, len(self.query_names), self.scores.shape[1], len(self.alpha.codes)
        t = self._supernodes8(i0, i1)
        st = self.lib.xp_g_edge_scores8(P(self.scores), E, i0, i1 - i0, P(t.sn8), P(t.se8),
                                        P(self.qevtab), self.R, S, nc, P(self.codes), P(self.nused), Qn)
        if st != 0:
            raise GPlaceException(f'Prefilter kernel failed ({st}) for edges {i0}..{i1}.')

    def _tile_prefilter_group(self, g):
        """TiledPlacer's per-tile prefilter with the Gmix kernel; the informative-site counts of
        the current query batch are computed once per batch in each process."""
        i0, i1 = g
        nq, _ = self._qb
        if getattr(self, '_bnused_src', None) is not self._bcodes:
            self._bnused = self._count_used(self._bcodes)
            self._bnused_src = self._bcodes
        t = self._supernodes8(i0, i1)
        st = self.lib.xp_g_edge_scores8(P(self.scores), self.max_edges, i0, i1 - i0, P(t.sn8), P(t.se8),
                                        P(self.qevtab), self.R, self.nsites, len(self.alpha.codes),
                                        P(self._bcodes), P(self._bnused), nq)
        if st != 0:
            raise GPlaceException(f'Prefilter kernel failed ({st}) for edges {i0}..{i1}.')

    def _query_evaluator(self, qi):
        """Placer._query_evaluator with the Gmix three-taxon tree."""
        L = self.lib
        t = self.tl()
        qrow = self.QM[qi]
        reind = np.nonzero(core.informative(qrow))[0].astype(np.int32)
        m = len(reind)
        qm = np.ascontiguousarray(self.alpha.lut[qrow][reind])
        tt = L.xp_gtt_new(m, self.R, P(self.statd), P(self.x), P(self.l), P(self.xit), P(self.rates), P(qm))
        if not tt:
            raise GPlaceException(f'Could not set up the likelihood for {self.query_names[qi]}.')
        res = np.empty(4)
        D, De, Pp, Pe = t.b[0], t.e[0], t.b[1], t.e[1]
        cur = [None]

        def set_edges(loc, pend, dist):
            if cur[0] != loc:
                self.get_d(loc, D, De); self.get_p(loc, Pp, Pe); cur[0] = loc
            L.xp_gtt_set(tt, P(D), P(De), P(Pp), P(Pe), P(reind), self.bl[loc], pend, dist)

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

        arrays = [qm, reind]                    # the engine holds pointers into these until freed

        def free():
            L.xp_tt_free(tt)
            arrays.clear()

        return evaluate, free


class GammaPlacer(GammaMixin, Placer):
    """Placer for gamma-model reference packages."""


class GammaTiledPlacer(GammaMixin, TiledPlacer):
    """TiledPlacer for gamma-model reference packages."""
