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

    def __init__(self, lib, R, mstatd=None):
        self._l, self.R, self.mstatd = lib, R, mstatd

    def __getattr__(self, name):
        return getattr(self._l, name)

    def xp_pairwise(self, dst, de, x, xe, y, ye, S):
        self._l.xp_g_pairwise(dst, de, x, xe, y, ye, S, self.R)

    def xp_statd_pairwise(self, dst, de, statd, a, ae, b, be, S):
        if self.mstatd is not None:              # mixture: every category's own frequencies
            self._l.xp_m_statd_pairwise(dst, de, P(self.mstatd), a, ae, b, be, S, self.R)
        else:
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
        self.mix = self._mixture_setup() if getattr(self.model, 'components', None) else None
        self.lib = _GammaLib(self.lib, self.R, self.mix['mstatd'] if self.mix else None)
        self.set_gamma_fix(False)
        # --full_model: an invariable-sites category (rate 0) evolves with the exact identity (per process)
        self.lib._l.xp_set_zero_rate_identity(1 if getattr(self.model, 'full_model', False) else 0)

    def root_proximal(self):
        """The root's proximal vector: ones (pplacer: likelihoods are the mean over the rate categories),
        or, with category weights w (--full_model: +I, +R), R * w_r in category r, so that this mean is
        the weighted sum. Every proximal vector, and so every likelihood, carries this factor once."""
        w = getattr(self.model, 'weights', None)
        if w is None:
            return 1.0
        v = np.empty((self.nsites, self.R, self.ns))
        v[:] = (np.asarray(w, dtype=float) * self.R)[None, :, None]
        return v.reshape(self.nsites, self.width)

    def set_gamma_fix(self, on):
        """pplacer's rescaling (False, default) or the underflow fix (True); per process, inherited
        by the forked workers."""
        self.gamma_fix = bool(on)
        self.lib._l.xp_set_gamma_fix(1 if on else 0)

    def _mixture_setup(self):
        """A mixture (--full_model): every component's eigen-system, and per category the arrays of its
        component (eigen-system, frequencies) for the engine; categories are grouped by component."""
        eig = [core.build_diagd(Rk, list(pik), self.lib)[:3] for Rk, pik in self.model.components]
        comp = self.model.cat_comp
        mx = np.ascontiguousarray([eig[k][0] for k in comp], dtype=float)
        ml = np.ascontiguousarray([eig[k][1] for k in comp], dtype=float)
        mxit = np.ascontiguousarray([eig[k][2] for k in comp], dtype=float)
        mstatd = np.ascontiguousarray([self.model.components[k][1] for k in comp], dtype=float)
        blocks, c0 = [], 0
        for c in range(1, len(comp) + 1):
            if c == len(comp) or comp[c] != comp[c0]:
                x, l, xit = (np.ascontiguousarray(a, dtype=float) for a in eig[comp[c0]])
                blocks.append((c0, c, x, l, xit)); c0 = c
        return {'mx': mx, 'ml': ml, 'mxit': mxit, 'mstatd': mstatd, 'blocks': blocks}

    def _multi_exp(self, tensor, bl):
        if self.mix is None:
            return super()._multi_exp(tensor, bl)
        for c0, c1, x, l, xit in self.mix['blocks']:     # each component's categories, its eigen-system
            self.lib.xp_multi_exp(P(tensor[c0:c1]), P(x), P(l), P(xit), P(self.rates[c0:c1]), c1 - c0,
                                  P(self.occupied[c0:c1]), bl)

    def evolve_t(self, dst, src, bl, tensor):
        """Gmix evolve_into: tensor_mul (gemmish, rate by rate) performs the same products and
        additions as the CAT evolve on nsites x R virtual sites where virtual site (s, r) uses
        rate r."""
        self._multi_exp(tensor, bl)
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

    def _batch_prefilter_group(self, g):
        """The prefilter of a query batch with the Gmix kernel (untiled and tiled); the
        informative-site counts of the batch are computed once per batch in each process."""
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
        if self.mix:
            mx = self.mix
            tt = L.xp_mtt_new(m, self.R, P(mx['mstatd']), P(mx['mx']), P(mx['ml']), P(mx['mxit']), P(self.rates), P(qm))
        else:
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
