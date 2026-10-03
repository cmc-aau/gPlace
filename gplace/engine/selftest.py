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

"""Bit-exactness self-tests of the compiled engine (`gplace check_install`).

1. libm: hypot/exp/log/pow/sin/log1p against values returned by the libm inside GTDB-Tk's pplacer binary.
2. Model: the WAG eigendecomposition against the one dumped from GTDB-Tk's pplacer binary.
3. Kernels: every SIMD kernel against its scalar reference, on random and extreme inputs, in the
   protein (20-state) and DNA (4-state) engines.
4. Gamma: the discrete-gamma rates (GSL 1.16) against values whose placements matched the stock
   pplacer, and the gamma exponent pull and prefilter against straightforward references.
"""

import os

import numpy as np

from gplace.engine import core
from gplace.engine.core import NS, P

M_LN2 = 0.6931471805599453              # C's M_LN2 as a double

# Gamma.discrete_gamma (n_cats, alpha) from this engine; placements with these rates were identical
# to the stock pplacer v1.1.alpha19 (including alpha 0.05, whose first rate is exactly 0).
DISCRETE_GAMMA = {
    (4, 0.5): ['0x1.1183319bd2ba1p-5', '0x1.01f63f0d9f629p-2', '0x1.a3fa3b000863cp-1', '0x1.727c9c97daafep+1'],
    (4, 0.05): ['0x0.0p+0', '0x1.1cfed712609bep-20', '0x1.5b4be6d114afbp-8', '0x1.ff525124a0bd1p+1'],
    (4, 0.8): ['0x1.87689c6631995p-4', '0x1.a0e7be2fe4901p-2', '0x1.e9f5fce96cd25p-1', '0x1.452a441c76acap+1'],
    (8, 0.7): ['0x1.b63f26ad1d218p-6', '0x1.f2304d3d5a3b0p-4', '0x1.0f7f3e502763ap-2', '0x1.d90346b10db3ap-2',
               '0x1.767676cf1b6bfp-1', '0x1.1e4d1eaf265eap+0', '0x1.bf87cf840461ap+0', '0x1.c36999db3815ep+1'],
    (4, 0.289): ['0x1.235776bb8ead2p-8', '0x1.8917b1edeaca8p-4', '0x1.275a02863c1f5p-1', '0x1.a94f1613a3da8p+1'],
}


def _bits(a):
    return np.ascontiguousarray(a).view(np.uint64)


def check_libm(ns=NS):
    L = core.lib(ns)
    n = bad = 0
    for line in open(os.path.join(core.DATA_DIR, 'stock_libm_vectors.tsv')):
        if line.startswith('#') or not line.strip():
            continue
        op, a, b, want = line.split()
        a, b, want = float.fromhex(a), float.fromhex(b), float.fromhex(want)
        got = {'h': lambda: L.xp_hypot(a, b), 'e': lambda: L.xp_exp(a), 'l': lambda: L.xp_log(a),
               'p': lambda: L.xp_pow(a, b), 's': lambda: L.xp_sin(a), 'q': lambda: L.xp_log1p(a)}[op]()
        n += 1
        bad += _bits(np.array([got]))[0] != _bits(np.array([want]))[0]
    return n, bad


def check_wag_eigen():
    """Same input matrix and bit-identical eigenvalues/eigenvectors as the stock binary."""
    L = core.lib()
    sec, vals = None, {'input': [], 'eval': [], 'evec': []}
    for line in open(os.path.join(core.DATA_DIR, 'stock_wag_eigen.txt')):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('['):
            sec = line[1:-1]; continue
        vals[sec].append(float.fromhex(line))
    R, pi = core.load_wag()
    _, _, _, M = core.build_diagd(R, pi)
    same_input = np.array_equal(_bits(M.ravel()), _bits(np.array(vals['input'])))
    Min = np.ascontiguousarray(np.array(vals['input']).reshape(NS, NS))
    ev = np.empty(NS); U = np.empty((NS, NS))
    L.xp_symmv(P(Min), P(ev), P(U), NS)
    same_eval = np.array_equal(_bits(ev), _bits(np.array(vals['eval'])))
    same_evec = np.array_equal(_bits(U.ravel()), _bits(np.array(vals['evec'])))
    return same_input, same_eval, same_evec


def _test_model(ns):
    """WAG for the protein engine; a GTR model for the DNA engine."""
    if ns == 20:
        return core.load_wag()
    from gplace.engine.models import gtr_exchangeabilities
    R = gtr_exchangeabilities({'ac': 1.27, 'ag': 2.94, 'at': 0.71, 'cg': 1.08, 'ct': 4.6, 'gt': 1.0})
    return R, [0.31, 0.18, 0.24, 0.27]


def check_kernels(seed=7, ns=NS):
    """Every kernel of the ns-state engine against its scalar reference."""
    L = core.lib(ns)
    NS = ns                                                          # noqa: N806 (shadows the module NS)
    rng = np.random.default_rng(seed)
    R, pi = _test_model(ns)
    x, l, xit, _ = core.build_diagd(R, pi, L)
    nr = 20
    rates = np.sort(rng.uniform(0.02, 5.0, nr))
    S = 1500
    cats = rng.integers(0, nr, S).astype(np.int32)
    occ = np.zeros(nr, np.uint8); occ[np.unique(cats)] = 1
    fails = []
    bls = np.concatenate([[0.0, 1e-8, 5e-324, 0.1, 2.0], rng.uniform(0, 3, 60), 10 ** rng.uniform(-9, 0.5, 60)])
    for bl in bls:
        tT = np.zeros((nr, NS, NS)); tR = np.zeros((nr, NS, NS))
        L.xp_multi_exp(P(tT), P(x), P(l), P(xit), P(rates), nr, P(occ), float(bl))
        L.xp_multi_exp_ref(P(tR), P(x), P(l), P(xit), P(rates), nr, P(occ), float(bl))
        if not np.array_equal(_bits(tT.transpose(0, 2, 1)), _bits(tR)):
            fails.append(f'multi_exp bl={bl!r}')
        src = rng.random((S, NS)) ** rng.uniform(1, 30)
        src[rng.random((S, NS)) < 0.05] = 0.0
        a = np.empty((S, NS)); b = np.empty((S, NS))
        L.xp_evolve(P(a), P(src), S, P(cats), P(tT))
        L.xp_evolve_ref(P(b), P(src), S, P(cats), P(tR))
        if not np.array_equal(_bits(a), _bits(b)):
            fails.append(f'evolve bl={bl!r}')
    nc = len(core.ALPHABETS[ns].codes)
    # exponent ranges: class-level trees; the full GTDB bacterial tree (~2^18 per site, totals
    # beyond int32); beyond the SIMD path's per-site bound (exact scalar fallback)
    se_ranges = [(-60, 1)] * 4 + [(-(1 << 21), -(1 << 17))] * 2 + [(-(1 << 24), 1)]
    for trial, (se_lo, se_hi) in enumerate(se_ranges):
        Q = int(rng.integers(1, 200)); ne = int(rng.integers(1, 9))
        qevtab = rng.random((nr, nc, NS)) * 0.2 + 1e-3
        sn8 = rng.random((8, S, NS)) ** rng.uniform(1, 8) * 1e-3
        se8 = rng.integers(se_lo, se_hi, (8, S)).astype(np.int64)
        codes255 = rng.integers(0, nc, (Q, S)).astype(np.uint8)
        gaps = rng.random((Q, S)) < 0.3
        codes255[gaps] = 255
        codes24 = codes255.copy(); codes24[gaps] = nc
        out8 = np.full((Q, 11), np.nan)
        L.xp_edge_scores8(P(out8), 11, 2, ne, P(sn8), P(se8), P(qevtab), P(cats), S, nc, P(codes24), Q)
        W = np.empty(S * nc)
        for lane in range(ne):
            ref = np.empty(Q)
            sn = np.ascontiguousarray(sn8[lane]); se = np.ascontiguousarray(se8[lane])
            L.xp_edge_scores(P(ref), 1, P(sn), P(se), P(qevtab), P(cats), S, nc, P(codes255), Q, P(W))
            if not np.array_equal(_bits(out8[:, 2 + lane]), _bits(ref)):
                fails.append(f'edge_scores8 trial {trial} lane {lane}')
    return len(bls) * 2 + len(se_ranges), fails


def check_discrete_gamma():
    """Returns the (n_cats, alpha) cases whose rates differ from DISCRETE_GAMMA."""
    from gplace.engine.models import discrete_gamma
    return [k for k, want in DISCRETE_GAMMA.items()
            if [v.hex() for v in discrete_gamma(*k)] != [float.fromhex(w).hex() for w in want]]


def _pull_reference(a, e, R, ns, min_twoexp, fix):
    """Gmix perhaps_pull_exponent in numpy: per site the largest frexp exponent of the rates'
    maxima (pplacer: frexp 0.0 = exponent 0; fix: all-zero rates skipped); exact scaling."""
    a = a.copy(); e = e.copy()
    for s in range(a.shape[0]):
        site = a[s].reshape(R, ns)
        mx = site.max(axis=1)
        exps = [int(np.frexp(m)[1]) for m in mx if not (fix and m == 0.0)]
        if not exps:
            continue
        t = max(exps)
        if t < min_twoexp:
            a[s] = (site * np.ldexp(0.5, -t + 1)).ravel()
            e[s] += t
    return a, e


def check_gamma_kernels(seed=11, ns=NS):
    """The Gmix exponent pull (with and without the underflow fix) and the Gmix prefilter of the
    ns-state engine against straightforward references (sequential sums, the engine's log)."""
    L = core.lib(ns)
    rng = np.random.default_rng(seed)
    fails, n = [], 0
    for R in (4, 8):
        S = 300
        a = rng.random((S, R * ns)) ** rng.uniform(1, 40) * 2.0 ** rng.integers(-700, 1, (S, 1))
        dead = rng.random((S, R)) < 0.3                       # rates underflowed to exactly 0.0
        a3 = a.reshape(S, R, ns); a3[dead] = 0.0
        e = rng.integers(-5000, 1, S).astype(np.int64)
        for fix in (0, 1):
            L.xp_set_gamma_fix(fix)
            got_a = np.ascontiguousarray(a.copy()); got_e = e.copy()
            L.xp_g_pull_exponent(P(got_a), P(got_e), S, R, -50)
            ref_a, ref_e = _pull_reference(a, e, R, ns, -50, fix)
            n += 1
            if not (np.array_equal(_bits(got_a), _bits(ref_a)) and np.array_equal(got_e, ref_e)):
                fails.append(f'g_pull_exponent R={R} fix={fix}')
        L.xp_set_gamma_fix(0)
        nc = len(core.ALPHABETS[ns].codes)
        Sp, Q, ne = 25, 6, int(rng.integers(1, 9))
        qevtab = rng.random((R, nc, ns)) * 0.2 + 1e-3
        sn8 = rng.random((8, Sp, R * ns)) * 1e-3
        se8 = rng.integers(-3000, 1, (8, Sp)).astype(np.int64)
        codes = rng.integers(0, nc + 1, (Q, Sp)).astype(np.uint8)       # nc = gap
        nused = (codes != nc).sum(axis=1).astype(np.int32)
        out = np.full((Q, 10), np.nan)
        L.xp_g_edge_scores8(P(out), 10, 1, ne, P(sn8), P(se8), P(qevtab), R, Sp, nc, P(codes), P(nused), Q)
        logR = L.xp_log(float(np.float32(R)))
        for lane in range(ne):
            for q in range(Q):
                ll, et = 0.0, 0
                for s in range(Sp):
                    c = codes[q, s]
                    if c == nc:
                        continue
                    util = 0.0
                    for r in range(R):
                        for st in range(ns):
                            util += float(qevtab[r, c, st]) * float(sn8[lane, s, r * ns + st])
                    ll += L.xp_log(util)
                    et += int(se8[lane, s])
                ll -= float(np.float32(nused[q])) * logR
                want = ll + (M_LN2 * (0. + float(et)))
                n += 1
                if _bits(np.array([out[q, 1 + lane]]))[0] != _bits(np.array([want]))[0]:
                    fails.append(f'g_edge_scores8 R={R} lane={lane} q={q}')
    return n, fails
