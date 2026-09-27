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

1. libm: hypot/exp/log against values returned by the libm inside GTDB-Tk's pplacer binary.
2. Model: the WAG eigendecomposition against the one dumped from GTDB-Tk's pplacer binary.
3. Kernels: every SIMD kernel against its scalar reference, on random and extreme inputs.
"""

import os

import numpy as np

from gplace.engine import core
from gplace.engine.core import NS, P


def _bits(a):
    return np.ascontiguousarray(a).view(np.uint64)


def check_libm():
    L = core.lib()
    n = bad = 0
    for line in open(os.path.join(core.DATA_DIR, 'stock_libm_vectors.tsv')):
        if line.startswith('#') or not line.strip():
            continue
        op, a, b, want = line.split()
        a, b, want = float.fromhex(a), float.fromhex(b), float.fromhex(want)
        got = L.xp_hypot(a, b) if op == 'h' else L.xp_exp(a) if op == 'e' else L.xp_log(a)
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


def check_kernels(seed=7):
    L = core.lib()
    rng = np.random.default_rng(seed)
    R, pi = core.load_wag()
    x, l, xit, _ = core.build_diagd(R, pi)
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
    nc = 24
    for trial in range(4):
        Q = int(rng.integers(1, 200)); ne = int(rng.integers(1, 9))
        qevtab = rng.random((nr, nc, NS)) * 0.2 + 1e-3
        sn8 = rng.random((8, S, NS)) ** rng.uniform(1, 8) * 1e-3
        se8 = rng.integers(-60, 1, (8, S)).astype(np.int64)
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
    return len(bls) * 2 + 4, fails
