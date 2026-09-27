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

"""Bindings to the compiled engine (libgplace_core.so) and the pplacer model set-up.

Everything here that produces a floating-point value either calls the engine or uses
IEEE-754 double operations in the same order as pplacer's OCaml code.
"""

import ctypes
import math
import os
import re

import numpy as np

from gplace.config.common import CONFIG, ENGINE_NAME
from gplace.exceptions import EngineNotBuilt

NS = 20
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# pplacer v1.1.alpha19 default preferences (prefs.ml), unchanged by GTDB-Tk's invocation
# `pplacer -m WAG -j <cpus> -c <refpkg> -o <out.json> <msa>`.
START_PEND = 0.1
MAX_PEND = 2.0
INITIAL_TOL = 0.01
FINAL_TOL = 1e-5
MAX_ITER = 200
KEEP_AT_MOST = 7
KEEP_FACTOR = 0.01
MAX_PITCHES = 40
MAX_STRIKES = 6
STRIKE_BOX = 3.0

XP_OK, XP_EXCEEDED, XP_GSL, XP_FINDSTART, XP_INVALIDSTART, XP_ASSERT = range(6)

V, I, D = ctypes.c_void_p, ctypes.c_int, ctypes.c_double


def P(a):
    """ctypes pointer to a numpy array's data (the array must stay referenced by the caller)."""
    return ctypes.c_void_p(a.ctypes.data)


def engine_path():
    """Location of the compiled engine: $GPLACE_ENGINE_DIR, else gplace/engine/lib/."""
    return os.path.join(CONFIG.ENGINE_DIR, ENGINE_NAME)


_lib = None


def lib():
    """Load libgplace_core.so once and declare its C signatures."""
    global _lib
    if _lib is not None:
        return _lib
    from gplace.engine.build import cpu_has_avx2
    if not cpu_has_avx2():
        raise EngineNotBuilt('This CPU does not report AVX2 and FMA support, which the gPlace engine requires.')
    path = engine_path()
    if not os.path.exists(path):
        raise EngineNotBuilt(f'The gPlace engine was not found at {path}. It is compiled when '
                             f'gPlace is installed (pip install .); to rebuild it: gplace build_engine')
    L = ctypes.CDLL(path)
    L.xp_init.restype = I
    L.xp_engine_version.restype = ctypes.c_char_p
    L.xp_exp.argtypes = [D]; L.xp_exp.restype = D
    L.xp_log.argtypes = [D]; L.xp_log.restype = D
    L.xp_hypot.argtypes = [D, D]; L.xp_hypot.restype = D
    L.xp_symmv.argtypes = [V, V, V, I]; L.xp_symmv.restype = I
    L.xp_dediag.argtypes = [V, V, V, V, I]
    L.xp_multi_exp.argtypes = [V, V, V, V, V, I, V, D]
    L.xp_multi_exp_ref.argtypes = [V, V, V, V, V, I, V, D]
    L.xp_evolve.argtypes = [V, V, I, V, V]
    L.xp_evolve_ref.argtypes = [V, V, I, V, V]
    L.xp_pairwise.argtypes = [V, V, V, V, V, V, I]
    L.xp_statd_pairwise.argtypes = [V, V, V, V, V, V, V, I]
    L.xp_pull_exponent.argtypes = [V, V, I, I]
    L.xp_total_twoexp.argtypes = [V, I]; L.xp_total_twoexp.restype = D
    L.xp_log_like3.argtypes = [V, V, V, V, I]; L.xp_log_like3.restype = D
    L.xp_edge_scores.argtypes = [V, ctypes.c_long, V, V, V, V, I, I, V, I, V]
    L.xp_edge_scores8.argtypes = [V, ctypes.c_long, ctypes.c_long, I, V, V, V, V, I, I, V, I]
    L.xp_edge_scores8.restype = I
    L.xp_tt_new.argtypes = [I, V, V, V, V, V, V, I, V]; L.xp_tt_new.restype = V
    L.xp_tt_set.argtypes = [V, V, V, V, V, V, D, D, D]
    L.xp_tt_optimize.argtypes = [V, D, D, I]; L.xp_tt_optimize.restype = I
    L.xp_tt_results.argtypes = [V, V]
    L.xp_tt_free.argtypes = [V]
    L.xp_trace_open.argtypes = [ctypes.c_char_p]
    L.xp_init()
    _lib = L
    return L


# ---------------------------------------------------------------------- alphabet --
PROT_CODE = 'ARNDCQEGHILKMFPSTWYV'
CODES = 'ARNDCQEGHILKMFPSTWYVBJZX'      # informative characters, in prefilter-code order


def _prot_map():
    """Prot_models.prot_map: one-hot amino acids, B=D+N, J=I+L, Z=E+Q, X/-/? = all ones."""
    m = {}
    for i, c in enumerate(PROT_CODE):
        v = np.zeros(NS); v[i] = 1.0; m[c] = v
    anyv = np.ones(NS)
    b = np.zeros(NS); b += m['D']; b += m['N']
    j = np.zeros(NS); j += m['I']; j += m['L']
    z = np.zeros(NS); z += m['E']; z += m['Q']
    m.update({'-': anyv, 'B': b, 'J': j, 'Z': z, '?': anyv, 'X': anyv})
    return m


PROT_MAP = _prot_map()
LUT = np.full((256, NS), np.nan)
for _ch, _v in PROT_MAP.items():
    LUT[ord(_ch)] = _v
VALID = np.zeros(256, bool)
for _ch in PROT_MAP:
    VALID[ord(_ch)] = True
GAP_CODE = len(CODES)
CODE_OF = np.full(256, GAP_CODE, np.uint8)
for _i, _c in enumerate(CODES):
    CODE_OF[ord(_c)] = _i


def informative(A):
    """Non-gap, non-'?' characters (pplacer's pre-mask and per-query site mask)."""
    return (A != ord('-')) & (A != ord('?'))


# --------------------------------------------------------------------------- WAG --
def load_wag():
    """WAG exchangeabilities R (symmetric, zero diagonal) and stationary frequencies, parsed
    from the verbatim pplacer v1.1.alpha19 strings (float() == OCaml float_of_string)."""
    section, rows, statd = None, [], None
    for line in open(os.path.join(DATA_DIR, 'wag.txt')):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line in ('[trans]', '[statd]'):
            section = line
            continue
        vals = [float(x) for x in re.split(r'[ \t]+', line)]
        if section == '[trans]':
            rows.append(vals)
        else:
            statd = vals
    n = len(rows) + 1
    R = [[0.0 if i == j else rows[max(i, j) - 1][min(i, j)] for j in range(n)] for i in range(n)]
    return R, statd


def build_diagd(R, pi):
    """Diagd.normed_of_exchangeable_pair, replicating pplacer's operation order.
    Returns (x, l, xit) with Q = x diag(l) xit^T (rate-normalised eigenvalues l)."""
    L = lib()
    n = len(pi)
    B = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = R[i][j]
            else:
                total = 0.0
                for k in range(n):
                    if k != i:
                        total = total + R[k][i] * pi[k]
                inv = 1.0 if pi[i] == 0.0 else 1.0 / pi[i]
                B[i][j] = -(total * inv)
    d_root = [math.sqrt(p) for p in pi]
    d_inv = [1.0 if d == 0.0 else 1.0 / d for d in d_root]
    M = np.ascontiguousarray(np.array([[d_root[i] * (B[i][j] * d_root[j]) for j in range(n)] for i in range(n)]))
    l = np.empty(n); u = np.empty((n, n))
    st = L.xp_symmv(P(M), P(l), P(u), n)
    if st != 0:
        raise RuntimeError(f'gsl_eigen_symmv failed with status {st}')
    x = np.array([[d_inv[i] * u[i][j] for j in range(n)] for i in range(n)])
    xit = np.array([[d_root[i] * u[i][j] for j in range(n)] for i in range(n)])
    Q = np.empty((n, n))
    L.xp_dediag(P(Q), P(x), P(l), P(xit), n)          # Q transposed; only the diagonal is used
    rate = 0.0
    for i in range(n):
        rate = rate - Q[i][i] * pi[i]
    scale = 1.0 / rate
    l = np.array([li * scale for li in l])
    return np.ascontiguousarray(x), np.ascontiguousarray(l), np.ascontiguousarray(xit), M


# ------------------------------------------------------------------------- FASTA --
def read_fasta(path):
    """(name, sequence) pairs; the name is the first whitespace-delimited header token."""
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0] if line[1:].strip() else '', []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf)))
    return out
