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

from gplace.config.common import CONFIG, ENGINE_NAMES, ENGINE_VARIANTS, engine_file
from gplace.exceptions import EngineNotBuilt, GPlaceException

NS = 20
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# pplacer v1.1.alpha19 default preferences (prefs.ml), unchanged by GTDB-Tk's invocation
# `pplacer -m WAG -j <cpus> -c <refpkg> -o <out.json> <msa>`; Search holds a run's values.
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


class Search(object):
    """pplacer's search settings (--ml-tolerance, --start-pend, --max-pend, --max-strikes,
    --strike-box, --max-pitches, --keep-at-most, --keep-factor, --no-pre-mask), checked as
    Prefs.check does. The final tolerance (1e-5) and the iteration limit (200) are fixed in pplacer."""

    def __init__(self, ml_tolerance=INITIAL_TOL, start_pend=START_PEND, max_pend=MAX_PEND,
                 max_strikes=MAX_STRIKES, strike_box=STRIKE_BOX, max_pitches=MAX_PITCHES,
                 keep_at_most=KEEP_AT_MOST, keep_factor=KEEP_FACTOR, no_pre_mask=False):
        self.ml_tolerance, self.start_pend, self.max_pend = float(ml_tolerance), float(start_pend), float(max_pend)
        self.max_strikes, self.strike_box, self.max_pitches = int(max_strikes), float(strike_box), int(max_pitches)
        self.keep_at_most, self.keep_factor = int(keep_at_most), float(keep_factor)
        self.no_pre_mask = bool(no_pre_mask)
        for k in ('ml_tolerance', 'start_pend', 'max_pend', 'strike_box', 'keep_factor'):
            if math.isnan(getattr(self, k)):
                raise GPlaceException(f'{k} must be a number, not nan.')
        if self.start_pend <= 0.:
            raise GPlaceException('Starting pendant branch length must be strictly positive.')
        if self.start_pend >= self.max_pend:
            raise GPlaceException('Starting pendant branch length must be strictly less than maximum '
                                  'pendant branch length.')
        # core.ml: max_strikes = 0 disables ball playing (every location is evaluated)
        if self.max_strikes == 0:
            self.t_max_pitches = self.t_max_strikes = math.inf
        else:
            self.t_max_pitches, self.t_max_strikes = self.max_pitches, self.max_strikes

    @classmethod
    def from_options(cls, options):
        """The settings of a parsed command line (attributes missing from it keep the defaults)."""
        names = ('ml_tolerance', 'start_pend', 'max_pend', 'max_strikes', 'strike_box', 'max_pitches',
                 'keep_at_most', 'keep_factor', 'no_pre_mask')
        return cls(**{k: getattr(options, k) for k in names if getattr(options, k, None) is not None})

    def changed(self):
        """[(name, value)] of the settings that differ from pplacer's defaults."""
        d = Search()
        return [(k, v) for k, v in vars(self).items() if not k.startswith('t_') and v != getattr(d, k)]

XP_OK, XP_EXCEEDED, XP_GSL, XP_FINDSTART, XP_INVALIDSTART, XP_ASSERT = range(6)

V, I, D = ctypes.c_void_p, ctypes.c_int, ctypes.c_double


def P(a):
    """ctypes pointer to a numpy array's data (the array must stay referenced by the caller)."""
    return ctypes.c_void_p(a.ctypes.data)


def engine_variant():
    """The engine variant this machine uses: $GPLACE_ENGINE_VARIANT (avx2, neon or portable) if
    set, else avx2 on x86-64 CPUs with AVX2 and FMA, neon on ARM64 (when that engine is installed),
    else portable (same results, scalar kernels)."""
    import platform
    from gplace.engine.build import cpu_has_avx2
    v = os.getenv('GPLACE_ENGINE_VARIANT', '').strip().lower()
    if v:
        if v not in ENGINE_VARIANTS:
            raise EngineNotBuilt(f'GPLACE_ENGINE_VARIANT must be one of {", ".join(ENGINE_VARIANTS)}, not {v!r}.')
        return v
    simd = os.path.exists(os.path.join(CONFIG.ENGINE_DIR, engine_file(NS, 'avx2')))
    if cpu_has_avx2() and simd:
        return 'avx2'
    if platform.machine().lower() in ('aarch64', 'arm64') and simd:
        return 'neon'
    return 'portable'


def engine_path(ns=NS, variant=None):
    """Location of the compiled engine for ns states (20: protein, 4: DNA) and a variant (default:
    engine_variant()): $GPLACE_ENGINE_DIR, else gplace/engine/lib/."""
    return os.path.join(CONFIG.ENGINE_DIR, engine_file(ns, variant or engine_variant()))


_libs = {}


def lib(ns=NS):
    """Load the engine for ns states once and declare its C signatures. Both engines can be loaded
    in the same process (separate libraries, local symbols)."""
    if ns in _libs:
        return _libs[ns]
    from gplace.engine.build import cpu_has_avx2
    variant = engine_variant()
    if variant == 'avx2' and not cpu_has_avx2():
        raise EngineNotBuilt('The avx2 engine needs an x86-64 CPU with AVX2 and FMA; this one has '
                             'neither (unset GPLACE_ENGINE_VARIANT to use the portable engine).')
    import platform
    if variant == 'neon' and platform.machine().lower() not in ('aarch64', 'arm64'):
        raise EngineNotBuilt('The neon engine is for ARM64 CPUs (unset GPLACE_ENGINE_VARIANT to use '
                             'this machine\'s engine).')
    if ns not in ENGINE_NAMES:
        raise EngineNotBuilt(f'No engine for {ns} states.')
    path = engine_path(ns)
    if not os.path.exists(path):
        raise EngineNotBuilt(f'The gPlace engine was not found at {path}. It is compiled when '
                             f'gPlace is installed (pip install .); to rebuild it: gplace build_engine')
    L = ctypes.CDLL(path)
    L.xp_n_states.restype = I
    if L.xp_n_states() != ns:
        raise EngineNotBuilt(f'{path} is an engine for {L.xp_n_states()} states, not {ns}; rebuild it with '
                             f'gplace build_engine.')
    L.xp_simd.restype = I
    if L.xp_simd() != {'portable': 0, 'avx2': 1, 'neon': 2}[variant]:
        raise EngineNotBuilt(f'{path} is not a {variant} engine; rebuild it with gplace build_engine.')
    L.xp_init.restype = I
    L.xp_engine_version.restype = ctypes.c_char_p
    L.xp_exp.argtypes = [D]; L.xp_exp.restype = D
    L.xp_log.argtypes = [D]; L.xp_log.restype = D
    L.xp_hypot.argtypes = [D, D]; L.xp_hypot.restype = D
    L.xp_pow.argtypes = [D, D]; L.xp_pow.restype = D
    L.xp_sin.argtypes = [D]; L.xp_sin.restype = D
    L.xp_log1p.argtypes = [D]; L.xp_log1p.restype = D
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
    # gamma mixture model (Gmix)
    L.xp_discrete_gamma.argtypes = [I, D, V]; L.xp_discrete_gamma.restype = I
    L.xp_gamma_gsl_error.argtypes = [ctypes.c_char_p, ctypes.c_size_t]; L.xp_gamma_gsl_error.restype = I
    L.xp_g_pairwise.argtypes = [V, V, V, V, V, V, I, I]
    L.xp_g_statd_pairwise.argtypes = [V, V, V, V, V, V, V, I, I]
    L.xp_m_statd_pairwise.argtypes = [V, V, V, V, V, V, V, I, I]
    L.xp_g_pull_exponent.argtypes = [V, V, I, I, I]
    L.xp_g_edge_scores8.argtypes = [V, ctypes.c_long, ctypes.c_long, I, V, V, V, I, I, I, V, V, I]
    L.xp_g_edge_scores8.restype = I
    L.xp_gtt_new.argtypes = [I, I, V, V, V, V, V, V]; L.xp_gtt_new.restype = V
    L.xp_mtt_new.argtypes = [I, I, V, V, V, V, V, V]; L.xp_mtt_new.restype = V
    L.xp_gtt_set.argtypes = [V, V, V, V, V, V, D, D, D]
    L.xp_set_gamma_fix.argtypes = [I]
    L.xp_set_zero_rate_identity.argtypes = [I]
    L.xp_init()
    _libs[ns] = L
    return L


# ---------------------------------------------------------------------- alphabet --
class Alphabet(object):
    """A sequence alphabet as pplacer uses it: character -> likelihood vector (char_map), the
    informative characters in prefilter-code order (codes), and lookup tables over bytes.
    '-' and '?' are the only non-informative characters (Alignment.informative)."""

    def __init__(self, name, char_map, codes):
        self.name, self.char_map, self.codes = name, char_map, codes
        self.ns = len(next(iter(char_map.values())))
        self.lut = np.full((256, self.ns), np.nan)
        self.valid = np.zeros(256, bool)
        for ch, v in char_map.items():
            self.lut[ord(ch)] = v
            self.valid[ord(ch)] = True
        self.gap_code = len(codes)
        self.code_of = np.full(256, self.gap_code, np.uint8)
        for i, c in enumerate(codes):
            self.code_of[ord(c)] = i


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


PROTEIN = Alphabet('protein', PROT_MAP, CODES)


def _nuc_map():
    """Nuc_models.nuc_map (order A C G T): IUPAC ambiguity codes as sums of one-hot vectors;
    N, X, '-' and '?' are all ones."""
    bases = 'ACGT'
    ambig = {'R': 'AG', 'Y': 'CT', 'M': 'AC', 'K': 'GT', 'S': 'CG', 'W': 'AT', 'H': 'ACT', 'B': 'CGT',
             'V': 'ACG', 'D': 'AGT', 'N': 'ACGT', '-': 'ACGT', '?': 'ACGT', 'X': 'ACGT'}
    m = {c: np.eye(4)[i] for i, c in enumerate(bases)}
    m.update({c: np.array([1.0 if b in bs else 0.0 for b in bases]) for c, bs in ambig.items()})
    return m


NUC_MAP = _nuc_map()
DNA = Alphabet('DNA', NUC_MAP, 'ACGTRYMKSWHBVDNX')
ALPHABETS = {20: PROTEIN, 4: DNA}


def informative(A):
    """Non-gap, non-'?' characters (pplacer's pre-mask and per-query site mask)."""
    return (A != ord('-')) & (A != ord('?'))


# --------------------------------------------------------------------------- WAG --
def load_wag():
    """WAG exchangeabilities and frequencies (see load_protein_model)."""
    return load_protein_model('WAG')


def load_protein_model(name):
    """Exchangeabilities R (symmetric, zero diagonal) and stationary frequencies of a protein model
    (WAG, LG or JTT), parsed from the verbatim pplacer v1.1.alpha19 strings in data/<name>.txt
    (float() == OCaml float_of_string)."""
    section, rows, statd = None, [], None
    with open(os.path.join(DATA_DIR, name.lower() + '.txt')) as fh:
        for line in fh:
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


def build_diagd(R, pi, L=None):
    """Diagd.normed_of_exchangeable_pair, replicating pplacer's operation order.
    Returns (x, l, xit) with Q = x diag(l) xit^T (rate-normalised eigenvalues l)."""
    L = L or lib()
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


# -------------------------------------------------------------------- alignments --
def _open_text(path):
    """A text file, gunzipped if its name ends in .gz (pplacer's MaybeZipped)."""
    if path.endswith('.gz'):
        import gzip
        return gzip.open(path, 'rt')
    return open(path)


def read_fasta(path):
    """(name, sequence) pairs; the name is the first whitespace-delimited header token. As
    pplacer's Fasta: ';' lines are comments, a .gz file is read gunzipped, and a sequence before
    the first name is an error."""
    out, name, buf = [], None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0] if line[1:].strip() else '', []
            elif line.startswith(';'):
                continue
            elif name is None:
                if line.strip():
                    raise GPlaceException(f'Syntax error in {path}: sequence before name.')
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf)))
    return out


_STOCKHOLM_LINE = re.compile(r'(?:|(# STOCKHOLM 1\.0)|(//)|(#=.+)|(([^ \t\n\r#]+)[ \t]+([^ \t\n\r#]+)))'
                             r'[ \t]*(?:[\r\n]+|$)')


def read_stockholm(path):
    """(name, sequence) pairs of a Stockholm alignment, as pplacer's Stockholm module reads it:
    one alignment ("# STOCKHOLM 1.0" ... "//"), #= markup ignored, the blocks of a name
    concatenated, '.' read as '-'; the result is in DESCENDING name order (pplacer folds a string
    map into a list), which is the order pplacer uses everywhere after reading."""
    def err(msg):
        raise GPlaceException(f'Syntax error in {path}: {msg}')

    seqs, state = {}, 'header'
    with _open_text(path) as fh:
        text = fh.read()
    for lineno, line in enumerate(text.split('\n'), 1):
        pos = 0
        while pos < len(line):
            m = _STOCKHOLM_LINE.match(line, pos)
            if not m or m.end() == pos:
                err(f'line {lineno} character {pos}')
            pos = m.end()
            kind = ('header' if m.group(1) else 'footer' if m.group(2) else 'markup' if m.group(3)
                    else 'seq' if m.group(4) else None)
            if kind is None:                         # blank
                continue
            if state == 'header':
                if kind != 'header':
                    err('something before header')
                state = 'body'
            elif state == 'body':
                if kind == 'header':
                    err('header unexpected')
                elif kind == 'seq':
                    seqs.setdefault(m.group(5), []).append(m.group(6).replace('.', '-'))
                elif kind == 'footer':
                    state = 'done'
            else:
                err('something after footer')
    if state != 'done':
        err("didn't reach footer by EOF")
    out = []
    for name in sorted(seqs):                        # StringMap.fold, consing: descending order
        seq = ''.join(seqs[name])
        if out and len(seq) != len(out[-1][1]):
            err(f"{name} (length {len(seq)}) doesn't match the previous sequence lengths ({len(out[-1][1])})")
        out.append((name, seq))
    return out[::-1]


def read_alignment(path):
    """A query alignment as pplacer reads it (Alignment.list_of_any_file): Stockholm if the name
    (without .gz) ends in .sto or .sth, otherwise FASTA."""
    base = path[:-3] if path.endswith('.gz') else path
    if base.endswith(('.sto', '.sth')):
        return read_stockholm(path)
    return read_fasta(path)
