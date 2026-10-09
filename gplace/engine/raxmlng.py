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

"""The complete model of a RAxML-NG fit (place --full_model), as RAxML-NG evaluates it.

RAxML-NG writes the fitted model to <prefix>.raxml.bestModel: its model string with every estimated
parameter written as a value (6 decimals, or --precision), e.g. GTR{1.27/3.38/.../1}+FU{0.41/...}+G4m{0.78}.
That file is RAxML-NG's own input (`raxml-ng --model <prefix>.raxml.bestModel`), and gPlace reads it the
way RAxML-NG does:

  * DNA models JC ... GTR (and DNA<symmetry> strings): the rates, normalised by the last one;
  * protein matrices: RAxML-NG's own tables (data/raxmlng_protein_models.json, from its source), or the
    fit's rates (PROTGTR{...}, <matrix>{...});
  * frequencies: the model's, +FU{...}, +FE, +FC (counted in the reference alignment, as RAxML-NG counts
    them); RAxML-NG rescales them to sum 1 when they are off by more than 1e-8;
  * +G<n>m / +G<n>a: RAxML-NG's discrete gamma (Yang's algorithm), mean or median rates;
  * +IU{p} / +IC (proportion counted in the alignment): a rate-0 category of weight p, the other rates
    divided by 1 - p. RAxML-NG counts a site as invariable only when one state fits every sequence;
  * +R<n>{rates}{weights}, normalised as RAxML-NG does (weights sum 1, mean rate 1);
  * LG4M and LG4X: four matrices, each with its own frequencies, paired with the four rate categories.

Results are then RAxML-NG's model, not pplacer's (even LG and +G4 differ in the last digits: other tables,
other gamma algorithm).
"""

import json
import logging
import os
import re

import numpy as np

from gplace.engine import core, yang_gamma
from gplace.exceptions import RefpkgMalformed

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'raxmlng_protein_models.json')
EPSILON = 1e-8                         # coraxlib's CORAX_MISC_EPSILON
MIN_FREQ = 1e-6                        # CORAX_EIGEN_MINFREQ: states at or below it leave the rate matrix
DNA_MODELS = {                         # rate symmetries (AC AG AT CG CT GT) and equal frequencies, as in coraxlib
    'JC': ('000000', True), 'K80': ('010010', True), 'F81': ('000000', False), 'HKY': ('010010', False),
    'TN93EF': ('010020', True), 'TN93': ('010020', False), 'K81': ('012210', True), 'K81UF': ('012210', False),
    'TPM2': ('010212', True), 'TPM2UF': ('010212', False), 'TPM3': ('012012', True), 'TPM3UF': ('012012', False),
    'TIM1': ('012230', True), 'TIM1UF': ('012230', False), 'TIM2': ('010232', True), 'TIM2UF': ('010232', False),
    'TIM3': ('012032', True), 'TIM3UF': ('012032', False), 'TVMEF': ('012314', True), 'TVM': ('012314', False),
    'SYM': ('012345', True), 'GTR': ('012345', False)}
DNA_ALIASES = {'TRNEF': 'TN93EF', 'TRN': 'TN93', 'TPM1': 'K81', 'TPM1UF': 'K81UF', 'TPM2EF': 'TPM2',
               'TPM3EF': 'TPM3', 'TIM1EF': 'TIM1', 'TIM2EF': 'TIM2', 'TIM3EF': 'TIM3'}
DNA_FIXED_RATES = {'JC', 'F81'}        # rates given by the model (equal), not fitted
_PROTEIN = None
_MIX = None


def protein_models():
    """RAxML-NG's protein matrices: {NAME (upper case): (rates [190, upper triangle], frequencies [20])}."""
    global _PROTEIN, _MIX
    if _PROTEIN is None:
        with open(DATA) as fh:
            d = json.load(fh)
        _PROTEIN = {k.upper(): ([float(x) for x in v['exchangeabilities']], [float(x) for x in v['frequencies']])
                    for k, v in d['models'].items()}
        _MIX = {k.upper(): v for k, v in d['mixtures'].items()}
    return _PROTEIN


def mixtures():
    """RAxML-NG's protein mixtures: {NAME: {'components': [matrix names], 'rates': 'gamma' or 'free'}}."""
    protein_models()
    return _MIX


# ------------------------------------------------------------------------------------------ gamma --
def gamma_rates(n, alpha, median=False):
    """RAxML-NG's rates of n equally likely gamma categories (mean, or median rescaled to mean 1)."""
    return yang_gamma.raxmlng_rates(n, alpha, median)


# ------------------------------------------------------------------------------------- the files --
def locate(path):
    """(bestModel, log) of a fit given any of its files (or its prefix); either may be None."""
    for ext in ('.raxml.bestModel', '.raxml.log', '.raxml.bestTree', '.raxml.rba', '.raxml.startTree'):
        if path.endswith(ext):
            path = path[:-len(ext)]
            break
    best, log = path + '.raxml.bestModel', path + '.raxml.log'
    return (best if os.path.isfile(best) else None), (log if os.path.isfile(log) else None)


def _values(s):
    return [float(x) for x in re.split(r'[/,\s]+', s.strip()) if x]


def parse_model(text, where):
    """A RAxML-NG model string (as RAxML-NG reads it): {'name', 'rates', 'freqs' (mode, values), 'pinv'
    (mode, value), 'gamma' (n, median, alpha), 'free' (n, rates, weights), 'other' (unsupported parts)}."""
    m = re.match(r'([^+{\[]*)', text)
    out = {'name': m.group(1).strip(), 'rates': None, 'freqs': None, 'pinv': None, 'gamma': None, 'free': None,
           'other': []}
    s = text[m.end():]
    param = re.compile(r'[{\[]([^}\]]*)[}\]]')
    p = param.match(s)
    if p:
        out['rates'] = _values(p.group(1))
        s = s[p.end():]
    for token in re.findall(r'\+([^+]*)', s):
        head, rest = token[:1].upper(), token[1:]
        params = [_values(v) if re.fullmatch(r'[-+0-9.eE/,\s]*', v) else v for v in param.findall(rest)]
        bare = param.sub('', rest).upper()
        if head == 'F' and (bare in ('', 'C', 'O', 'E', 'U') or not bare):
            mode = 'U' if params else (bare or 'C')
            out['freqs'] = (mode, params[0] if params else None)
        elif head == 'I' and bare in ('', 'O', 'C', 'U'):
            mode = 'U' if params else (bare or 'O')
            out['pinv'] = (mode, params[0][0] if params else None)
        elif head == 'G' and re.fullmatch(r'\d*[AM]?', bare):
            n = int(re.match(r'\d*', bare).group(0) or 4)
            out['gamma'] = (n, bare.endswith('A'), params[0][0] if params else None)
        elif head == 'R' and re.fullmatch(r'\d*', bare):
            n = int(bare or 4)
            out['free'] = (n, params[0] if params else None, params[1] if len(params) > 1 else None)
        elif head == 'B' and bare in ('U', '') and params and params[0] == [1.0]:
            pass                                        # branch length scaler 1
        else:
            out['other'].append('+' + token)
    if any(isinstance(v, str) for v in (out['freqs'] or (None, None))[1:]) or isinstance(out['rates'], str):
        out['other'].append('parameters from a file')
    return out


def _read_log(log):
    """(program with version, final log-likelihood, printed empirical frequencies, printed empirical p-inv,
    number of alignment columns RAxML-NG used)."""
    if not log:
        return 'RAxML-NG', None, None, None, None
    with open(log, errors='replace') as fh:
        text = fh.read()
    v = re.search(r'^RAxML-NG v\. (\S+)', text, re.MULTILINE)
    logl = re.findall(r'^Final LogLikelihood:\s*(\S+)', text, re.MULTILINE)
    freqs = re.findall(r'Base frequencies \(empirical\):\s*([0-9.eE\s+-]+?)\s*$', text, re.MULTILINE)
    pinv = re.findall(r'P-inv \(empirical\):\s*(\S+)', text)
    sites = re.findall(r'^Alignment sites / patterns:\s*(\d+)', text, re.MULTILINE)
    return (f'RAxML-NG {v.group(1)}' if v else 'RAxML-NG', float(logl[-1]) if logl else None,
            _values(freqs[-1]) if freqs else None, float(pinv[-1]) if pinv else None,
            int(sites[-1]) if sites else None)


def _set_freqs(freqs):
    """coraxlib's corax_set_frequencies: rescaled to sum 1 when the sum is off by more than 1e-8."""
    total = 0.0
    for f in freqs:
        total += f
    return [f / total for f in freqs] if abs(total - 1.0) > EPSILON else list(freqs)


def _symmetric(rates, ns):
    """The ns x ns exchangeability matrix of rates in coraxlib's order (upper triangle, row by row)."""
    R = [[0.0] * ns for _ in range(ns)]
    k = 0
    for i in range(ns):
        for j in range(i + 1, ns):
            R[i][j] = R[j][i] = rates[k]
            k += 1
    return R


def _user_rates(values, sym, where):
    """Model::set_user_srates: the unique rates divided by the last fitted one, spread over the symmetry."""
    n_uniq = (max(sym) + 1) if sym else len(values)
    if len(values) != n_uniq:
        raise RefpkgMalformed(f'--full_model: {where}: {len(values)} substitution rates, expected {n_uniq}.')
    last = values[sym[-1]] if sym else values[-1]
    if not last > 0:
        raise RefpkgMalformed(f'--full_model: {where}: the last substitution rate is {last}.')
    values = [r / last for r in values]
    return [values[k] for k in sym] if sym else values


def read_fit(best_model, log=None):
    """The fitted model: {'model': the model string, 'program', 'alphabet', 'matrix', 'R', 'freqs' (None when
    counted from the alignment: freqs_mode 'empirical'), 'alpha', 'median', 'n_cats', 'p_invar' (None: counted,
    p_invar_mode 'empirical'), 'free' (weights, rates), 'components' (LG4M, LG4X), 'paired', 'rates',
    'categories', 'logl' (RAxML-NG's log-likelihood of the tree, from the log), 'sites' (partition length)}."""
    with open(best_model, errors='replace') as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    if not lines:
        raise RefpkgMalformed(f'--full_model: {best_model} is empty.')
    if len(lines) > 1:
        raise RefpkgMalformed(f'--full_model: {best_model}: partitioned models ({len(lines)} partitions) are not '
                              f'supported yet (see the FAQ).')
    model, _, part = lines[0].partition(',')
    model = model.strip()
    rng = re.fullmatch(r'\s*[^=]*=\s*1\s*-\s*(\d+)\s*', part)
    if part.strip() and not rng:
        raise RefpkgMalformed(f'--full_model: {best_model}: the model covers the columns "{part.split("=")[-1].strip()}"'
                              f'; gPlace places with one model over the whole alignment (partitions are not '
                              f'supported yet).')
    program, logl, log_freqs, log_pinv, log_sites = _read_log(log)
    pm = parse_model(model, best_model)
    if pm['other']:
        raise RefpkgMalformed(f'--full_model: {model} has parts gPlace does not place with yet: '
                              f'{", ".join(pm["other"])} (supported: DNA and protein matrices, frequencies, +G, +I, +R).')
    name = pm['name'].upper()
    fit = {'model': model, 'program': program, 'source': best_model, 'logl': logl,
           'sites': int(rng.group(1)) if rng else None, 'log_freqs': log_freqs, 'log_pinv': log_pinv,
           'log_sites': log_sites,
           'components': None, 'paired': False, 'median': False}
    proteins = protein_models()
    mix = mixtures().get(name)
    sym = None
    if mix or name in proteins or name == 'PROTGTR':
        fit['alphabet'], ns = core.PROTEIN, 20
        if mix:
            if pm['rates'] is not None:
                raise RefpkgMalformed(f'--full_model: {model}: user rates for a mixture.')
            comps = [proteins[c.upper()] for c in mix['components']]
        elif pm['rates'] is not None:
            comps = [(_user_rates(pm['rates'], None, model), (proteins.get(name) or (None, None))[1])]
        elif name == 'PROTGTR':
            raise RefpkgMalformed(f'--full_model: {model}: PROTGTR without its rates; give the .raxml.bestModel file.')
        else:
            comps = [proteins[name]]
    else:
        name = DNA_ALIASES.get(name, name)
        if name in DNA_MODELS:
            sym_str, equal = DNA_MODELS[name]
        elif name.startswith('DNA') and len(name) == 9:
            sym_str, equal = name[3:], False
        elif name.startswith('PROT'):
            raise RefpkgMalformed(f'--full_model: {model}: protein models with custom rate symmetries are not '
                                  f'supported yet.')
        elif name.startswith(('BIN', 'MULTI', 'GT', 'GP')):
            raise RefpkgMalformed(f'--full_model: {model}: binary, multi-state and genotype data are not supported '
                                  f'yet (gPlace places DNA and protein).')
        else:
            raise RefpkgMalformed(f'--full_model: {model}: unknown RAxML-NG model "{pm["name"]}".')
        low = min(sym_str)
        sym = [ord(c) - ord(low) for c in sym_str]
        fit['alphabet'], ns = core.DNA, 4
        if pm['rates'] is not None:
            rates = _user_rates(pm['rates'], sym, model)
        elif name in DNA_FIXED_RATES:
            rates = [1.0] * 6
        else:
            raise RefpkgMalformed(f'--full_model: {model}: no substitution rates; give the .raxml.bestModel file '
                                  f'(the fitted rates are written there).')
        comps = [(rates, [0.25] * 4 if equal else None)]
    # frequencies: the option, else the model's
    mode, values = pm['freqs'] or ('model', None)
    if mode == 'U':
        if len(values) != ns:
            raise RefpkgMalformed(f'--full_model: {model}: {len(values)} frequencies, expected {ns}.')
        comps = [(r, values) for r, _ in comps]
    elif mode == 'E':
        comps = [(r, [1.0 / ns] * ns) for r, _ in comps]
    elif mode == 'C':
        comps = [(r, None) for r, _ in comps]
    elif mode == 'O':
        raise RefpkgMalformed(f'--full_model: {model}: ML frequencies without their values; give the '
                              f'.raxml.bestModel file.')
    elif any(f is None for _, f in comps):
        raise RefpkgMalformed(f'--full_model: {model}: no frequencies; give the .raxml.bestModel file.')
    fit['freqs_mode'] = 'empirical' if mode == 'C' else 'given'
    fit['matrices'] = [(_symmetric(r, ns), None if f is None else _set_freqs(f)) for r, f in comps]
    fit['matrix'] = pm['name'] + ('{...}' if pm['rates'] is not None else '')
    # rates
    pinv = pm['pinv']
    if pinv and pinv[0] == 'O':
        raise RefpkgMalformed(f'--full_model: {model}: +I without its value; give the .raxml.bestModel file.')
    fit['p_invar_mode'] = 'empirical' if pinv and pinv[0] == 'C' else 'given'
    fit['p_invar'] = None if pinv and pinv[0] == 'C' else (pinv[1] if pinv else 0.0)
    fit.update(alpha=None, free=None, n_cats=1)
    rate_het = mix['rates'] if mix else None
    if pm['gamma']:
        n, median, alpha = pm['gamma']
        if alpha is None:
            raise RefpkgMalformed(f'--full_model: {model}: +G without alpha; give the .raxml.bestModel file.')
        fit.update(alpha=alpha, median=median, n_cats=n)
    elif pm['free']:
        n, rates, weights = pm['free']
        if rates is None:
            raise RefpkgMalformed(f'--full_model: {model}: +R without its rates; give the .raxml.bestModel file.')
        if len(rates) != n or (weights is not None and len(weights) != n):
            raise RefpkgMalformed(f'--full_model: {model}: {n} free-rate categories but {len(rates)} rates and '
                                  f'{len(weights or [])} weights.')
        if weights is None:
            weights = [1.0 / n] * n
        else:                                           # Model::init_model_opts: weights sum 1, mean rate 1
            total = 0.0
            for w in weights:
                total += w
            weights = [w / total for w in weights]
        mean = 0.0
        for r, w in zip(rates, weights):
            mean += r * w
        fit.update(free=(weights, [r / mean for r in rates]), n_cats=n)
    elif rate_het:
        raise RefpkgMalformed(f'--full_model: {model}: the mixture\'s {rate_het} rates are missing; give the '
                              f'.raxml.bestModel file.')
    if mix:
        K = len(mix['components'])
        if fit['n_cats'] != K:
            raise RefpkgMalformed(f'--full_model: {model}: {K} matrices with {fit["n_cats"]} rate categories.')
        if fit['p_invar'] != 0.0:
            raise RefpkgMalformed(f'--full_model: {model}: invariable sites (+I) in a mixture are not supported yet.')
        fit['components'] = [{'name': c, 'matrix': c} for c in mix['components']]
        fit['paired'] = True
    fit['rates'] = (('invariant+' if fit['p_invar'] != 0.0 else '')
                    + ('gamma' if fit['alpha'] is not None else 'free-rate' if fit['free'] else '')).strip('+') or 'none'
    fit['categories'] = fit['n_cats'] + (1 if fit['p_invar'] != 0.0 else 0)
    return fit


# -------------------------------------------------------------------- what needs the alignment --
_MASKS = {}


def _state_masks(alphabet):
    """coraxlib's character map (corax_map_aa / corax_map_nt) as bit masks of states by byte; 0: not a
    character RAxML-NG reads, all bits: a gap."""
    if alphabet.name not in _MASKS:
        ns = alphabet.ns
        mask = np.zeros(256, np.uint32)
        if alphabet is core.PROTEIN:
            for i, c in enumerate(core.PROT_CODE):
                mask[ord(c)] = 1 << i
            idx = core.PROT_CODE.index
            for c, states in (('B', 'ND'), ('Z', 'QE'), ('J', 'IL')):
                mask[ord(c)] = sum(1 << idx(s) for s in states)
            gaps = 'X-?.*'
        else:
            for c, v in core.NUC_MAP.items():
                mask[ord(c)] = sum(1 << i for i in range(4) if v[i])
            mask[ord('U')] = mask[ord('T')]
            gaps = 'NOX-?.'
        for c in gaps:
            mask[ord(c)] = (1 << ns) - 1
        for c in range(ord('A'), ord('Z') + 1):         # lower case as upper case
            mask[c + 32] = mask[c]
        _MASKS[alphabet.name] = mask
    return _MASKS[alphabet.name]


def alignment_stats(ref_align, alphabet, gap_columns=False):
    """(empirical frequencies, proportion of invariable columns) of the reference alignment, as coraxlib's
    corax_msa_compute_stats counts them: gaps (every state) left out, an ambiguous character shared among its
    states; a column is invariable when exactly one state fits every sequence. RAxML-NG first drops the
    columns that are gaps in every sequence (unless run with --force): gap_columns=True keeps them."""
    ns, full = alphabet.ns, (1 << alphabet.ns) - 1
    lut = _state_masks(alphabet)
    A = np.frombuffer(''.join(s for _, s in ref_align).encode('latin-1'), np.uint8).reshape(len(ref_align), -1)
    masks = lut[A]
    if not masks.all():
        i, j = map(int, np.argwhere(masks == 0)[0])
        raise RefpkgMalformed(f'--full_model: character "{chr(A[i, j])}" of {ref_align[i][0]} (column {j + 1}) is '
                              f'not one RAxML-NG reads; the reference alignment is not the one RAxML-NG fitted.')
    codes, numbers = np.unique(A[masks != full], return_counts=True)
    counts = [0.0] * ns
    for code, number in zip(codes, numbers):
        m = int(lut[code])
        share = 1.0 / bin(m).count('1')
        for k in range(ns):
            if m >> k & 1:
                counts[k] += int(number) * share
    known = int(numbers.sum())
    common = np.bitwise_and.reduce(masks, axis=0)
    single = (common != 0) & ((common & (common - 1)) == 0)
    columns = A.shape[1] if gap_columns else int((common != full).sum())
    return [c / known for c in counts], int(single.sum()) / float(columns)


def complete(fit, ref_align):
    """The fit with what RAxML-NG counts in the alignment (+FC frequencies, +IC proportion), checked against
    RAxML-NG's log; and the alignment length checked against the model's columns."""
    if fit.get('_complete'):
        return fit
    fit = dict(fit)
    log = logging.getLogger('timestamp')
    if fit['sites'] and ref_align and len(ref_align[0][1]) != fit['sites']:
        raise RefpkgMalformed(f'--full_model: {fit["source"]}: the model covers columns 1-{fit["sites"]}, the '
                              f'reference alignment has {len(ref_align[0][1])}.')
    if fit['freqs_mode'] == 'empirical' or fit['p_invar'] is None:
        freqs, pinv = alignment_stats(ref_align, fit['alphabet'],       # the log says how many columns it used
                                      gap_columns=fit['log_sites'] == len(ref_align[0][1]))
        if fit['freqs_mode'] == 'empirical':
            fit['matrices'] = [(R, _set_freqs(freqs)) for R, _ in fit['matrices']]
            if fit['log_freqs'] and max(abs(a - b) for a, b in zip(freqs, fit['log_freqs'])) > 6e-7:
                log.warning(f'--full_model: the frequencies counted in the reference alignment differ from those '
                            f'RAxML-NG counted (its log): is it the alignment RAxML-NG fitted?')
        if fit['p_invar'] is None:
            fit['p_invar'] = pinv
            if fit['log_pinv'] is not None and abs(pinv - fit['log_pinv']) > 6e-7:
                log.warning(f'--full_model: the proportion of invariable columns of the reference alignment ({pinv}) '
                            f'differs from RAxML-NG\'s ({fit["log_pinv"]}): is it the alignment RAxML-NG fitted?')
    low = [k for _, f in fit['matrices'] for k, x in enumerate(f) if x <= MIN_FREQ]
    if low:
        raise RefpkgMalformed(f'--full_model: {fit["model"]}: frequencies at or below {MIN_FREQ} (states '
                              f'{", ".join(fit["alphabet"].codes[k] for k in sorted(set(low)))}), which RAxML-NG '
                              f'leaves out of the rate matrix; not supported yet.')
    fit['R'], fit['freqs'] = fit['matrices'][0]
    if fit['components']:
        fit['components'] = [dict(c, R=R, freqs=f, weight=None, rate=1.0)
                             for c, (R, f) in zip(fit['components'], fit['matrices'])]
    fit['_complete'] = True
    return fit


def categories(fit):
    """(rates, weights, component of each category or None) as RAxML-NG weighs them: +I a rate-0 category of
    weight p (the others' weights times 1 - p, their rates divided by 1 - p when p > 1e-8)."""
    p = fit['p_invar']
    if fit['free']:
        weights, rates = list(fit['free'][0]), list(fit['free'][1])
    elif fit['alpha'] is not None:
        n = fit['n_cats']
        rates, weights = gamma_rates(n, fit['alpha'], fit['median']), [1.0 / n] * n
    else:
        rates, weights = [1.0], [1.0]
    comp = list(range(len(rates))) if fit['components'] else None
    if p:
        if p > EPSILON:
            rates = [r / (1.0 - p) for r in rates]
        rates, weights = [0.0] + rates, [p] + [w * (1.0 - p) for w in weights]
    return rates, weights, comp


def describe(fit):
    parts = []
    if fit['components']:
        parts.append(f'{len(fit["components"])} matrices, each with its own frequencies and rate category')
    if fit['p_invar'] is None or fit['p_invar']:
        parts.append('invariable sites' + (f' ({fit["p_invar"]})' if fit['p_invar'] else ' (counted)'))
    if fit['alpha'] is not None:
        parts.append(f'gamma ({fit["n_cats"]} categories, {"median" if fit["median"] else "mean"} rates, '
                     f'alpha {fit["alpha"]})')
    if fit['free']:
        parts.append(f'free rates ({fit["n_cats"]} categories)')
    freqs = 'frequencies counted from the reference alignment' if fit['freqs_mode'] == 'empirical' else \
        'the fitted frequencies'
    return f'{fit["model"]} ({fit["program"]}): {", ".join(parts) or "no rate variation"}, {freqs}'
