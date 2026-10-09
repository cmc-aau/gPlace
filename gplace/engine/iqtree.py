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

"""The complete model of an IQ-TREE fit (place --full_model), at IQ-TREE's full precision.

IQ-TREE writes every fitted parameter to its checkpoint (<prefix>.ckp.gz) with ten significant
digits, and the model's name and rate-model details to <prefix>.iqtree. pplacer reads only part of
such a fit (fitted.py); with --full_model gPlace places with all of it:

  * every reversible DNA model (JC ... GTR, 6-digit codes): the checkpoint's six GTR-form rates;
  * every reversible protein matrix of IQ-TREE (data/iqtree_protein_models.json, from its source);
  * the frequencies IQ-TREE used (model, +F, +FO, +FQ, ...), from the checkpoint;
  * rates: none, +G (mean or median rates, computed as IQ-TREE computes them: yang_gamma.iqtree_rates),
    +I, +I+G, +R and +I+R: rate categories with weights (an invariable-sites category has rate 0);
  * protein mixtures (C10..C60, EX2, EX3, EHO, UL2, UL3, LG4M, LG4X, MIX{...}, <matrix>+C20, ...):
    every component's matrix (IQ-TREE's definitions), frequencies, weight and rate, crossed with the
    rate categories (+G, +R) or paired with them (*G4, *R4).

Results are then IQ-TREE's model, not pplacer's: they differ from pplacer's, which knows no such
models (even LG differs: IQ-TREE normalises LG's published frequencies).
"""

import gzip
import logging
import json
import os
import re

from gplace.engine import core, yang_gamma
from gplace.exceptions import RefpkgMalformed

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'iqtree_protein_models.json')
_PROTEIN = None


_MIX = None


def protein_models():
    """IQ-TREE's reversible protein matrices: {NAME (upper case): (exchangeabilities [190], frequencies [20])}."""
    global _PROTEIN, _MIX
    if _PROTEIN is None:
        with open(DATA) as fh:
            d = json.load(fh)
        _PROTEIN = {k.upper(): ([float(x) for x in v['exchangeabilities']], [float(x) for x in v['frequencies']])
                    for k, v in d['models'].items()}
        _MIX = {k.upper(): v for k, v in d.get('mixtures', {}).items()}
    return _PROTEIN


def mixture_definitions():
    """IQ-TREE's built-in mixture definitions: {NAME (upper case): expression}."""
    protein_models()
    return _MIX


def _exchangeabilities(name, where):
    """The 20 x 20 exchangeability matrix of an IQ-TREE protein matrix (zero diagonal)."""
    table = protein_models().get(name.upper())
    if table is None:
        raise RefpkgMalformed(f'--full_model: {where}: protein matrix {name} is not one of IQ-TREE\'s reversible '
                              f'built-in matrices that gPlace knows.')
    exch, _ = table
    R = [[0.0] * 20 for _ in range(20)]
    k = 0
    for i in range(1, 20):                              # PAML: lower triangle, row by row
        for j in range(i):
            R[i][j] = R[j][i] = exch[k]
            k += 1
    return R


def read_checkpoint(path):
    """The sections of an IQ-TREE checkpoint: {section: {key: value string}} (top-level keys under '')."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', errors='replace') as fh:
        lines = fh.read().split('\n')
    out, section = {'': {}}, ''
    for line in lines:
        if not line.strip() or line.startswith('---'):
            continue
        if not line.startswith(' '):
            key, _, value = line.partition(':')
            if value.strip():
                out[''][key.strip()] = value.strip()
                section = ''
            else:
                section = key.strip()
                out.setdefault(section, {})
        else:
            key, _, value = line.strip().partition(':')
            out.setdefault(section, {})[key.strip()] = value.strip()
    return out


def _floats(s):
    return [float(x) for x in re.split(r'[,\s]+', s.strip()) if x]


def locate(path):
    """(checkpoint, .iqtree) of a fit given either file (or its prefix); either may be None."""
    for ext in ('.ckp.gz', '.iqtree', '.treefile', '.log'):
        if path.endswith(ext):
            path = path[:-len(ext)]
            break
    ckp, iqt = path + '.ckp.gz', path + '.iqtree'
    return (ckp if os.path.isfile(ckp) else None), (iqt if os.path.isfile(iqt) else None)


def read_fit(checkpoint, iqtree_file=None):
    """The fitted model: {'model': name as fitted, 'matrix', 'alphabet' (core.DNA or core.PROTEIN),
    'R' (exchangeabilities, 20 x 20 or 4 x 4, zero diagonal), 'freqs', 'rates' ('none' or 'gamma'),
    'alpha', 'n_cats', 'logl' (IQ-TREE's log-likelihood of the tree), 'program'}."""
    ck = read_checkpoint(checkpoint)
    text = ''
    if iqtree_file:
        with open(iqtree_file, errors='replace') as fh:
            text = fh.read()
    m = re.search(r'^Model of substitution:\s*(\S+)', text, re.MULTILINE)
    model = m.group(1) if m else None
    if model is None:                                   # the command line, when there is no .iqtree
        cmd = ck.get('iqtree', {}).get('command', '')
        m = re.search(r'-m\s+(\S+)', cmd)
        model = m.group(1) if m else None
    if not model or model.upper() in ('MFP', 'TEST', 'TESTONLY', 'MF', 'TESTNEW', 'MFP+MERGE'):
        raise RefpkgMalformed(f'{checkpoint}: no fitted model name; give the .iqtree file next to the checkpoint.')
    version = ck.get('iqtree', {}).get('version', '')
    if any(s.startswith('PartitionModel') for s in ck) or re.search(r'^Partition file name', text, re.MULTILINE):
        raise RefpkgMalformed(f'--full_model: {checkpoint}: partitioned models (IQ-TREE -p/-q/-Q) are not supported yet '
                              f'(see the FAQ).')
    fit = {'model': model, 'program': f'IQ-TREE {version}'.strip(), 'source': checkpoint}
    if any(s.startswith('ModelMixture') for s in ck):
        return _read_mixture(ck, text, fit, checkpoint)
    parts = model.split('+')
    base = parts[0].split('{')[0]
    rest = [p.split('{')[0] for p in parts[1:]]
    unsupported = [p for p in rest if not re.fullmatch(r'G\d*|I|R\d*|F|FO|FQ|FU|FRY|FWS|FMK|F\d{4}', p)]
    if model.upper().startswith(('MIX', 'FMIX')) or unsupported or '*' in model:
        raise RefpkgMalformed(f'--full_model: {model} has parts gPlace does not place with yet: '
                              f'{", ".join(unsupported) or model} (supported: reversible matrices, frequencies, +G, +I, +R).')
    if 'ModelDNA' in ck:
        rates = _floats(ck['ModelDNA']['rates'])
        if len(rates) != 6:
            raise RefpkgMalformed(f'{checkpoint}: {len(rates)} DNA rates, expected 6 (GTR form).')
        freqs = _floats(ck['ModelDNA']['state_freq'])
        R = [[0.0] * 4 for _ in range(4)]
        for (i, j), r in zip(((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)), rates):
            R[i][j] = R[j][i] = r
        fit.update(alphabet=core.DNA, matrix=base, R=R, freqs=freqs)
    elif 'ModelProtein' in ck:
        R = _exchangeabilities(base, checkpoint)
        freqs = _floats(ck['ModelProtein']['state_freq'])
        fit.update(alphabet=core.PROTEIN, matrix=base, R=R, freqs=freqs)
    else:
        raise RefpkgMalformed(f'{checkpoint}: no DNA or protein model in the checkpoint '
                              f'(sections: {", ".join(s for s in ck if s)}).')
    if len(fit['freqs']) != fit['alphabet'].ns:
        raise RefpkgMalformed(f'{checkpoint}: {len(fit["freqs"])} frequencies, expected {fit["alphabet"].ns}.')
    total = sum(fit['freqs'])                           # IQ-TREE rescales them to sum 1 (LG's sum to 1 + 3e-8)
    fit['freqs'] = [f / total for f in fit['freqs']]
    g = next((p for p in rest if re.fullmatch(r'G\d*', p)), None)
    r = next((p for p in rest if re.fullmatch(r'R\d*', p)), None)
    inv = 'I' in rest
    sec = lambda prefix: next((ck[s] for s in ck if s.startswith(prefix)), None)
    fit.update(alpha=None, p_invar=0.0, free=None, median=False)
    if g:
        gamma = sec('RateGamma')
        if gamma is None or 'gamma_shape' not in gamma:
            raise RefpkgMalformed(f'{checkpoint}: {model} but no gamma shape in the checkpoint.')
        fit.update(alpha=float(gamma['gamma_shape']), n_cats=int(g[1:] or 4), median=_median(text, ck))
        if inv:
            fit['p_invar'] = float(gamma['p_invar'])
    elif r:
        free = sec('RateFree')
        if free is None or 'prop' not in free:
            raise RefpkgMalformed(f'{checkpoint}: {model} but no free rates in the checkpoint.')
        fit['free'] = (_floats(free['prop']), _floats(free['rates']))
        fit['n_cats'] = len(fit['free'][0])
        if inv:
            fit['p_invar'] = float(free['p_invar'])
    else:
        fit['n_cats'] = 1
        if inv:
            fit['p_invar'] = float((sec('RateInvar') or {}).get('p_invar', 'nan'))
    fit['rates'] = ('invariant+' if inv else '') + ('gamma' if g else 'free-rate' if r else 'none' if not inv else '')
    fit['rates'] = fit['rates'].rstrip('+') or 'none'
    fit['categories'] = fit['n_cats'] + (1 if inv else 0)
    cs = ck.get('CandidateSet', {})
    fit['logl'] = float(cs['0'].split()[0]) if '0' in cs else None
    return fit


def _median(text, ck):
    """Median gamma rates (--gamma-median): the .iqtree file says so, else the command line in the checkpoint."""
    if 'computed as MEDIAN' in text or 'computed as MEAN' in text:
        return 'computed as MEDIAN' in text
    return bool(re.search(r'(^|\s)(--gamma-median|-gmedian)(\s|$)', ck.get('iqtree', {}).get('command', '')))


def describe(fit):
    if fit.get('components'):
        K = len(fit['components'])
        rates = (f'gamma ({fit["n_cats"]} categories, {"median" if fit.get("median") else "mean"} rates, '
                 f'alpha {fit["alpha"]})' if fit['alpha'] is not None
                 else f'free rates ({fit["n_cats"]} categories)' if fit['free'] else 'no rate variation')
        how = 'each paired with one rate category' if fit['paired'] else 'each under every rate category'
        return (f'{fit["model"]} ({fit["program"]}): mixture of {K} components ({how}), {rates}; '
                f'{fit["categories"]} categories')
    parts = []
    if fit['p_invar']:
        parts.append(f'invariable sites ({fit["p_invar"]})')
    if fit['alpha'] is not None:
        parts.append(f'gamma ({fit["n_cats"]} categories, {"median" if fit.get("median") else "mean"} rates, '
                     f'alpha {fit["alpha"]})')
    if fit['free']:
        parts.append(f'free rates ({fit["n_cats"]} categories)')
    return (f'{fit["model"]} ({fit["program"]}): {fit["matrix"]}, {", ".join(parts) or "no rate variation"}, '
            f'the fitted frequencies')


def categories(fit):
    """(rates, weights) of the fit's rate categories, as IQ-TREE weighs them: +G IQ-TREE's gamma rates
    (yang_gamma.iqtree_rates), weights 1/n; +I a rate-0 category of weight p, the gamma rates then averaging
    1/(1 - p) (RateGammaInvar), weights (1 - p)/n; +I alone rate 1/(1 - p), weight 1 - p; free rates and their
    proportions are the checkpoint's (with +I they already sum to 1 - p and average 1)."""
    p = fit['p_invar']
    if fit['free']:
        props, rates = fit['free']
        rates, weights = list(rates), list(props)
    elif fit['alpha'] is not None:
        n = fit['n_cats']
        rates = yang_gamma.iqtree_rates(n, fit['alpha'], fit.get('median', False), 1.0 / (1.0 - p) if p else 1.0)
        weights = [(1.0 - p) / n] * n
    else:
        rates, weights = [1.0 / (1.0 - p) if p else 1.0], [1.0 - p]
    if p:
        rates, weights = [0.0] + rates, [p] + weights
    return rates, weights


def _profile_weights():
    """{PROFILE: (rate, weight)} of the frequency profiles of IQ-TREE's FMIX definitions (C10pi1, ...)."""
    out = {}
    for expr in mixture_definitions().values():
        m = re.search(r'FMIX\{(.*)\}', expr)
        for entry in (m.group(1).split(',') if m else []):
            f = entry.split(':')
            if len(f) == 3:
                out[f[0].upper()] = (float(f[1]), float(f[2]))
    return out


def _rate_categories(ck, text, model, checkpoint):
    """(rates, weights, alpha, n) of a mixture's rate model: none, +G (gamma, computed by the caller:
    rates None), +R (the checkpoint's)."""
    sec = lambda prefix: next((ck[s] for s in ck if s.startswith(prefix)), None)
    if re.search(r'[+*]I(\W|$)', model) or 'Invar' in (re.search(r'^Model of rate heterogeneity:\s*(.*)$', text,
                                                                     re.MULTILINE) or [None, ''])[1]:
        raise RefpkgMalformed(f'--full_model: {model}: invariable sites (+I) in a mixture are not supported yet.')
    gamma, free = sec('RateGamma'), sec('RateFree')
    if gamma is not None:
        m = re.search(r'Gamma with (\d+) categories', text)
        return None, None, float(gamma['gamma_shape']), int(m.group(1)) if m else 4
    if free is not None:
        props, rates = _floats(free['prop']), _floats(free['rates'])
        return rates, props, None, len(props)
    return [1.0], [1.0], None, 1


def _read_mixture(ck, text, fit, checkpoint):
    """A protein mixture: components (matrix, frequencies, weight, rate) and how they combine with
    the rate categories, from the .iqtree component table, the checkpoint (frequencies, fitted
    weights) and IQ-TREE's definitions (matrices, fixed weights and rates)."""
    mix = next(ck[s] for s in ck if s.startswith('ModelMixture'))
    if any('!ModelDNA!' in k for k in mix):
        raise RefpkgMalformed(f'--full_model: {checkpoint}: DNA mixtures are not supported yet (see the FAQ).')
    m = re.search(r'^Mixture model of substitution:\s*(\S+)', text, re.MULTILINE)
    if not m:
        raise RefpkgMalformed(f'--full_model: {checkpoint}: a mixture needs its .iqtree file (component table).')
    expr = m.group(1)
    rows = re.findall(r'^\s+(\d+)\s+(\S+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+(\S+)\s*$',
                      text[m.end():text.find('Model of rate heterogeneity', m.end())], re.MULTILINE)
    if not rows:
        raise RefpkgMalformed(f'--full_model: {checkpoint}: no component table for the mixture {expr}.')
    K = len(rows)
    freqs = {int(k[len('Component'):].split('!')[0]): _floats(v) for k, v in mix.items()
             if k.startswith('Component') and k.endswith('state_freq')}
    if sorted(freqs) != list(range(1, K + 1)):
        raise RefpkgMalformed(f'--full_model: {checkpoint}: {K} components in the table but frequencies for '
                              f'{sorted(freqs)} in the checkpoint.')
    paired = '*' in expr
    # the definition: a built-in name (EX2, C10, LG4M*G4, ...) or an explicit MIX{...}/FMIX{...}
    head = re.split(r'[+*]', expr)[0] if not expr.startswith(('MIX{', 'FMIX{')) else None
    defn = mixture_definitions().get(head.upper()) if head else None
    if defn is None and head and '+' in expr:                    # <matrix>+C20...: the profile part's definition
        for part in expr.split('+')[1:]:
            defn = mixture_definitions().get(re.split(r'[*{]', part)[0].upper())
            if defn:
                break
    items = re.search(r'F?MIX\{(.*)\}', defn or expr)
    def_rate, def_weight = [1.0] * K, None
    if items:
        entries = items.group(1).split(',')
        if len(entries) == K:
            fields = [e.split(':') for e in entries]
            def_rate = [float(f[1]) if len(f) > 1 and f[1] else 1.0 for f in fields]
            if all(len(f) > 2 for f in fields):
                def_weight = [float(f[2]) for f in fields]
    profiles = [re.search(r'\+F([A-Za-z0-9_.]+)$', r[4]) for r in rows]   # e.g. LG+FC10pi1: profile C10pi1
    fixed = _profile_weights()
    if 'prop' in mix:
        weights = _floats(mix['prop'])
    elif def_weight is not None:
        weights = def_weight
    elif paired:
        weights = [1.0] * K
    elif all(p and p.group(1).upper() in fixed for p in profiles):   # <matrix>+C10...: the profiles' own weights
        def_rate = [fixed[p.group(1).upper()][0] for p in profiles]
        weights = [fixed[p.group(1).upper()][1] for p in profiles]
    else:
        weights = [float(r[3]) for r in rows]
        logging.getLogger('timestamp').warning(
            f'--full_model: {expr}: mixture weights as printed in the .iqtree file (4 decimals).')
    total = sum(weights)
    weights = [w / total for w in weights] if not paired else weights
    scale = sum(w * r for w, r in zip(weights, def_rate)) if not paired else 1.0
    comp_rate = [r / scale for r in def_rate]                    # IQ-TREE: mean component rate 1
    comps = []
    for (no, name, rate, weight, params), k in zip(rows, range(1, K + 1)):
        matrix = re.split(r'\+F', params)[0]
        total_f = sum(freqs[k])
        comps.append({'name': name, 'matrix': matrix, 'R': _exchangeabilities(matrix, checkpoint),
                      'freqs': [f / total_f for f in freqs[k]], 'weight': weights[k - 1],
                      'rate': comp_rate[k - 1]})
    rates, rweights, alpha, n = _rate_categories(ck, text, expr, checkpoint)
    if paired and n != K:
        raise RefpkgMalformed(f'--full_model: {expr}: {K} components paired with {n} rate categories.')
    fit.update(alphabet=core.PROTEIN, matrix=expr, R=comps[0]['R'], freqs=comps[0]['freqs'], components=comps,
               paired=paired, alpha=alpha, n_cats=n, p_invar=0.0, median=alpha is not None and _median(text, ck),
               free=None if rates is None or n == 1 else (rweights, rates),
               rates='gamma' if alpha is not None else 'free-rate' if n > 1 else 'none',
               categories=K * (1 if paired else n))
    cs = ck.get('CandidateSet', {})
    fit['logl'] = float(cs['0'].split()[0]) if '0' in cs else None
    return fit


def mixture_categories(fit):
    """(rates, weights, component index) of every category of a mixture, component by component."""
    if fit['alpha'] is not None:
        g_rates = yang_gamma.iqtree_rates(fit['n_cats'], fit['alpha'], fit.get('median', False))
        g_weights = [1.0 / fit['n_cats']] * fit['n_cats']
    else:
        g_weights, g_rates = fit['free'] if fit['free'] else ([1.0], [1.0])
    rates, weights, comp = [], [], []
    for k, c in enumerate(fit['components']):
        if fit['paired']:
            rates.append(c['rate'] * g_rates[k]); weights.append(g_weights[k]); comp.append(k)
        else:
            for r, w in zip(g_rates, g_weights):
                rates.append(c['rate'] * r); weights.append(c['weight'] * w); comp.append(k)
    return rates, weights, comp
