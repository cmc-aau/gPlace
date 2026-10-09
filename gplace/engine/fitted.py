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

"""Models fitted by IQ-TREE or RAxML-NG, compared with the model pplacer places with.

`taxit create --tree-stats` reads an IQ-TREE .iqtree file or a RAxML-NG log into the package's
phylo_model and keeps the file as tree_stats. pplacer takes only the substitution matrix, a gamma
model (categories, alpha), GTR rates and an empirical-frequencies flag from it; it drops the rest of
the fitted model without a word (+I, free rates, ML-estimated frequencies, ...) and stops with
"Not_found" when no gamma model is left. gPlace places exactly as pplacer does; this module says
what was dropped (warnings) and why a package cannot be used (errors, in the cases where pplacer
stops too).
"""

import logging
import re

from gplace.exceptions import RefpkgMalformed

PPLACER_MATRICES = {'WAG', 'LG', 'JTT', 'GTR'}
_warned = set()


def fit_program(path):
    """'IQ-TREE' or 'RAxML-NG' when path is such a file (as taxit guesses it), else None."""
    try:
        with open(path, errors='replace') as fh:
            for line in fh:
                if line.startswith('RAxML-NG'):
                    return 'RAxML-NG'
                if line.startswith('IQ-TREE'):
                    return 'IQ-TREE'
                if (line.startswith('FastTree') or line.startswith('This is RAxML')
                        or line.startswith('You are using RAxML') or 'PhyML' in line):
                    return None
    except OSError:
        return None
    return None


def _search(rx, text):
    m = re.search(rx, text, re.MULTILINE)
    return m.group(1).strip() if m else None


def read_fit(path):
    """The model an IQ-TREE .iqtree file or a RAxML-NG log reports, as a dict: program (with
    version), model (as written), matrix, freqs ('model', 'empirical', 'ML', 'equal' or the file's
    words), rates ('none', 'gamma', 'invariant', 'invariant+gamma', 'free-rate',
    'invariant+free-rate' or the file's words), n_cats, median (gamma categories by median rates),
    partitions. None for other files."""
    program = fit_program(path)
    if program is None:
        return None
    with open(path, errors='replace') as fh:
        text = fh.read()
    fit = {'partitions': 1, 'median': False, 'n_cats': None}
    if program == 'IQ-TREE':
        fit['program'] = _search(r'^(IQ-TREE [^\s]+)', text) or 'IQ-TREE'
        fit['model'] = _search(r'^Model of substitution:\s*(\S+)', text)
        f = (_search(r'^State frequencies:\s*\((.*)\)', text) or '').lower()
        fit['freqs'] = ('empirical' if 'empirical' in f else 'ML' if 'maximum likelihood' in f or 'estimated' in f
                        else 'equal' if 'equal' in f else 'model' if f == 'model' else f or 'model')
        r = _search(r'^Model of rate heterogeneity:\s*(.*)$', text) or ''
        n = re.search(r'with (\d+) categories', r)
        fit['n_cats'] = int(n.group(1)) if n else None
        low = r.lower()
        inv = low.startswith('invar')
        fit['rates'] = ('none' if low.startswith('uniform') or not low
                        else 'invariant+gamma' if inv and 'gamma' in low
                        else 'invariant+free-rate' if inv and 'freerate' in low
                        else 'invariant' if inv else 'gamma' if 'gamma' in low
                        else 'free-rate' if 'freerate' in low else r)
        fit['median'] = 'computed as MEDIAN' in text
        fit['partitions'] = 2 if re.search(r'^\s*Partition model|^\s*Edge-linked|^\s*Edge-unlinked', text,
                                           re.MULTILINE | re.IGNORECASE) else 1
    else:
        fit['program'] = _search(r'^(RAxML-NG v\. [^\s]+)', text) or 'RAxML-NG'
        models = re.findall(r'^Model: (\S+)', text, re.MULTILINE)
        fit['model'] = models[0] if models else None
        fit['partitions'] = max(1, len(re.findall(r'^Partition \d+:', text, re.MULTILINE)))
        f = (_search(r'Base frequencies \(([^)]*)\)', text) or '').lower()
        fit['freqs'] = {'empirical': 'empirical', 'ml': 'ML', 'model': 'model', 'equal': 'equal'}.get(f, f or 'model')
        parts = (fit['model'] or '').split('+')[1:]
        g = next((p for p in parts if re.fullmatch(r'G\d*[ma]?(\{.*\})?', p)), None)
        inv = any(re.fullmatch(r'I[OCU]?(\{.*\})?', p) for p in parts)
        rr = next((p for p in parts if re.fullmatch(r'R\d*(\{.*\})?', p)), None)
        if g:
            n = re.match(r'G(\d*)', g).group(1)
            fit['n_cats'] = int(n) if n else 4
            fit['median'] = bool(re.match(r'G\d*a', g))
        elif rr:
            n = re.match(r'R(\d*)', rr).group(1)
            fit['n_cats'] = int(n) if n else 4
        fit['rates'] = ('invariant+gamma' if inv and g else 'invariant+free-rate' if inv and rr
                        else 'invariant' if inv else 'gamma' if g else 'free-rate' if rr else 'none')
    fit['matrix'] = (fit['model'] or '').split('+')[0].split('{')[0] or None
    return fit


def _describe_pplacer(pm):
    """The model pplacer builds from a phylo_model, in words."""
    subs = pm.get('subs_model')
    g = pm.get('gamma') or {}
    rates = (f"gamma ({g.get('n_cats')} categories, alpha {g.get('alpha')})" if pm.get('ras_model') == 'gamma'
             else 'CAT' if pm.get('ras_model') == 'Price-CAT' else 'no rate model')
    if subs == 'GTR':
        freqs = 'base frequencies counted from the reference alignment'
    elif pm.get('empirical_frequencies'):
        freqs = 'amino-acid frequencies counted from the reference alignment'
    else:
        freqs = f"{subs}'s own amino-acid frequencies"
    return f'{subs} + {rates}, {freqs}'


def dropped(fit, pm):
    """The parts of a fitted model pplacer does not place with (empty when it uses all of it)."""
    out = []
    matrix, subs = fit.get('matrix'), pm.get('subs_model')
    rates = fit.get('rates')
    if fit.get('partitions', 1) > 1:
        out.append('the partitions (pplacer uses one model for the whole alignment)')
    if 'invariant' in rates:
        out.append('the invariant sites (+I)')
    if 'free-rate' in rates:
        out.append('the free rates (+R)')
    if rates not in ('none', 'gamma', 'invariant', 'invariant+gamma', 'free-rate', 'invariant+free-rate'):
        out.append(f'the rate model ({rates})')
    if fit.get('median') and pm.get('ras_model') == 'gamma':
        out.append('the median gamma rates (pplacer uses mean rates)')
    counted = subs == 'GTR' or bool(pm.get('empirical_frequencies'))
    f = fit.get('freqs')
    if f == 'ML':
        out.append('the ML-estimated frequencies (pplacer ' + ('counts them from the reference alignment)'
                                                                if counted else f"uses {subs}'s own)"))
    elif f == 'empirical' and not counted:
        out.append(f"the frequencies counted from the alignment (pplacer uses {subs}'s own: taxit records "
                   f"{fit['program'].split()[0]}'s counted frequencies as model frequencies)")
    elif f == 'model' and counted and subs != 'GTR':
        out.append(f"{matrix}'s own frequencies (pplacer counts them from the reference alignment)")
    elif f not in ('ML', 'empirical', 'model'):
        out.append(f'the {f} frequencies (pplacer ' + ('counts them from the reference alignment)' if counted
                                                       else f"uses {subs}'s own)"))
    return out


def _refit_hint(fit, pm):
    """Why pplacer cannot use the package, for a fit that pplacer's model reading rejects."""
    matrix, rates, model = fit.get('matrix'), fit.get('rates'), fit.get('model')
    if matrix not in PPLACER_MATRICES:
        if matrix in ('HKY', 'K80', 'K2P', 'F81', 'JC', 'JC69', 'TN93', 'TN', 'SYM', 'TIM', 'TVM', 'K81', 'TPM2',
                      'TPM3', 'TIM2', 'TIM3', 'TPM2u', 'TPM3u', 'TIM2e', 'TIM3e', 'TVMe', 'K3P', 'K3Pu', 'TNe'):
            return (f'{matrix} is a special case of GTR, which pplacer has: refit with GTR+G4 '
                    f'({"GTR+F+G4" if fit["program"].startswith("IQ") else "GTR+G4"}).')
        return 'pplacer knows WAG, LG and JTT (protein) and GTR (DNA): refit with one of them, with +G4.'
    if rates == 'none':
        return f'it has no rate variation across sites, and pplacer needs a gamma model: refit with {matrix}+G4.'
    if 'free-rate' in rates:
        return f'pplacer has no free-rate model (+R): refit with {matrix}+G4.'
    if rates == 'invariant+gamma':
        return (f'pplacer has no invariant-sites model (+I), and taxit records no rate model for +I+G: '
                f'refit with {matrix}+G4.')
    if rates == 'invariant':
        return f'pplacer has no invariant-sites model (+I) and needs a gamma model: refit with {matrix}+G4.'
    return f'refit with {matrix}+G4.'


def check(pm, stats_path, subject, validate):
    """Explain a model against the fit it was made from (stats_path; nothing for other files);
    subject names the tree ("The tree of <package>"). validate(pm) is the check pplacer's model
    reading does; when it rejects the model, the error says why and what to refit with. Otherwise
    the parts of the fit pplacer drops are logged once per file as a warning."""
    fit = read_fit(stats_path) if stats_path else None
    try:
        validate(pm)
    except RefpkgMalformed as e:
        if fit is None or not fit.get('model'):
            raise
        raise RefpkgMalformed(f'{e} {subject} was fitted with {fit["model"]} ({fit["program"]}); '
                              f'{_refit_hint(fit, pm)}') from None
    if fit is None or not fit.get('model') or stats_path in _warned:
        return fit
    _warned.add(stats_path)
    lost = dropped(fit, pm)
    if lost:
        logging.getLogger('timestamp').warning(
            f'{subject} was fitted with {fit["model"]} ({fit["program"]}). pplacer, and so gPlace, '
            f'place with {_describe_pplacer(pm)}; not used: {"; ".join(lost)}.')
    return fit
