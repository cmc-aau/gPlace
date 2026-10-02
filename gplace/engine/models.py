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

"""The evolutionary model of a reference package, chosen as pplacer v1.1.alpha19 chooses it.

For a refpkg, pplacer reads the model from its phylo_model JSON (refpkg.ml; the -m option is
ignored):
  * ras_model "Price-CAT" -> Gcat_model (FastTree CAT: one rate per site, from "Price-CAT"),
    ras_model "gamma"     -> Gmix_model (every site under all discrete-gamma rates, averaged);
  * subs_model WAG, LG or JTT (protein; Prot_models) or GTR (DNA; "subs_rates" ac..gt);
  * "empirical_frequencies": frequencies counted from the reference alignment
    (Alignment.emper_freq) instead of the model's; GTR always uses empirical frequencies
    (Gstar_support.seqtype_and_trans_statd_of_info).
"""

import json
import os

import numpy as np

from gplace.engine import core
from gplace.exceptions import RefpkgMalformed

PROTEIN_MODELS = ('WAG', 'LG', 'JTT')
GTR_KEYS = ('ac', 'ag', 'at', 'cg', 'ct', 'gt')


class Model(object):
    """subs_model, ras ('cat' or 'gamma'), alphabet, exchangeabilities R, frequencies statd, rates,
    site categories (CAT), gamma (alpha, n_cats) and empirical (bool)."""

    def __init__(self, subs_model, ras, alphabet, R, statd, rates, site_categories, empirical, gamma=None):
        self.subs_model, self.ras, self.alphabet = subs_model, ras, alphabet
        self.R, self.statd, self.rates = R, statd, rates
        self.site_categories, self.empirical, self.gamma = site_categories, empirical, gamma

    @property
    def ns(self):
        return self.alphabet.ns

    @property
    def n_rates(self):
        return len(self.rates)

    @property
    def width(self):
        """Doubles per site in a likelihood vector: states, times the rates for gamma."""
        return self.ns * (self.n_rates if self.ras == 'gamma' else 1)

    def describe(self):
        ras = (f'gamma ({self.n_rates} rates, alpha {self.gamma[0]:g})' if self.ras == 'gamma'
               else f'CAT ({self.n_rates} rates)')
        freqs = 'empirical' if self.empirical else 'model'
        return f'{self.subs_model} + {ras}, {freqs} frequencies'


def l1_norm(v):
    """Linear_utils.lp_norm 1.: x = 0; x = x +. v.{i} ** 1. (sequential); x ** 1. (= x)."""
    x = 0.0
    for a in v:
        x = x + float(a)
    return x


def emper_freq(align, alphabet, chunk=1 << 20):
    """Alignment.emper_freq over the (uppercased, unmasked) reference alignment: every character
    except '-' and '?' adds its likelihood vector scaled by 1 /. l1_norm, sequence by sequence in
    alignment order; the total is then scaled the same way. The additions are sequential per state
    (np.add.accumulate), as Gsl_vector.add does them: ambiguity codes add fractions (X: 1/20 to
    every amino acid), so the order of the additions shows in the result."""
    ns = alphabet.ns
    lut = np.full((256, ns), np.nan)
    for ch, v in alphabet.char_map.items():
        if ch in '-?':
            continue
        w = np.array(v, dtype=float)
        lut[ord(ch)] = w * (1.0 / l1_norm(w))
    total = np.zeros(ns)
    for name, seq in align:
        b = np.frombuffer(seq.encode(), np.uint8)
        b = b[(b != ord('-')) & (b != ord('?'))]
        bad = np.isnan(lut[b, 0])
        if bad.any():
            raise RefpkgMalformed(f"'{chr(b[bad][0])}' not a known base in {name}!")
        for i in range(0, len(b), chunk):
            total = np.add.accumulate(np.vstack([total[None, :], lut[b[i:i + chunk]]]), axis=0)[-1]
    return list(total * (1.0 / l1_norm(total)))


def gtr_exchangeabilities(subs_rates):
    """Nuc_models.b_of_trans_vector: symmetric 4x4 (order A C G T) from ac ag at cg ct gt."""
    try:
        v = [float(subs_rates[k]) for k in GTR_KEYS]
    except (KeyError, TypeError, ValueError):
        raise RefpkgMalformed('GTR specified but the substitution rates (subs_rates: ac ag at cg ct gt) '
                              'are missing or invalid.')
    R = [[0.0] * 4 for _ in range(4)]
    for k, (i, j) in enumerate([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]):
        R[i][j] = R[j][i] = v[k]
    return R


def discrete_gamma(n_cats, alpha, lib=None):
    """Gamma.discrete_gamma n_cats alpha: Yang (1994) mean rates of n_cats equally likely
    categories, scaled to mean 1, with GSL 1.16's gamma functions (in the engine)."""
    if int(n_cats) < 1 or not alpha > 0:
        raise RefpkgMalformed(f'Invalid gamma model: n_cats {n_cats}, alpha {alpha}.')
    L = lib or core.lib()
    out = np.empty(int(n_cats))
    if L.xp_discrete_gamma(int(n_cats), float(alpha), core.P(out)) != 0:
        raise RefpkgMalformed(f'Invalid gamma model: n_cats {n_cats}, alpha {alpha}.')
    return [float(v) for v in out]


def read_phylo_model(refpkg_dir):
    files = json.load(open(os.path.join(refpkg_dir, 'CONTENTS.json'))).get('files', {})
    if 'phylo_model' not in files:
        raise RefpkgMalformed(f'Reference package {refpkg_dir} has no "phylo_model" file.')
    return json.load(open(os.path.join(refpkg_dir, files['phylo_model'])))


def alphabet_of(pm):
    """The alphabet a phylo_model JSON needs (GTR: DNA; WAG/LG/JTT: protein)."""
    subs = pm.get('subs_model')
    if subs == 'GTR':
        return core.DNA
    if subs in PROTEIN_MODELS:
        return core.PROTEIN
    raise RefpkgMalformed(f'Unsupported substitution model "{subs}": pplacer v1.1.alpha19 knows WAG, '
                          f'LG and JTT (protein) and GTR (DNA).')


def width_of(pm):
    """Doubles per site of the likelihood vectors for a phylo_model JSON (memory planning)."""
    ns = alphabet_of(pm).ns
    if pm.get('ras_model') == 'gamma':
        return ns * int(pm.get('gamma', {}).get('n_cats', 4))
    return ns


def model_from_refpkg(pm, ref_align, lib=None):
    """The Model of a refpkg from its phylo_model JSON (pm) and its reference alignment
    [(name, sequence)] (uppercased; the model is built before pplacer's pre-masking)."""
    alphabet = alphabet_of(pm)
    subs, ras = pm.get('subs_model'), pm.get('ras_model')
    if subs == 'GTR':
        R = gtr_exchangeabilities(pm.get('subs_rates'))
        statd, empirical = emper_freq(ref_align, alphabet), True
    else:
        R, statd = core.load_protein_model(subs)
        empirical = bool(pm.get('empirical_frequencies'))
        if empirical:
            statd = emper_freq(ref_align, alphabet)
    if ras == 'Price-CAT':
        cat = pm.get('Price-CAT') or {}
        try:
            rates = [float(v) for v in cat['Rates']]
            site_categories = np.array([int(v) - 1 for v in cat['SiteCategories']], dtype=np.int32)
        except (KeyError, TypeError, ValueError):
            raise RefpkgMalformed('Price-CAT model without valid "Rates" and "SiteCategories".')
        return Model(subs, 'cat', alphabet, R, statd, rates, site_categories, empirical)
    if ras == 'gamma':
        g = pm.get('gamma') or {}
        try:
            n_cats, alpha = int(g['n_cats']), float(g['alpha'])
        except (KeyError, TypeError, ValueError):
            raise RefpkgMalformed('gamma model without valid "n_cats" and "alpha".')
        rates = discrete_gamma(n_cats, alpha, lib)
        return Model(subs, 'gamma', alphabet, R, statd, rates, None, empirical, (alpha, n_cats))
    raise RefpkgMalformed(f'Unsupported rate model "{ras}": pplacer knows "Price-CAT" and "gamma".')
