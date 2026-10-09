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

"""Yang's discrete gamma rates (1994) as the fitting programs compute them (place --full_model).

IQ-TREE (model/rategamma.cpp) and RAxML-NG (coraxlib, model/gamma.c) both use the functions of Yang's PAML
(tools.c): the percentage points of the chi-square distribution and the incomplete gamma ratio, to their own
accuracy (1e-8). Their category rates differ from GSL's (pplacer's, models.discrete_gamma) by 1e-8 to 1e-7
relative, and much more at small alpha, where GSL's lowest rates underflow to 0. The functions here are
computed step by step as there (bit-identical results); the two programs then differ only in how they
rescale the rates (raxmlng_rates, iqtree_rates).
"""

import math

from gplace.exceptions import RefpkgMalformed


def ln_gamma(alpha):
    """ln(Gamma(alpha)), Pike & Hill (1966)."""
    x, f = alpha, 0.0
    if x < 7.0:
        f, z = 1.0, alpha - 1.0
        while True:
            z = z + 1.0
            if not z < 7.0:
                break
            f *= z
        x = z
        f = -math.log(f)
    z = 1 / (x * x)
    return (f + (x - 0.5) * math.log(x) - x + .918938533204673
            + (((-.000595238095238 * z + .000793650793651) * z - .002777777777778) * z + .083333333333333) / x)


def incomplete_gamma(x, alpha, ln_gamma_alpha):
    """I(x, alpha), Bhattacharjee (1970) AS 32; -1 on error."""
    p, g = alpha, ln_gamma_alpha
    accurate, overflow = 1e-8, 1e30
    if x == 0:
        return 0.0
    if x < 0 or p <= 0:
        return -1.0
    factor = math.exp(p * math.log(x) - x - g)
    if not (x > 1 and x >= p):                          # series expansion
        gin, term, rn = 1.0, 1.0, p
        while True:
            rn += 1
            term *= x / rn
            gin += term
            if not term > accurate:
                break
        gin *= factor / p
        return gin
    a = 1 - p                                           # continued fraction
    b = a + x + 1
    term = 0.0
    pn = [1.0, x, x + 1, x * b, 0.0, 0.0]
    gin = pn[2] / pn[3]
    while True:
        a += 1
        b += 2
        term += 1
        an = a * term
        for i in range(2):
            pn[i + 4] = b * pn[i + 2] - an * pn[i]
        if pn[5] != 0:
            rn = pn[4] / pn[5]
            dif = abs(gin - rn)
            if dif <= accurate and dif <= accurate * rn:
                return 1 - factor * gin
            gin = rn
        for i in range(4):
            pn[i] = pn[i + 2]
        if not abs(pn[4]) < overflow:
            for i in range(4):
                pn[i] /= overflow


def point_normal(prob):
    """z with P(X < z) = prob, X ~ N(0, 1), Odeh & Evans (1974) AS 70."""
    a0, a1, a2, a3 = -.322232431088, -1.0, -.342242088547, -.0204231210245
    a4, b0, b1 = -.453642210148e-4, .0993484626060, .588581570495
    b2, b3, b4 = .531103462366, .103537752850, .0038560700634
    p = prob
    p1 = p if p < 0.5 else 1 - p
    if p1 < 1e-20:
        return -9999.0
    y = math.sqrt(math.log(1 / (p1 * p1)))
    z = y + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0) / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0)
    return -z if p < 0.5 else z


def point_chi2(prob, v):
    """z with P(X < z) = prob, X ~ Chi2(v), Best & Roberts (1975) AS 91; -1 on error."""
    e, aa, p = .5e-6, .6931471805, prob
    if p < .000002 or p > .999998 or v <= 0:
        return -1.0
    g = ln_gamma(v / 2)
    xx = v / 2
    c = xx - 1
    if v < -1.24 * math.log(p):
        ch = math.pow(p * xx * math.exp(g + xx * aa), 1 / xx)
        if ch - e < 0:
            return ch
    elif v <= .32:
        ch = 0.4
        a = math.log(1 - p)
        while True:
            q = ch
            p1 = 1 + ch * (4.67 + ch)
            p2 = ch * (6.73 + ch * (6.66 + ch))
            t = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch * (13.32 + 3 * ch)) / p2
            ch -= (1 - math.exp(a + g + .5 * ch + c * aa) * p2 / p1) / t
            if abs(q / ch - 1) - .01 <= 0:
                break
    else:
        x = point_normal(p)
        p1 = 0.222222 / v
        ch = v * math.pow(x * math.sqrt(p1) + 1 - p1, 3.0)
        if ch > 2.2 * v + 6:
            ch = -2 * (math.log(1 - p) - c * math.log(.5 * ch) + g)
    while True:
        q = ch
        p1 = .5 * ch
        t = incomplete_gamma(p1, xx, g)
        if t < 0.0:
            return -1.0
        p2 = p - t
        t = p2 * math.exp(xx * aa + g + p1 - c * math.log(ch))
        b = t / ch
        a = 0.5 * t - b * c
        s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420
        s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520
        s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520
        s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040
        s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520
        s6 = (120 + c * (346 + 127 * c)) / 5040
        ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))
        if not abs(q / ch - 1) > e:
            return ch


def mean_rates(n, alpha):
    """Yang's mean rates of n equally likely categories (Eq. 9-10), times n as both programs compute them."""
    lnga1 = ln_gamma(alpha + 1)
    probs = [point_chi2((i + 1.0) / n, 2.0 * alpha) / (2.0 * alpha) for i in range(n - 1)]
    probs = [incomplete_gamma(x * alpha, alpha + 1, lnga1) for x in probs]
    return ([probs[0] * n] + [(probs[i] - probs[i - 1]) * n for i in range(1, n - 1)]
            + [(1 - probs[n - 2]) * n])


def _check(n, alpha, program):
    if not alpha > 0 or n < 1:
        raise RefpkgMalformed(f'--full_model: invalid gamma model ({n} categories, alpha {alpha}) in the {program} fit.')


def raxmlng_rates(n, alpha, median=False):
    """RAxML-NG's rates (corax_compute_gamma_cats): mean, or median rescaled to mean 1."""
    if alpha < 0.02 or n < 1:
        raise RefpkgMalformed(f'--full_model: invalid gamma model ({n} categories, alpha {alpha}): RAxML-NG needs '
                              f'alpha >= 0.02.')
    if n == 1:
        return [1.0]
    factor = alpha / alpha * n
    if median:
        middle = 1.0 / (2.0 * n)
        rates = [point_chi2(float(i * 2 + 1) * middle, 2.0 * alpha) / (2.0 * alpha) for i in range(n)]
        t = 0.0
        for r in rates:
            t += r
        return [r * (factor / t) for r in rates]
    return mean_rates(n, alpha)


def iqtree_rates(n, alpha, median=False, start=1.0):
    """IQ-TREE's rates (RateGamma::computeRates): mean, or median rescaled to mean 1; then rescaled so that
    their sum is the sum of the rates they replace, which start at `start` (1, or 1 / (1 - p) with +I:
    RateGammaInvar), so that with +I they average 1 / (1 - p)."""
    _check(n, alpha, 'IQ-TREE')
    if n == 1:
        return [1.0]
    cur = 0.0
    for _ in range(n):
        cur += start
    if median:
        rates = []
        for cat in range(n):
            r = point_chi2((2.0 * cat + 1) / (2.0 * n), 2.0 * alpha) / (2.0 * alpha)
            rates.append(-r if r < 0.0 else r)
        total = 0.0
        for r in rates:
            total += r
        rates = [r * n / total for r in rates]
    else:
        rates = mean_rates(n, alpha)
    new = 0.0
    for r in rates:
        new += r
    if new != cur:
        rates = [r * (cur / new) for r in rates]
    return rates
