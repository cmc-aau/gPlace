#!/usr/bin/env python3
"""Extract RAxML-NG's built-in protein models from its source (coraxlib) into
gplace/engine/data/raxmlng_protein_models.json, verbatim (numbers as written): every matrix of
src/corax/model/modutil_aa.c with its tables in src/corax/model/aa_rates.c (190 exchangeabilities, upper
triangle row by row: A-R, A-N, ..., Y-V; 20 frequencies; amino acids ARNDCQEGHILKMFPSTWYV) and the mixtures
LG4M (gamma rates) and LG4X (free rates).

    python scripts/make_raxmlng_models.py <raxml-ng source dir> <version>

PROTGTR has no tables (its rates are always the fit's)."""
import json, os, re, sys

src, version = sys.argv[1], sys.argv[2]
model_dir = os.path.join(src, 'libs', 'coraxlib', 'src', 'corax', 'model')
strip = lambda t: re.sub(r'//[^\n]*', '', re.sub(r'/\*.*?\*/', '', t, flags=re.DOTALL))
rates_c = strip(open(os.path.join(model_dir, 'aa_rates.c')).read())
modutil = strip(open(os.path.join(model_dir, 'modutil_aa.c')).read())
NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def table(symbol):
    """The numbers of a table (symbol, or symbol[k] of a [4][n] table)."""
    m = re.fullmatch(r'(\w+)(?:\[(\d)\])?', symbol)
    name, k = m.group(1), m.group(2)
    t = re.search(r'const double ' + name + r'\s*(\[\d+\])+\s*=\s*\{(.*?)\};', rates_c, re.DOTALL)
    if t is None:
        raise SystemExit(f'table {name} not found')
    body = t.group(2)
    if k is not None:                                   # [4][n]: the k-th brace block
        blocks = re.findall(r'\{([^{}]*)\}', body)
        body = blocks[int(k)]
    return re.findall(NUM, body)


models = {}
for m in re.finditer(r'const corax_subst_model_t M_\w+\s*=\s*\{\s*"([^"]+)",\s*20,\s*([\w\[\]]+),\s*([\w\[\]]+)',
                     modutil):
    name, rates, freqs = m.groups()
    if rates == 'NULL':
        continue                                        # PROTGTR
    r, f = table(rates), table(freqs)
    if len(r) != 190 or len(f) != 20:
        raise SystemExit(f'{name}: {len(r)} rates, {len(f)} frequencies')
    models[name] = {'exchangeabilities': r, 'frequencies': f}
mixtures = {}
for m in re.finditer(r'const corax_subst_model_t \*(\w+)_matrices\[\]\s*=\s*\{([^}]*)\};', modutil):
    comps = [c.strip().lstrip('&').replace('M_', '', 1) for c in m.group(2).split(',') if c.strip()]
    mix = re.search(r'const corax_mixture_model_t M_(\w+)\s*=\s*\{\s*"([^"]+)",\s*(\d+),\s*\(corax_subst_model_t \*\*\)'
                    + m.group(1) + r'_matrices,\s*NULL,\s*NULL,\s*CORAX_UTIL_MIXTYPE_(\w+)\s*\}', modutil)
    names = [next(n for n, v in models.items() if n.replace('-', '').upper() == c.upper()) for c in comps]
    mixtures[mix.group(2)] = {'components': names, 'rates': {'GAMMA': 'gamma', 'FREE': 'free'}[mix.group(4)]}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gplace', 'engine', 'data',
                   'raxmlng_protein_models.json')
with open(out, 'w') as fh:
    json.dump({'source': f'RAxML-NG {version}, coraxlib src/corax/model/aa_rates.c and modutil_aa.c; 190 '
                         'exchangeabilities (upper triangle, row by row: A-R, A-N, ..., Y-V) then 20 frequencies, '
                         'order ARNDCQEGHILKMFPSTWYV',
               'models': models, 'mixtures': mixtures}, fh, indent=1)
print(f'{len(models)} matrices, {len(mixtures)} mixtures ({", ".join(mixtures)}) written to {out}')
