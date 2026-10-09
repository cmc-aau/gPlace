#!/usr/bin/env python3
"""Extract IQ-TREE's built-in protein models from its source into gplace/engine/data/iqtree_protein_models.json,
verbatim (numbers as written): the reversible matrices of model/modelprotein.cpp and model/modelmixture.cpp
(exchangeabilities and frequencies, PAML layout, amino acids ARNDCQEGHILKMFPSTWYV), the frequency profiles
(C10..C60: "frequency NAME = ...") and the mixture definitions ("model NAME = MIX{...}", FMIX{...}, ...).

    python scripts/make_iqtree_models.py <iqtree3 source dir> <version>

Non-reversible matrices (NQ.*: a full 20 x 20 matrix) are left out."""
import json, os, re, sys

src, version = sys.argv[1], sys.argv[2]
models, skipped, profiles, mixtures = {}, [], {}, {}
NUM = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$')
for source in ('modelprotein.cpp', 'modelmixture.cpp'):
    text = re.sub(r'\[.*?\]', ' ', open(os.path.join(src, 'model', source)).read(), flags=re.DOTALL)
    for m in re.finditer(r'^(model|frequency) ([A-Za-z0-9._]+)\s*=\s*(.*?);', text, re.MULTILINE | re.DOTALL):
        kind, name, body = m.group(1), m.group(2), m.group(3)
        words = body.split()
        if kind == 'frequency':
            if len(words) == 20 and all(NUM.match(w) for w in words):
                profiles[name] = words
            continue
        if words and all(NUM.match(w) for w in words):
            if len(words) == 210:
                models[name] = {'exchangeabilities': words[:190], 'frequencies': words[190:]}
            else:
                skipped.append(f'{name} ({len(words)} numbers)')
        else:
            mixtures[name] = ' '.join(words)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gplace', 'engine', 'data', 'iqtree_protein_models.json')
with open(out, 'w') as fh:
    json.dump({'source': f'IQ-TREE {version}, model/modelprotein.cpp and model/modelmixture.cpp (GPL-2.0-or-later); '
                         'PAML layout: 190 exchangeabilities (lower triangle, rows R..V) then 20 frequencies, order '
                         'ARNDCQEGHILKMFPSTWYV',
               'models': models, 'frequencies': profiles, 'mixtures': mixtures}, fh, indent=1)
print(f'{len(models)} reversible matrices, {len(profiles)} frequency profiles, {len(mixtures)} mixture definitions '
      f'written to {out}')
print('left out:', ', '.join(skipped) or 'none')
