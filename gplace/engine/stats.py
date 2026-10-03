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

"""The model of a statistics file (pplacer's -s, or a reference package without a phylo_model),
as pplacer v1.1.alpha19 reads it (parse_stats.ml, Gmix_model.init_of_prefs): a RAxML 7 info file
(including RAxML's re-estimated `-f e` format), a PhyML 3.0 stats file or a JSON model file. The
model is always a gamma mixture; what the file does not set comes from the command line (pplacer's
preferences: --gamma-cats, --gamma-alpha, --model-freqs).

The parsers are line-for-line ports, quirks included: lines are matched from their first
character with pplacer's regular expressions, RAxML 7.2.3-style files always get 4 gamma
categories, PhyML files only set the categories and alpha when their lines read exactly as PhyML
3.0 writes them, and later RAxML/PhyML versions mostly fail as they do in pplacer.
"""

import json
import logging
import math
import re

from gplace.exceptions import GPlaceException

KNOWN_RAXML = ['7.0.4', '7.2.3', '7.2.5', '7.2.6', '7.2.7']
KNOWN_PHYML = ['3.0', '3.0_246M']
GTR_KEYS = ('ac', 'ag', 'at', 'cg', 'ct', 'gt')


class StatsParsingError(GPlaceException):
    pass


def _ocaml_float(s):
    """OCaml float_of_string: '_' dropped, hexadecimal with 0x, else C strtod over the whole string."""
    t = s.replace('_', '')
    try:
        if not t or t != t.strip():
            raise ValueError
        if re.match(r'[+-]?0[xX]', t):
            return float.fromhex(t)
        return float(t)
    except ValueError:
        raise StatsParsingError(f"float_of_string failed on: '{s}'")


def _safe_float(s):
    """safe_float_of_string: every space, tab, newline and carriage return removed first."""
    return _ocaml_float(re.sub(r'[ \n\t\r]*', '', s))


def _ocaml_int(s):
    """int_of_string on the digits PhyML's category line captured ([0-9]*)."""
    if not s:
        raise StatsParsingError('int_of_string')
    return int(s)


def _match(rx, line):
    """Str.string_match rx line 0."""
    return re.compile(rx).match(line)


def _str_split_space(s):
    """Str.split (Str.regexp "[ ]") s: one leading and the trailing delimiter are dropped, empty
    fields between two spaces are kept."""
    start = 1 if s.startswith(' ') else 0
    out = []
    while start < len(s):
        pos = s.find(' ', start)
        if pos < 0:
            out.append(s[start:])
            break
        out.append(s[start:pos])
        start = pos + 1
    return out


def _find_beginning(rx, lines):
    """File_parsing.find_beginning: (first matching line and its match, the lines after it)."""
    for i, line in enumerate(lines):
        m = _match(rx, line)
        if m:
            return m, lines[i + 1:]
    raise LookupError


def _partition_list(rx, lines):
    """File_parsing.partition_list: sublists, each new one starting at a line matching rx (the
    first sublist starts at the first line whatever it is)."""
    if not lines:
        return [[]]
    parts = [[lines[0]]]
    for line in lines[1:]:
        if _match(rx, line):
            parts.append([line])
        else:
            parts[-1].append(line)
    return parts


def _assert_and_extract_float(beginning, s):
    for i, c in enumerate(beginning):
        if i >= len(s):
            raise StatsParsingError('index out of bounds')
        if c != s[i]:
            raise StatsParsingError("assert_and_extract_float: didn't match template!")
    return _safe_float(s[len(beginning):])


def _check_version(program, version, known):
    if version not in known:
        logging.getLogger('timestamp').warning(
            f'your stats file is from {program} {version}; {program} has been tested with the following '
            f'versions: {"; ".join(known)}')


# ------------------------------------------------------------------------- RAxML --
def _single_partition(lines):
    parts = _partition_list(r'Partition:', lines)
    if len(parts) < 2:
        raise StatsParsingError("couldn't find a partition line")
    if len(parts) > 2:
        raise StatsParsingError('too many partitions. Only one is allowed.')
    return parts[1]


def _raxml_7_2_3(lines, p):
    part = _single_partition(lines)
    try:
        m, rest = _find_beginning(r'Substitution Matrix: (.*)', part)
    except LookupError:
        raise StatsParsingError("couldn't find substitution matrix line")
    p['model_name'] = m.group(1)
    try:
        a, _ = _find_beginning(r'.*alpha\[0\]: ([^ ]*) (.*)', rest)
    except LookupError:
        raise StatsParsingError("couldn't find alpha line")
    rate_info = a.group(2)
    p['gamma_n_cat'] = 4                                 # raxml gamma always 4 categories
    p['gamma_alpha'] = _safe_float(a.group(1))
    r = _match(r'rates\[0\] ac ag at cg ct gt: (.*)', rate_info)
    if r:
        if p['model_name'] != 'GTR':
            raise StatsParsingError('have rates but model is not GTR! GTR is only allowed nucleotide model.')
        return [_safe_float(x) for x in _str_split_space(r.group(1))]
    return None                                          # pplacer's `== "GTR"` check never fires


def _raxml_re_estimated(lines, p):
    part = _single_partition(lines)
    try:
        m, rest = _find_beginning(r'DataType: (.*)', part)
    except LookupError:
        raise StatsParsingError("couldn't find data type line")
    data_type = m.group(1)
    if not rest:
        raise StatsParsingError('unexpected end of file after subs matrix line')
    sm = _match(r'Substitution Matrix: (.*)', rest[0])
    if not sm:
        raise StatsParsingError("couldn't match substitution matrix line")
    p['model_name'] = sm.group(1)
    try:
        a, final = _find_beginning(r'alpha: (.*)', rest)
    except LookupError:
        raise StatsParsingError("couldn't find alpha line")
    p['gamma_n_cat'] = 4
    p['gamma_alpha'] = _safe_float(a.group(1))
    if data_type == 'DNA':
        if not final:
            raise StatsParsingError('unexpected end of file after alpha line')
        rate_lines = final[1:7]
        templates = ['rate A <-> C:', 'rate A <-> G:', 'rate A <-> T:', 'rate C <-> G:', 'rate C <-> T:',
                     'rate G <-> T:']
        if len(rate_lines) != len(templates):
            raise StatsParsingError('Array.map2: arrays must have the same length')
        return [_assert_and_extract_float(t, s) for t, s in zip(templates, rate_lines)]
    if data_type == 'AA':
        return None
    raise StatsParsingError(f'data type unknown: {data_type}')


def _raxml(version, lines, p):
    _check_version('RAxML', version, KNOWN_RAXML)
    if version in KNOWN_RAXML and any(_match(r'.*raxmlHPC.*-f e (.*)', l) for l in lines):
        return _raxml_re_estimated(lines, p)
    return _raxml_7_2_3(lines, p)


# ------------------------------------------------------------------------- PhyML --
def _phyml(version, lines, p):
    _check_version('phyml', version, KNOWN_PHYML)
    try:
        m, after_m = _find_beginning(r'. Model of .* substitution:[ \t]+([A-Z]+)', lines)
    except LookupError:
        raise StatsParsingError("couldn't find model line!")
    p['model_name'] = m.group(1)
    try:
        g, after_g = _find_beginning(r'. Discrete gamma model:[ \t]*(.*)', after_m)
    except LookupError:
        raise StatsParsingError('Not_found')
    if g.group(1) == 'Yes':
        if len(after_g) < 2:
            raise StatsParsingError('not enough lines after gamma!')
        c = _match(r'  - Number of categories:[ \t]*([0-9]*)', after_g[0])
        if c:
            p['gamma_n_cat'] = _ocaml_int(c.group(1))
        a = _match(r'  - Gamma shape parameter:[ \t]*([0-9.]*)', after_g[1])
        if a:
            p['gamma_alpha'] = _safe_float(a.group(1))
    try:
        _, nuc = _find_beginning(r'. Nucleotides frequencies:', after_g)
    except LookupError:
        return None                                      # pplacer's `== "GTR"` check never fires
    if p['model_name'] != 'GTR':
        raise StatsParsingError('model must be GTR if we are using nucleotides')
    if not nuc:
        raise StatsParsingError('tl')
    after_freq = nuc[1:][4:]
    rate_lines = after_freq[3:][:6]
    templates = ['  A <-> C', '  A <-> G', '  A <-> T', '  C <-> G', '  C <-> T', '  G <-> T']
    if len(rate_lines) != len(templates):
        raise StatsParsingError('map2: arrays must have the same length')
    return [_assert_and_extract_float(t, s) for t, s in zip(templates, rate_lines)]


# -------------------------------------------------------------------------- JSON --
def _json(path, p):
    try:
        root = json.load(open(path))
    except (OSError, ValueError) as e:
        raise StatsParsingError(f'{path}: {e}')

    def find(obj, key, kind):
        if not isinstance(obj, dict) or key not in obj:
            raise StatsParsingError(f'missing field {key}')
        v = obj[key]
        ok = {'string': isinstance(v, str), 'int': isinstance(v, int) and not isinstance(v, bool),
              'float': isinstance(v, float)}.get(kind, True)       # Simple_json: 1 is not a float
        if not ok:
            raise StatsParsingError(f'field {key} is not a {kind}')
        return v

    subs = find(root, 'subs_model', 'string')
    p['model_name'] = subs
    ras = find(root, 'ras_model', None)
    if ras is None:
        pass
    elif ras == 'gamma':
        g = find(root, 'gamma', None)
        p['gamma_n_cat'] = find(g, 'n_cats', 'int')
        p['gamma_alpha'] = find(g, 'alpha', 'float')
    elif isinstance(ras, str):
        raise StatsParsingError(f'Unknown_ras_model("{ras}")')
    else:
        raise StatsParsingError('Wrong_field_type("ras_model")')
    datatype = find(root, 'datatype', 'string')
    if datatype == 'DNA':
        if subs != 'GTR':
            raise StatsParsingError(f'Unsupported_model("{subs}")')
        rates = find(root, 'subs_rates', None)
        return [find(rates, k, 'float') for k in GTR_KEYS]
    if datatype == 'AA':
        if subs != 'WAG':
            raise StatsParsingError(f"We don't currently support {subs}. Unsupported_model(\"{subs}\")")
        return None
    raise StatsParsingError(f'Unknown_datatype("{datatype}")')


# -------------------------------------------------------------------------- main --
def _lines(path):
    """The lines as OCaml's input_line returns them (split at '\\n' only)."""
    with open(path, newline='') as fh:
        text = fh.read()
    lines = text.split('\n')
    if lines and lines[-1] == '':
        lines.pop()
    return lines


def parse_stats(path, gamma_n_cat=1, gamma_alpha=1.0, model_name='LG'):
    """Parse_stats.parse_stats: (model_name, gamma_n_cat, gamma_alpha, GTR rates or None), each
    setting not given by the file keeping the value passed in."""
    p = {'model_name': model_name, 'gamma_n_cat': gamma_n_cat, 'gamma_alpha': gamma_alpha}
    try:
        if path.endswith('.json'):
            trans = _json(path, p)
        else:
            lines = _lines(path)
            v = next((m.group(1) for m in (_match(r'.* RAxML version ([^ ]+)', l) for l in lines) if m), None)
            if v is not None:
                trans = _raxml(v, lines, p)
            else:
                v = next((m.group(1) for m in (_match(r'[ \t]*---  PhyML v([^ ]+)', l) for l in lines) if m),
                         None)
                if v is None:
                    raise StatsParsingError('is this a RAxML v7 info or PHYML v3 statistics file? The header '
                                            "didn't match.")
                trans = _phyml(v, lines, p)
    except StatsParsingError as e:
        raise StatsParsingError(f'Problem parsing info or stats file {path}: {e}')
    except OSError as e:
        raise StatsParsingError(f'Problem parsing info or stats file {path}: {e}')
    return p['model_name'], p['gamma_n_cat'], p['gamma_alpha'], trans


def phylo_model_of_stats(path, empirical_freqs=True, gamma_n_cat=1, gamma_alpha=1.0, model_name='LG'):
    """The phylo_model dict (models.model_from_refpkg) pplacer builds from a stats file
    (Gmix_model.init_of_prefs): always a gamma mixture; empirical frequencies unless --model-freqs
    (GTR always uses the alignment's frequencies)."""
    name, n_cat, alpha, trans = parse_stats(path, gamma_n_cat, gamma_alpha, model_name)
    pm = {'subs_model': name, 'ras_model': 'gamma', 'empirical_frequencies': bool(empirical_freqs),
          'gamma': {'n_cats': n_cat, 'alpha': alpha}, 'source': path}
    if trans is not None:
        if len(trans) != 6:
            raise StatsParsingError(f'Problem parsing info or stats file {path}: {len(trans)} GTR rates, '
                                    f'expected 6 (ac ag at cg ct gt).')
        pm['subs_rates'] = dict(zip(GTR_KEYS, trans))
    elif name == 'GTR':
        raise StatsParsingError('GTR specified but no substitution rates given.')
    if not (isinstance(n_cat, int) and n_cat > 0) or math.isnan(alpha):
        raise StatsParsingError(f'Invalid gamma model from {path}: {n_cat} categories, alpha {alpha}.')
    return pm
