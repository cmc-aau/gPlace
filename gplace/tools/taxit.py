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

"""``taxit create``, ``taxit check`` and ``taxit info`` (taxtastic v1.0.1, GPL-3.0,
https://github.com/fhcrc/taxtastic): reference packages for pplacer.

Ported from taxtastic's ``refpkg.py``, ``utils.py`` and ``subcommands/{create,check,info}.py`` so
that a refpkg comes out as taxit writes it: the same files, the same ``CONTENTS.json`` (key order,
log, md5 sums, ``indent=4``) and the same ``phylo_model.json`` parsed from the RAxML, RAxML-NG,
IQ-TREE, FastTree or PhyML log. As with taxit, the creation date and the file name of
``phylo_model`` (a temporary file name) differ from run to run. With a taxonomy, seq_info and tree,
taxit reroots the tree with pplacer's ``rppr reroot`` when ``rppr`` is on PATH, and otherwise warns
and skips it; so does this.
"""

import copy
import csv
import errno
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

from gplace.engine.core import read_fasta
from gplace.engine.newick import parse_newick
from gplace.fileio import read_text

FORMAT_VERSION = '1.1'
MANIFEST = 'CONTENTS.json'


def _warn(msg):
    print(msg, file=sys.stderr)


# --------------------------------------------------------- stats-file parsers (utils.py) ---

def try_set_fields(d, regex, text, hook=lambda x: x):
    v = re.search(regex, text, re.MULTILINE)
    if v:
        d.update(dict([(key, hook(val)) for key, val in v.groupdict().items()]))
    return d


class InvalidLogError(ValueError):
    pass


# https://github.com/amkozlov/raxml-ng/wiki/Input-data#single-mode
DNA = ['JC', 'K80', 'F81', 'HKY', 'TN93ef', 'TN93', 'K81', 'K81uf', 'TPM2',
       'TPM2uf', 'TPM3', 'TPM3uf', 'TIM1', 'TIM1uf', 'TIM2', 'TIM2uf',
       'TIM3', 'TIM3uf', 'TVMef', 'TVM', 'SYM', 'GTR']
PROT = ['Blosum62', 'cpREV', 'Dayhoff', 'DCMut', 'DEN', 'FLU', 'HIVb',
        'HIVw', 'JTT', 'JTT-DCMut', 'LG', 'mtART', 'mtMAM', 'mtREV',
        'mtZOA', 'PMB', 'rtREV', 'stmtREV', 'VT', 'WAG', 'LG4M', 'LG4X',
        'PROTGTR']


def parse_raxmlng(handle):
    """Parse RAxML-NG's log (taxtastic: protein models get datatype 'RNA', kept as is)."""
    s = handle.read()
    result = {}
    try_set_fields(result, r'(?P<program>RAxML-NG v. [\d\.]+)', s)
    try_set_fields(result, r'Model: (?P<subs_model>\w+)', s)
    try_set_fields(result, r'raxml-ng.*?--data-type (?P<datatype>\w+)', s)
    if 'datatype' not in result:
        if result['subs_model'] in DNA:
            result['datatype'] = 'DNA'
        elif result['subs_model'] in PROT:
            result['datatype'] = 'RNA'  # AA specified in --data-type
        else:
            result['datatype'] = ''  # There's more but stop here for now
    result['empirical_frequencies'] = re.search(r'Base frequencies \(ML\)', s) is not None
    rates = {}
    try_set_fields(rates, r'Substitution rates \(ML\): (?P<ac>\d+\.\d+) (?P<ag>\d+\.\d+) (?P<at>\d+\.\d+) '
                          r'(?P<cg>\d+\.\d+) (?P<ct>\d+\.\d+) (?P<gt>\d+\.\d+)', s, hook=float)
    if len(rates) > 0:
        result['subs_rates'] = rates
    gamma = {}
    try_set_fields(gamma, r'Rate heterogeneity: GAMMA \((?P<n_cats>\d+) cats, mean\),  alpha: (?P<alpha>\d+\.\d+)',
                   s, hook=float)
    if gamma:
        gamma['n_cats'] = int(gamma['n_cats'])
        result['gamma'] = gamma
        result['ras_model'] = 'gamma'
    return result


def parse_iqtree(handle):
    """Parse IQ-TREE's .iqtree summary."""
    s = handle.read()
    result = {}
    try_set_fields(result, r'(?P<program>IQ-TREE [\d\.]+)', s)
    try_set_fields(result, r'Model of substitution: (?P<model>[\w\+]+)\n', s)
    if 'model' in result:
        models = set(result['model'].split('+'))
        subs_model = list(set(DNA).intersection(models))
        if len(subs_model) > 0:
            result['subs_model'] = subs_model[0]
            result['datatype'] = 'DNA'
        else:
            subs_model = list(set(PROT).intersection(models))
            if len(subs_model) > 0:
                result['subs_model'] = subs_model[0]
                result['datatype'] = 'AA'
    result['empirical_frequencies'] = re.search(
        r'State frequencies: \(empirical counts from alignment\)\n', s) is not None
    re_rate = re.compile(r'(?P<first>[A-Z])\-(?P<second>[A-Z]): (?P<rate>\d+\.\d+)\n')
    rates = {f'{f.lower()}{t.lower()}': float(r) for f, t, r in re_rate.findall(s)}
    if len(rates) > 0:
        result['subs_rates'] = rates
    frequencies = {b: float(f) for (b, f) in re.findall(r'pi\((?P<base>[ACTG])\) = (\d+\.\d+)\n', s)}
    if len(frequencies) > 0:
        result['frequencies'] = frequencies
    gamma = {}
    try_set_fields(gamma, r'Model of rate heterogeneity: Gamma with (?P<n_cats>\d+) categories\n'
                          r'Gamma shape alpha: (?P<alpha>\d+\.\d+)\n', s, hook=float)
    if gamma:
        gamma['n_cats'] = int(gamma['n_cats'])
        result['gamma'] = gamma
        result['ras_model'] = 'gamma'
    return result


def parse_raxml(handle):
    """Parse RAxML's info file."""
    s = ''.join(handle.readlines())
    result = {}
    try_set_fields(result, r'(?P<program>RAxML version [0-9.]+)', s)
    try_set_fields(result, r'DataType: (?P<datatype>DNA|RNA|AA)', s)
    result['empirical_frequencies'] = (
        result['datatype'] != 'AA' or
        re.search('empirical base frequencies', s, re.IGNORECASE) is not None)
    try_set_fields(result, r'Substitution Matrix: (?P<subs_model>\w+)', s)
    rates = {}
    if result['datatype'] != 'AA':
        try_set_fields(rates, (r"rates\[0\] ac ag at cg ct gt: "
                               r"(?P<ac>[0-9.]+) (?P<ag>[0-9.]+) (?P<at>[0-9.]+) "
                               r"(?P<cg>[0-9.]+) (?P<ct>[0-9.]+) (?P<gt>[0-9.]+)"), s, hook=float)
        try_set_fields(rates, r'rate A <-> C: (?P<ac>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate A <-> G: (?P<ag>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate A <-> T: (?P<at>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate C <-> G: (?P<cg>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate C <-> T: (?P<ct>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate G <-> T: (?P<gt>[0-9.]+)', s, hook=float)
        if len(rates) > 0:
            result['subs_rates'] = rates
    result['gamma'] = {'n_cats': 4}
    try_set_fields(result['gamma'], r"alpha[\[\]0-9]*: (?P<alpha>[0-9.]+)", s, hook=float)
    result['ras_model'] = 'gamma'
    return result


JTT_MODEL = 'ML Model: Jones-Taylor-Thorton, CAT approximation with 20 rate categories'
WAG_MODEL = 'ML Model: Whelan-And-Goldman, CAT approximation with 20 rate categories'
LG_MODEL = 'ML Model: Le-Gascuel 2008, CAT approximation with 20 rate categories'


def parse_fasttree(fobj):
    """Parse FastTree's log."""
    data = {
        'empirical_frequencies': True,
        'datatype': 'DNA',
        'subs_model': 'GTR',
        'ras_model': 'Price-CAT',
        'Price-CAT': {},
    }
    for line in fobj:
        if not line.strip():
            continue
        splut = line.split()
        if splut[0] == 'FastTree':
            data['program'] = line.strip()
        elif splut[0] == 'Rates':
            data['Price-CAT']['Rates'] = list(map(float, splut[1:]))
        elif splut[0] == 'SiteCategories':
            data['Price-CAT']['SiteCategories'] = list(map(int, splut[1:]))
        elif splut[0] == 'NCategories':
            data['Price-CAT']['n_cats'] = int(splut[1])
        elif splut[0] == 'GTRRates':
            data['subs_rates'] = dict(list(zip(['ac', 'ag', 'at', 'cg', 'ct', 'gt'], list(map(float, splut[1:])))))
        elif line.strip() == JTT_MODEL:
            data['subs_model'] = 'JTT'
            data['datatype'] = 'AA'
            data['empirical_frequencies'] = False
        elif line.strip() == WAG_MODEL:
            data['subs_model'] = 'WAG'
            data['datatype'] = 'AA'
            data['empirical_frequencies'] = False
        elif line.strip() == LG_MODEL:
            data['subs_model'] = 'LG'
            data['datatype'] = 'AA'
            data['empirical_frequencies'] = False
    if data['subs_model'] == 'GTR' and 'subs_rates' not in data:
        raise InvalidLogError("GTR model, but no substitution rates found!")
    return data


def parse_phyml(fobj, frequency_type=None):
    """Parse PhyML's stats file."""
    if frequency_type not in (None, 'empirical', 'model'):
        raise ValueError("Unknown frequency_type: {0}".format(frequency_type))
    s = ''.join(fobj)
    result = {'gamma': {}}
    try_set_fields(result, r'---\s*(?P<program>PhyML.*?)\s*---', s)
    try_set_fields(result['gamma'], r'Number of categories:\s+(?P<n_cats>\d+)', s, hook=int)
    try_set_fields(result['gamma'], r'Gamma shape parameter:\s+(?P<alpha>\d+\.\d+)', s, hook=float)
    result['ras_model'] = 'gamma'
    if 'nucleotides' in s:
        result['datatype'] = 'DNA'
        try_set_fields(result, r'Model of nucleotides substitution:\s+(?P<subs_model>\w+)', s)
        rates = {}
        try_set_fields(rates, r'A <-> C\s+(?P<ac>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'A <-> G\s+(?P<ag>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'A <-> T\s+(?P<at>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'C <-> G\s+(?P<cg>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'C <-> T\s+(?P<ct>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'G <-> T\s+(?P<gt>\d+\.\d+)', s, hook=float)
        if rates:
            result['subs_rates'] = rates
        # PhyML doesn't record whether empirical base frequencies were used, or ML estimates were
        # made; taxtastic sets empirical.
        result['empirical_frequencies'] = True
    elif 'amino acids' in s:
        result['datatype'] = 'AA'
        if not frequency_type:
            raise ValueError("frequency type required for PhyML AA models.")
        result['empirical_frequencies'] = frequency_type == 'empirical'
        try_set_fields(result, r'Model of amino acids substitution:\s+(?P<subs_model>\w+)', s)
    else:
        raise ValueError('Could not determine if alignment is AA or DNA')
    return result


def parse_stockholm(fobj):
    """Names of a Stockholm alignment (ValueError without the // terminator)."""
    names = {}
    found_eof = False
    for line in fobj:
        line = line.strip()
        if line == '//':
            found_eof = True
        elif line.startswith('#') or not line.strip():
            continue
        else:
            name, __ = line.split(None, 1)
            names[name] = None
    if not found_eof:
        raise ValueError('Invalid Stockholm format: no file terminator')
    return list(names.keys())


def guess_stats_type(stats_file):
    with open(stats_file) as fobj:
        for line in fobj:
            if line.startswith('FastTree'):
                return 'FastTree'
            elif line.startswith('This is RAxML') or line.startswith('You are using RAxML'):
                return 'RAxML'
            elif line.startswith('RAxML-NG'):
                return 'RAxML-NG'
            elif 'PhyML' in line:
                return 'PhyML'
            elif line.startswith("IQ-TREE"):
                return 'IQTREE'
    raise ValueError("couldn't guess log type for %r" % (stats_file,))


def phylo_model_of(stats_type, stats_file, frequency_type=None):
    """Refpkg.update_phylo_model's parsing: the phylo_model dictionary."""
    if frequency_type not in (None, 'model', 'empirical'):
        raise ValueError('Unknown frequency type: "{0}"'.format(frequency_type))
    if frequency_type and stats_type not in (None, 'PhyML'):
        raise ValueError('Frequency type should only be specified for PhyML alignments.')
    if stats_type is None:
        stats_type = guess_stats_type(stats_file)
    parser = {'RAxML': parse_raxml, 'RAxML-NG': parse_raxmlng, 'FastTree': parse_fasttree,
              'IQTREE': parse_iqtree,
              'PhyML': lambda h: parse_phyml(h, frequency_type=frequency_type)}.get(stats_type)
    if parser is None:
        raise ValueError('invalid log type: %r' % (stats_type,))
    with open(stats_file) as h:
        return parser(h)


# ---------------------------------------------------------------------- refpkg (refpkg.py) ---

def md5file(path):
    with open(path, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def manifest_template():
    return {'metadata': {'create_date': time.strftime('%Y-%m-%d %H:%M:%S'), 'format_version': FORMAT_VERSION},
            'files': {}, 'md5': {}, 'log': [], 'rollback': None, 'rollforward': None}


class Refpkg(object):
    """The parts of taxtastic's Refpkg that create, check and info use."""

    def __init__(self, path, create=False):
        self.path = os.path.abspath(path)
        if not os.path.exists(path):
            if create:
                os.mkdir(path)
                self.contents = manifest_template()
                self.save()
            else:
                raise ValueError("Reference package {0} does not exist.".format(path))
        self.load()

    # -- manifest
    def save(self):
        with open(os.path.join(self.path, MANIFEST), 'w') as h:
            json.dump(self.contents, h, indent=4)
            h.write('\n')

    def load(self):
        try:
            fobj = open(os.path.join(self.path, MANIFEST))
        except IOError as e:
            if e.errno == errno.ENOENT:
                raise ValueError("couldn't find manifest file in %s" % (self.path,))
            elif e.errno == errno.ENOTDIR:
                raise ValueError("%s is not a directory" % (self.path,))
            raise
        with fobj:
            self.contents = json.load(fobj)
        self.contents.setdefault('log', [])
        self.contents.setdefault('rollback', None)
        self.contents.setdefault('rollforward', None)
        error = self.is_invalid()
        if error:
            raise ValueError("Refpkg is invalid: %s" % error)

    # -- files
    def resource_path(self, resource):
        if resource not in self.contents['files']:
            raise ValueError("No such resource %r in refpkg" % (resource,))
        return os.path.join(self.path, self.contents['files'][resource])

    def add_file(self, key, path):
        """Refpkg.update_file: copy path in (a unique name if its basename is taken) and record
        its md5 sum."""
        filename = os.path.basename(path)
        base, ext = os.path.splitext(filename)
        if os.path.exists(os.path.join(self.path, filename)):
            with tempfile.NamedTemporaryFile(dir=self.path, prefix=base, suffix=ext) as tf:
                filename = os.path.basename(tf.name)
        shutil.copyfile(path, os.path.join(self.path, filename))
        self.contents['files'][key] = filename
        self.contents['md5'][key] = md5file(path)

    # -- transactions (the state they leave in CONTENTS.json)
    def start_transaction(self):
        self.rollback_state = copy.deepcopy(self.contents)

    def commit_transaction(self, log):
        rollback = self.rollback_state
        rollback.pop('log')
        rollback.pop('rollforward')
        self.contents['log'].insert(0, log)
        self.contents['rollback'] = rollback
        self.contents['rollforward'] = None
        self.save()

    def strip(self):
        """Remove rollback/rollforward and every file the manifest does not name."""
        self.load()
        current = set(self.contents['files'].values())
        to_delete = set(os.listdir(self.path)).difference(current)
        to_delete.discard(MANIFEST)
        for f in to_delete:
            os.unlink(os.path.join(self.path, f))
        self.contents['rollback'] = None
        self.contents['rollforward'] = None
        self.contents['log'].insert(0, 'Stripped refpkg (removed %d files)' % len(to_delete))
        self.save()

    # -- validation
    def is_invalid(self):
        c = self.contents
        for k in ['metadata', 'files', 'md5']:
            if k not in c:
                return "Manifest file missing key %s" % k
            if not isinstance(c[k], dict):
                return "Key %s in manifest did not refer to a dictionary" % k
        if 'rollback' not in c:
            return "Manifest file missing key rollback"
        if not isinstance(c['rollback'], dict) and c["rollback"] is not None:
            return ("Key rollback in manifest did not refer to a dictionary or None, found %s") % str(c['rollback'])
        if 'rollforward' not in c:
            return "Manifest file missing key rollforward"
        if c['rollforward'] is not None:
            if not isinstance(c['rollforward'], list):
                return "Key rollforward was not a list, found %s" % str(c['rollforward'])
            elif len(c['rollforward']) != 2:
                return "Key rollforward had wrong length, found %d" % len(c['rollforward'])
            elif not isinstance(c['rollforward'][0], str):
                return "Key rollforward's first entry was not a string, found %s" % str(c['rollforward'][0])
            elif not isinstance(c['rollforward'][1], dict):
                return "Key rollforward's second entry was not a dict, found %s" % str(c['rollforward'][1])
        if "log" not in c:
            return "Manifest file missing key 'log'"
        if not isinstance(c['log'], list):
            return "Key 'log' in manifest did not refer to a list"
        if c['files'].keys() != c['md5'].keys():
            return ("Files and MD5 sums in manifest do not match (files: %s, MD5 sums: %s)") % \
                (list(c['files'].keys()), list(c['md5'].keys()))
        for key, filename in c['files'].items():
            expected = c['md5'][key]
            found = md5file(self.resource_path(key))
            if found != expected:
                return ("File %s referred to by key %s did not match its MD5 sum (found: %s, expected %s)") % \
                    (filename, key, found, expected)
        return False

    def is_ill_formed(self):
        """Refpkg.is_ill_formed: is_invalid plus the content checks of ``taxit check``. Mismatched
        names are listed sorted (taxit lists them in set order, which changes from run to run)."""
        m = self.is_invalid()
        if m:
            return m
        for k in ('aln_fasta', 'aln_sto', 'seq_info', 'tree', 'taxonomy', 'phylo_model'):
            if k not in self.contents['files']:
                return "RefPkg has no key " + k
        with open(self.resource_path('aln_fasta')) as f:
            if not f.readline().startswith('>'):
                return 'aln_fasta file is not valid FASTA.'
        fasta_names = {name for name, _ in read_fasta(self.resource_path('aln_fasta'))}
        with open(self.resource_path('seq_info')) as f:
            lines = list(csv.reader(f))
            headers = set(lines[0])
            for req_header in 'seqname', 'tax_id':
                if req_header not in headers:
                    return "seq_info is missing {0}".format(req_header)
            if len({len(line) for line in lines}) > 1:
                return "some lines in seq_info differ in field cout"
            csv_names = {line[0] for line in lines[1:]}
        with open(self.resource_path('aln_sto')) as f:
            try:
                sto_names = set(parse_stockholm(f))
            except ValueError:
                return 'aln_sto file is not valid Stockholm.'
        try:
            t = parse_newick(read_text(self.resource_path('tree')))
            tree_names = {t.label[v] for v in t.leaves()}
        except Exception:
            return 'tree file is not valid Newick.'
        d = fasta_names.symmetric_difference(sto_names)
        if len(d) != 0:
            return "Names in aln_fasta did not match aln_sto.  Mismatches: " + ', '.join(sorted(str(x) for x in d))
        d = fasta_names.symmetric_difference(csv_names)
        if len(d) != 0:
            return "Names in aln_fasta did not match seq_info.  Mismatches: " + ', '.join(sorted(str(x) for x in d))
        d = fasta_names.symmetric_difference(tree_names)
        if len(d) != 0:
            return "Names in aln_fasta did not match nodes in tree.  Mismatches: " + ', '.join(sorted(str(x) for x in d))
        with open(self.resource_path('taxonomy')) as f:
            if len({len(line) for line in csv.reader(f)}) > 1:
                return "Taxonomy is invalid: not all lines had the same number of fields."
        with open(self.resource_path('phylo_model')) as f:
            try:
                json.load(f)
            except ValueError:
                return "phylo_model is not valid JSON."
        return False


def _has_rppr(rppr):
    """utils.has_rppr: rppr runs (it exits non-zero without arguments)."""
    try:
        subprocess.run([rppr], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
        raise
    return True


# ------------------------------------------------------------------------------ commands ---

CREATE_FILES = ['aln_fasta', 'aln_sto', 'mask', 'profile', 'seq_info', 'taxonomy', 'tree', 'tree_stats',
                'readme', 'model']


def create(package_name, locus, clobber=False, author=None, description=None, package_version=None,
           stats_type=None, frequency_type=None, reroot=True, rppr='rppr', **files):
    """taxit create; files: the CREATE_FILES keys (and tree_stats). Returns the exit status."""
    if clobber and os.path.isdir(package_name):
        try:
            shutil.rmtree(package_name)
        except Exception:
            _warn("Failed: Could not delete %s" % package_name)
            return 1
    elif clobber and os.path.exists(package_name):
        try:
            os.unlink(package_name)
        except Exception:
            _warn("Failed: Could not delete %s" % package_name)
            return 1
    elif not clobber and os.path.exists(package_name):
        _warn('Failed: {0} exists.'.format(package_name))
        return 1

    r = Refpkg(package_name, create=True)
    r.start_transaction()
    r.contents['metadata']['locus'] = locus
    for key, value in (('description', description), ('author', author), ('package_version', package_version)):
        if value:
            r.contents['metadata'][key] = value
    if files.get('tree_stats'):
        model = phylo_model_of(stats_type, files['tree_stats'], frequency_type)
        tf = tempfile.NamedTemporaryFile(prefix='phylo_model', suffix='.json', delete=False)
        tf.close()
        try:
            with open(tf.name, 'w') as fh:
                json.dump(model, fh, indent=4)
            r.add_file('phylo_model', tf.name)
        finally:
            os.unlink(tf.name)
    for key in CREATE_FILES:
        if files.get(key):
            r.add_file(key, files[key])
    r.commit_transaction('Loaded initial files into empty refpkg')
    r.strip()

    if reroot and files.get('taxonomy') and files.get('seq_info') and files.get('tree'):
        if _has_rppr(rppr):
            r.start_transaction()
            tf = tempfile.NamedTemporaryFile(prefix='tree', suffix='.tre', delete=False)
            tf.close()
            try:
                subprocess.check_call([rppr or 'rppr', 'reroot', '-c', r.path, '-o', tf.name])
                r.add_file('tree', tf.name)
            finally:
                os.unlink(tf.name)
            r.commit_transaction('Rerooted')
        else:
            _warn('"%s" not found. Skipping rerooting' % rppr)
    return 0


def check(refpkg):
    """taxit check: prints the problem and returns 1, or returns 0."""
    r = Refpkg(refpkg)
    msg = r.is_ill_formed()
    if msg:
        print(msg)
        return 1
    return 0


def info(refpkg, seq_names=False, tally=False, lengths=False):
    """taxit info [-n | -t | -l]."""
    pkg = Refpkg(refpkg)
    with open(pkg.resource_path('seq_info')) as seq_info:
        snames = [row['seqname'] for row in csv.DictReader(seq_info)]
    if seq_names:
        print('\n'.join(snames))
    elif tally:
        counts = {}
        with open(pkg.resource_path('taxonomy')) as taxtab, open(pkg.resource_path('seq_info')) as si:
            taxdict = {row['tax_id']: row for row in csv.DictReader(taxtab)}
            tax_ids = [d['tax_id'] for d in csv.DictReader(si)]
        for tax_id in tax_ids:
            counts[tax_id] = counts.get(tax_id, 0) + 1
        rows = [(taxdict[tax_id]['tax_name'], tax_id, count) for tax_id, count in list(counts.items())]
        csv.writer(sys.stdout, quoting=csv.QUOTE_NONNUMERIC).writerows(sorted(rows))
    elif lengths:
        writer = csv.writer(sys.stdout)
        writer.writerow(["seqname", "length"])
        for name, seq in read_fasta(pkg.resource_path('aln_fasta')):
            writer.writerow([name, len(seq.replace('-', ''))])
    else:
        print('number of sequences:', len(snames))
        print('package components\n', '\n'.join(sorted(pkg.contents['files'].keys())))
    return 0
