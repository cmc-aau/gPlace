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

"""The reference side of a placement run, as pplacer assembles it (pplacer_run.run_file, Refpkg):
a reference package (-c) and/or files given on the command line (-t tree, -r reference alignment,
-s statistics file), the files prefixed with the -d directory. Files given on the command line
replace the package's.

  * reference alignment: aln_fasta (FASTA), else the package's aln_sto (Stockholm), else none:
    the reference sequences must then be in the query file;
  * model: the package's phylo_model JSON, else a statistics file (tree_stats or -s; stats.py)
    read with the model settings of the command line (--gamma-cats, --gamma-alpha,
    --model-freqs). As in pplacer, those settings and -s do not change a phylo_model.
"""

import logging
import os

from gplace.engine import core
from gplace.exceptions import RefpkgMalformed, FileNotFound
from gplace.fileio import read_json, read_text


class ModelSettings(object):
    """pplacer's model preferences (-m, --model-freqs, --gamma-cats, --gamma-alpha): used only for
    a statistics file, for what the file does not set."""

    def __init__(self, model='LG', model_freqs=False, gamma_cats=1, gamma_alpha=1.0):
        self.model, self.model_freqs = model, bool(model_freqs)
        self.gamma_cats, self.gamma_alpha = int(gamma_cats), float(gamma_alpha)

    @classmethod
    def from_options(cls, options):
        kw = {k: getattr(options, k) for k in ('model', 'model_freqs', 'gamma_cats', 'gamma_alpha')
              if getattr(options, k, None) is not None}
        return cls(**kw)


class Reference(object):
    """Files of one reference: refpkg (directory, optional), tree / ref_aln / stats (paths that
    override the package's), ref_dir (prefix for those three, as pplacer's -d)."""

    def __init__(self, refpkg=None, tree=None, ref_aln=None, stats=None, ref_dir=None, settings=None):
        self.refpkg = refpkg
        self.settings = settings or ModelSettings()
        self.files = {}
        if refpkg:
            cpath = os.path.join(refpkg, 'CONTENTS.json')
            if not os.path.isfile(cpath):
                raise RefpkgMalformed(f'No CONTENTS.json in reference package: {refpkg}')
            for k, v in read_json(cpath).get('files', {}).items():
                self.files[k] = os.path.join(refpkg, v)
        prefix = ''
        if ref_dir:
            if not os.path.isdir(ref_dir):
                raise FileNotFound(f'Reference directory does not exist: {ref_dir}')
            prefix = ref_dir if ref_dir.endswith('/') else ref_dir + '/'
        for key, v in (('tree', tree), ('aln_fasta', ref_aln), ('tree_stats', stats)):
            if v:
                self.files[key] = prefix + v
        if 'tree' not in self.files:
            raise RefpkgMalformed('please specify a reference tree with --tree or --refpkg' if not refpkg else
                                  f'Reference package {refpkg} has no "tree" file.')
        self._pm = None

    @classmethod
    def of(cls, ref):
        """A Reference for a Reference or a reference package directory."""
        return ref if isinstance(ref, Reference) else cls(refpkg=ref)

    def describe(self):
        if self.refpkg:
            return self.refpkg
        return f'tree {self.files["tree"]}'

    def _path(self, key):
        p = self.files[key]
        if not os.path.isfile(p):
            raise FileNotFound(f'{key} file does not exist: {p}' if not self.refpkg else
                               f'Reference package {self.refpkg} has no "{key}" file ({p}).')
        return p

    def has(self, key):
        return key in self.files

    def tree_text(self):
        return read_text(self._path('tree')).strip()

    def alignment(self):
        """The reference alignment [(name, sequence)] in pplacer's order (Stockholm: descending
        names), not uppercased; None if there is none (references in the query file)."""
        if 'aln_fasta' in self.files:
            return core.read_fasta(self._path('aln_fasta'))
        if 'aln_sto' in self.files:
            return core.read_stockholm(self._path('aln_sto'))
        return None

    def phylo_model(self):
        """The model as a phylo_model dict (models.model_from_refpkg)."""
        if self._pm is None:
            if 'phylo_model' in self.files:
                self._pm = read_json(self._path('phylo_model'))
            elif 'tree_stats' in self.files:
                from gplace.engine.stats import phylo_model_of_stats
                if self.refpkg:
                    logging.getLogger('timestamp').warning(
                        'using a statistics file directly is deprecated in pplacer; the reference package '
                        'has no phylo_model, so the model comes from its tree_stats file (as in pplacer).')
                s = self.settings
                self._pm = phylo_model_of_stats(self._path('tree_stats'), not s.model_freqs, s.gamma_cats,
                                                s.gamma_alpha, s.model)
            else:
                raise RefpkgMalformed('please specify a tree model with --stats or --refpkg' if not self.refpkg
                                      else f'Reference package {self.refpkg} has no "phylo_model" or '
                                           f'"tree_stats" file.')
        return self._pm

    def tax_equipped(self):
        return 'taxonomy' in self.files and 'seq_info' in self.files
