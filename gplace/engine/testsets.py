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

"""The bundled test sets (gplace/tests/data), for `gplace test` and the unit tests.

The WAG set is GTDB-Tk's own output. The other models (tests/data/models, see its README) have
expected outputs from the stock pplacer v1.1.alpha19; the mini.refpkg variants share its tree and
alignment and are assembled into work_dir.
"""

import json
import os
from collections import namedtuple

TEST_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests', 'data')
MODELS = os.path.join(TEST_DATA, 'models')
MINI = os.path.join(TEST_DATA, 'mini.refpkg')

TestSet = namedtuple('TestSet', 'name refpkg queries expected fix_gamma tile_edges')

_MINI_VARIANTS = ('LG', 'JTT', 'WAG_emp', 'LG_emp_g1.3x4', 'WAG_g0.5x4', 'JTT_g0.7x8', 'WAG_g0.05x4')


def wag_set():
    return TestSet('WAG (GTDB-Tk)', MINI, os.path.join(TEST_DATA, 'queries.fasta'),
                   os.path.join(TEST_DATA, 'expected_pplacer.jplace'), False, 40)


def _mini_variant(tag, work_dir):
    """mini.refpkg with the phylo_model of models/mini_<tag>.phylo_model.json."""
    out = os.path.join(work_dir, f'mini_{tag}.refpkg')
    os.makedirs(out, exist_ok=True)
    files = json.load(open(os.path.join(MINI, 'CONTENTS.json')))['files']
    contents = {'tree': os.path.join(MINI, files['tree']), 'aln_fasta': os.path.join(MINI, files['aln_fasta']),
                'phylo_model': os.path.join(MODELS, f'mini_{tag}.phylo_model.json')}
    with open(os.path.join(out, 'CONTENTS.json'), 'w') as fh:     # absolute paths: nothing copied
        json.dump({'metadata': {'format_version': '1.1'}, 'files': contents}, fh, indent=2)
    return out


def model_sets(work_dir):
    """The other-model test sets; mini.refpkg variants are written to work_dir."""
    q_mini = os.path.join(TEST_DATA, 'queries.fasta')
    sets = [TestSet(f'mini {tag}', _mini_variant(tag, work_dir), q_mini,
                    os.path.join(MODELS, f'mini_{tag}.pplacer.jplace'), False, 40) for tag in _MINI_VARIANTS]
    sets.append(TestSet('mini WAG_g0.05x4 --fix_gamma_underflow', _mini_variant('WAG_g0.05x4', work_dir), q_mini,
                        os.path.join(MODELS, 'mini_WAG_g0.05x4.fix.jplace'), True, 40))
    q_hiv = os.path.join(MODELS, 'like_gtr_queries.fasta')
    sets += [TestSet('HIV GTR + gamma (pplacer tests)', os.path.join(MODELS, 'like_gtr.refpkg'), q_hiv,
                     os.path.join(MODELS, 'like_gtr.pplacer.jplace'), False, 6),
             TestSet('actin WAG + gamma (pplacer tests)', os.path.join(MODELS, 'like_wag.refpkg'),
                     os.path.join(MODELS, 'like_wag_queries.fasta'),
                     os.path.join(MODELS, 'like_wag.pplacer.jplace'), False, 6),
             TestSet('HIV GTR + CAT (FastTree)', os.path.join(MODELS, 'dna_cat.refpkg'), q_hiv,
                     os.path.join(MODELS, 'dna_cat.pplacer.jplace'), False, 6)]
    return sets
