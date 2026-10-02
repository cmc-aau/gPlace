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

"""The right placer for a reference package: CAT models -> Placer / TiledPlacer, gamma models ->
GammaPlacer / GammaTiledPlacer."""

from gplace.engine import models
from gplace.engine.placer import Placer
from gplace.engine.tiled import TiledPlacer


def is_gamma(refpkg_dir):
    return models.read_phylo_model(refpkg_dir).get('ras_model') == 'gamma'


def new_placer(refpkg_dir, queries, cpus=1, fix_gamma_underflow=False, discard_nonoverlapped=False):
    """A loaded placer for refpkg_dir (untiled). fix_gamma_underflow only affects gamma models;
    discard_nonoverlapped drops queries that do not overlap the reference alignment (pplacer's
    --discard-nonoverlapped) instead of stopping."""
    if is_gamma(refpkg_dir):
        from gplace.engine.gamma import GammaPlacer
        placer = GammaPlacer(refpkg_dir, queries, cpus, discard_nonoverlapped=discard_nonoverlapped)
        placer.set_gamma_fix(fix_gamma_underflow)
        return placer
    return Placer(refpkg_dir, queries, cpus, discard_nonoverlapped=discard_nonoverlapped)


def tiled_placer(placer, tile_edges, **kw):
    """Tile a loaded placer (TiledPlacer.from_placer, keeping its gamma setting)."""
    from gplace.engine.gamma import GammaMixin, GammaTiledPlacer
    if isinstance(placer, GammaMixin):
        tp = GammaTiledPlacer.from_placer(placer, tile_edges, **kw)
        tp.set_gamma_fix(placer.gamma_fix)
        return tp
    return TiledPlacer.from_placer(placer, tile_edges, **kw)
