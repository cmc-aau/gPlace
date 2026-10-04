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

"""The right placer for a reference (package): CAT models -> Placer / TiledPlacer, gamma models ->
GammaPlacer / GammaTiledPlacer."""

from gplace.engine.reference import Reference
from gplace.engine.placer import Placer
from gplace.engine.tiled import TiledPlacer


def is_gamma(refpkg_dir):
    """Whether the model of a reference (package directory or reference.Reference) is gamma."""
    return Reference.of(refpkg_dir).phylo_model().get('ras_model') == 'gamma'


def new_placer(refpkg_dir, queries, cpus=1, fix_gamma_underflow=False, discard_nonoverlapped=False,
               search=None, always_refine=False, ranked_edges=None):
    """A loaded placer (untiled) for refpkg_dir, a reference package directory or a
    reference.Reference. fix_gamma_underflow only affects gamma models; discard_nonoverlapped drops
    queries that do not overlap the reference alignment (pplacer's --discard-nonoverlapped) instead
    of stopping; search: core.Search (default: pplacer's); always_refine: refine the model as
    pplacer's --always-refine (CAT models; also done when the query file holds the references);
    ranked_edges: keep each query's first ranked edges instead of a queries x edges score table."""
    refpkg_dir = Reference.of(refpkg_dir)
    if is_gamma(refpkg_dir):
        from gplace.engine.gamma import GammaPlacer
        placer = GammaPlacer(refpkg_dir, queries, cpus, discard_nonoverlapped=discard_nonoverlapped,
                             search=search, always_refine=always_refine, ranked_edges=ranked_edges)
        placer.set_gamma_fix(fix_gamma_underflow)
        return placer
    return Placer(refpkg_dir, queries, cpus, discard_nonoverlapped=discard_nonoverlapped, search=search,
                  always_refine=always_refine, ranked_edges=ranked_edges)


def tiled_placer(placer, tile_edges, **kw):
    """Tile a loaded placer (TiledPlacer.from_placer, keeping its gamma setting)."""
    from gplace.engine.gamma import GammaMixin, GammaTiledPlacer
    if isinstance(placer, GammaMixin):
        tp = GammaTiledPlacer.from_placer(placer, tile_edges, **kw)
        tp.set_gamma_fix(placer.gamma_fix)
        return tp
    return TiledPlacer.from_placer(placer, tile_edges, **kw)
