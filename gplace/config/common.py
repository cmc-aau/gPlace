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

import logging
import os
from functools import lru_cache

ENGINE_NAME = 'libgplace_core.so'          # 20 states: protein models
ENGINE_NAME_NT = 'libgplace_core_nt.so'    # 4 states: DNA (GTR)
ENGINE_NAMES = {20: ENGINE_NAME, 4: ENGINE_NAME_NT}
# "avx2": AVX2 kernels, x86-64 CPUs with AVX2 and FMA; "neon": NEON kernels, ARM64; "portable":
# scalar kernels, any 64-bit CPU. All give the same results. The SIMD engine of an architecture
# (avx2 on x86-64, neon on ARM64) keeps the names above; the portable one is *_portable.so.
ENGINE_VARIANTS = ('avx2', 'neon', 'portable')


def engine_file(ns, variant='avx2'):
    name = ENGINE_NAMES[ns]
    return name if variant != 'portable' else name.replace('.so', '_portable.so')


class __GPlaceCommonConfig:

    def __init__(self):
        # Internal settings used for logging.
        self.LOG_TASK = 21
        self.logger = logging.getLogger('timestamp')

    @property
    def ENGINE_DIR(self):
        """Directory holding the engine libraries: $GPLACE_ENGINE_DIR, else gplace/engine/lib."""
        raw = os.getenv('GPLACE_ENGINE_DIR')
        if raw:
            return os.path.expandvars(os.path.expanduser(raw))
        return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'engine', 'lib')

    @property
    def BUILD_DIR(self):
        """Scratch directory for building GSL 1.16 and the engine."""
        raw = os.getenv('GPLACE_BUILD_DIR')
        if raw:
            return os.path.expandvars(os.path.expanduser(raw))
        return os.path.join(os.path.expanduser('~'), '.cache', 'gplace', 'build')

    @property
    def CACHE_DIR(self):
        """Per-user cache (saved --cpus tuning): $GPLACE_CACHE_DIR, else ~/.cache/gplace."""
        raw = os.getenv('GPLACE_CACHE_DIR')
        if raw:
            return os.path.expandvars(os.path.expanduser(raw))
        return os.path.join(os.path.expanduser('~'), '.cache', 'gplace')

    @property
    def GSL_URLS(self):
        """Where gsl-1.16.tar.gz is downloaded from, in order (GNU's mirror redirector, GNU's own server,
        kernel.org's GNU mirror); its SHA-256 is checked, so any of them will do."""
        return ('https://ftpmirror.gnu.org/gsl/gsl-1.16.tar.gz',
                'https://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz',
                'https://mirrors.kernel.org/gnu/gsl/gsl-1.16.tar.gz')

    @property
    def GSL_SHA256(self):
        return '73bc2f51b90d2a780e6d266d43e487b3dbd78945dd0b04b14ca5980fe28d2f53'


@lru_cache(maxsize=1)
def __get_config():
    return __GPlaceCommonConfig()


CONFIG = __get_config()
