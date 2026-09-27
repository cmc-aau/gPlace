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

"""Compile the gPlace placement engine (libgplace_core.so) as part of the package build.

Environment variables (all optional):
  GPLACE_GSL_TARBALL   path to gsl-1.16.tar.gz (offline builds; otherwise downloaded, SHA-256 checked)
  GPLACE_GSL_PREFIX    existing static, -fPIC GSL 1.16 build (lib/libgsl.a, include/gsl)
  GPLACE_BUILD_DIR     scratch directory for the GSL build (default ~/.cache/gplace/build; reused)
  GPLACE_BUILD_JOBS    parallel make jobs for GSL (default 8)
"""

import os
import sys

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution

try:
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:                                  # setuptools < 70.1
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

ROOT = os.path.dirname(os.path.abspath(__file__))


class BuildPyWithEngine(build_py):
    def run(self):
        super().run()
        sys.path.insert(0, ROOT)
        from gplace.engine.build import build_engine
        if getattr(self, 'editable_mode', False):
            out_dir = os.path.join(ROOT, 'gplace', 'engine', 'lib')
        else:
            out_dir = os.path.join(self.build_lib, 'gplace', 'engine', 'lib')
        print('[gPlace] Compiling the placement engine (GSL 1.16 + glibc-2.19 hypot + CORE-MATH exp/log)',
              flush=True)
        build_engine(out_dir=out_dir,
                     gsl_prefix=os.environ.get('GPLACE_GSL_PREFIX') or None,
                     gsl_tarball=os.environ.get('GPLACE_GSL_TARBALL') or None,
                     work_dir=os.environ.get('GPLACE_BUILD_DIR') or None,
                     jobs=int(os.environ.get('GPLACE_BUILD_JOBS', '8')))


class BdistWheelPy3(_bdist_wheel):
    """The engine is loaded with ctypes and does not depend on the Python version: tag the
    wheel py3-none-<platform> so one wheel serves every supported Python."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, plat = super().get_tag()
        return 'py3', 'none', plat


class BinaryDistribution(Distribution):
    """The wheel contains a compiled library: tag it for the platform."""

    def has_ext_modules(self):
        return True


setup(cmdclass={'build_py': BuildPyWithEngine, 'bdist_wheel': BdistWheelPy3}, distclass=BinaryDistribution)
