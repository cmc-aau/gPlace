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

"""Build libgplace_core.so.

Bit-identity with GTDB-Tk's pplacer requires the same numerical libraries as its statically
linked binary: GSL 1.16 (whose eigen-solver differs from GSL 2.x in the last bits), glibc 2.19's
hypot() (called by GSL; glibc >= 2.35 changed it) and correctly-rounded exp/log (glibc 2.19's
are; CORE-MATH provides the same correctly-rounded results). Every object is compiled with
-ffp-contract=off so no multiply-add is fused implicitly.
"""

import hashlib
import logging
import os
import shutil
import subprocess
import tarfile
import urllib.request

from gplace.config.common import CONFIG, ENGINE_NAME
from gplace.exceptions import GPlaceExit

CSRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc')
FLAGS = ['-O2', '-fPIC', '-ffp-contract=off']
SIMD = ['-mavx2', '-mfma']          # AVX2 kernels; -mfma only for CORE-MATH's explicit fma() calls


def _say(logger, msg):
    if logger and logger.handlers:
        logger.info(msg)
    else:
        print(f'[gPlace] {msg}', flush=True)


def _run(cmd, cwd=None, logger=None):
    if logger:
        logger.debug(' '.join(cmd))
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise GPlaceExit(f'Command failed ({p.returncode}): {" ".join(cmd)}\n{p.stdout[-4000:]}')
    return p.stdout


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def cpu_has_avx2():
    try:
        with open('/proc/cpuinfo') as fh:
            for line in fh:
                if line.startswith('flags'):
                    return ' avx2' in line and ' fma' in line
    except OSError:
        pass
    return False


def gsl_version(prefix):
    h = os.path.join(prefix, 'include', 'gsl', 'gsl_version.h')
    if not os.path.exists(h):
        return None
    for line in open(h):
        if line.startswith('#define GSL_VERSION'):
            return line.split('"')[1]
    return None


def build_gsl116(work_dir, tarball=None, jobs=8, logger=None):
    """Static, position-independent GSL 1.16 in <work_dir>/gsl116; returns the prefix."""
    prefix = os.path.join(work_dir, 'gsl116')
    if gsl_version(prefix) == '1.16' and os.path.exists(os.path.join(prefix, 'lib', 'libgsl.a')):
        return prefix
    os.makedirs(work_dir, exist_ok=True)
    if tarball is None:
        tarball = os.path.join(work_dir, 'gsl-1.16.tar.gz')
        if not os.path.exists(tarball):
            _say(logger, f'Downloading {CONFIG.GSL_URL}')
            urllib.request.urlretrieve(CONFIG.GSL_URL, tarball)
    digest = _sha256(tarball)
    if digest != CONFIG.GSL_SHA256:
        raise GPlaceExit(f'Unexpected SHA-256 for {tarball}: {digest} (expected {CONFIG.GSL_SHA256}).')
    src = os.path.join(work_dir, 'gsl-1.16')
    if os.path.exists(src):
        shutil.rmtree(src)
    with tarfile.open(tarball) as tf:
        tf.extractall(work_dir)
    _say(logger, 'Building GSL 1.16 (static, -fPIC, -ffp-contract=off); this takes a few minutes.')
    env_cflags = 'CFLAGS=' + ' '.join(['-O2', '-ffp-contract=off'])
    _run(['./configure', '--prefix=' + prefix, '--disable-shared', '--with-pic', env_cflags],
         cwd=src, logger=logger)
    _run(['make', f'-j{jobs}'], cwd=src, logger=logger)
    _run(['make', 'install'], cwd=src, logger=logger)
    return prefix


def build_engine(out_dir=None, gsl_prefix=None, gsl_tarball=None, work_dir=None, cc='gcc', jobs=8):
    """Compile libgplace_core.so into out_dir (default: $GPLACE_ENGINE_DIR or gplace/engine/lib)."""
    logger = logging.getLogger('timestamp')
    if not cpu_has_avx2():
        raise GPlaceExit('This CPU does not report AVX2 and FMA support, which the gPlace engine requires.')
    out_dir = out_dir or CONFIG.ENGINE_DIR
    work_dir = work_dir or CONFIG.BUILD_DIR
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    if gsl_prefix:
        ver = gsl_version(gsl_prefix)
        if ver != '1.16':
            raise GPlaceExit(f'GSL at {gsl_prefix} is version {ver}; gPlace needs GSL 1.16 '
                             f'(built static with -fPIC).')
    else:
        gsl_prefix = build_gsl116(work_dir, gsl_tarball, jobs, logger)
    libgsl = os.path.join(gsl_prefix, 'lib', 'libgsl.a')
    libcblas = os.path.join(gsl_prefix, 'lib', 'libgslcblas.a')
    for f in (libgsl, libcblas):
        if not os.path.exists(f):
            raise GPlaceExit(f'Missing {f}: GSL 1.16 must be built with --disable-shared --with-pic.')

    obj = os.path.join(work_dir, 'obj')
    os.makedirs(obj, exist_ok=True)
    cm = os.path.join(CSRC, 'core-math')
    gl = os.path.join(CSRC, 'glibc')
    _run([cc] + FLAGS + SIMD + ['-c', '-o', os.path.join(obj, 'cr_exp.o'), os.path.join(cm, 'exp.c')], logger=logger)
    _run([cc] + FLAGS + SIMD + ['-I' + cm, '-c', '-o', os.path.join(obj, 'cr_log.o'), os.path.join(cm, 'log.c')], logger=logger)
    _run([cc] + FLAGS + ['-I' + gl, '-D__ieee754_hypot=hypot', '-c', '-o', os.path.join(obj, 'hypot219.o'),
                         os.path.join(gl, 'e_hypot.c')], logger=logger)
    target = os.path.join(out_dir, ENGINE_NAME)
    tmp = target + '.tmp'
    _run([cc] + FLAGS + SIMD + ['-shared', '-DXP_CR', '-I' + os.path.join(gsl_prefix, 'include'),
                                '-Wl,-Bsymbolic', '-o', tmp, os.path.join(CSRC, 'xp.c'),
                                os.path.join(obj, 'cr_exp.o'), os.path.join(obj, 'cr_log.o'),
                                os.path.join(obj, 'hypot219.o'), libgsl, libcblas, '-lm'], logger=logger)
    os.replace(tmp, target)
    _say(logger, f'Engine written to {target}')
    return target
