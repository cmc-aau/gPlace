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

"""Build the engine: libgplace_core.so (20 states, protein) and libgplace_core_nt.so (4 states, DNA),
compiled from gplace/engine/csrc/xp.c and linked with the third-party C it needs (GSL 1.16, built
here; CORE-MATH; glibc 2.19's hypot and log1p). Variants, all with the same results: on x86-64
"avx2" (AVX2 kernels) and "portable" (scalar, *_portable.so; CPUs without AVX2); on ARM64 "neon"
(NEON kernels) and "portable" (scalar, -DXP_FORCE_SCALAR); elsewhere "portable" only. The loader
picks one per CPU (core.engine_variant). Linux and macOS (Apple clang; Apple Silicon and Intel).

Bit-identity with GTDB-Tk's pplacer requires the same numerical libraries as its statically
linked binary: GSL 1.16 (whose eigen-solver differs from GSL 2.x in the last bits), glibc 2.19's
hypot() and log1p() (called by GSL; glibc >= 2.35 changed hypot) and correctly-rounded
exp/log/pow/sin (glibc 2.19's are; CORE-MATH provides the same correctly-rounded results, also for
GSL's own calls). Every object is compiled with -ffp-contract=off, so no multiply-add is fused
implicitly.
"""

import hashlib
import logging
import os
import platform
import shutil
import subprocess
import tarfile
import urllib.request

from gplace.config.common import CONFIG, ENGINE_NAMES, engine_file
from gplace.exceptions import GPlaceExit

CSRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc')   # the engine and third-party C
# -Werror=implicit-function-declaration: a builtin the compiler does not know (e.g. __builtin_roundeven
# before gcc 10) would otherwise compile as a call to a missing function and only fail when the engine
# is loaded
FLAGS = ['-O2', '-fPIC', '-ffp-contract=off', '-Werror=implicit-function-declaration']
SIMD = ['-mavx2', '-mfma']          # AVX2 kernels; -mfma only for CORE-MATH's explicit fma() calls
VARIANT_FLAGS = {'avx2': SIMD, 'neon': [], 'portable': []}    # ARM64: NEON is part of ARMv8-A
# GSL's own libm calls (its gamma functions use exp, log, pow, sin) must give glibc 2.19's results,
# which are correctly rounded: bind them to CORE-MATH inside the library. The other libm functions
# GSL's object files reference are never called by gPlace; they are bound to traps that abort, so
# no inexact system libm function can be used silently (tests/test_engine.py checks the imports).
LIBM_CR = {'exp': 'cr_exp', 'log': 'cr_log', 'pow': 'cr_pow', 'sin': 'cr_sin'}
LIBM_TRAPS = ('cos', 'cosh', 'sinh', 'tan', 'tanh', 'sincos', 'expm1', 'atan', 'atan2', 'asin',
              'acos', 'acosh', 'atanh', 'asinh', 'log10', 'log2', 'exp2', 'cbrt', 'erf', 'erfc', 'lgamma',
              'tgamma', 'fma')     # log1p and hypot: glibc 2.19's own code (csrc/glibc), linked in
LIBM_MAP = ([f'-Wl,--defsym={k}={v}' for k, v in LIBM_CR.items()]
            + [f'-Wl,--defsym={k}=xp_libm_trap_{k}' for k in LIBM_TRAPS])
# macOS: Apple's linker has no --defsym; csrc/libm_bind.c defines these functions inside the engine
# instead (the same bindings and traps; see that file).
# The portable engine is compiled without -mfma, so on x86-64 CORE-MATH's fma() calls go to libm:
# fma is correctly rounded by definition (IEEE 754), so the system's is exact like sqrt. On ARM64
# fma() is an instruction.
LIBM_MAP_PORTABLE = ([f'-Wl,--defsym={k}={v}' for k, v in LIBM_CR.items()]
                     + [f'-Wl,--defsym={k}=xp_libm_trap_{k}' for k in LIBM_TRAPS if k != 'fma'])


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
    """True on an x86-64 CPU that reports AVX2 and FMA (the avx2 engine)."""
    if platform.machine().lower() not in ('x86_64', 'amd64'):
        return False
    if platform.system() == 'Darwin':                    # Intel Macs: the CPU features from sysctl
        sysctl = '/usr/sbin/sysctl' if os.path.exists('/usr/sbin/sysctl') else 'sysctl'
        try:
            out = subprocess.run([sysctl, '-n', 'machdep.cpu.features', 'machdep.cpu.leaf7_features'],
                                 capture_output=True, text=True).stdout.upper().split()
        except OSError:
            return False
        return 'AVX2' in out and 'FMA' in out
    try:
        with open('/proc/cpuinfo') as fh:
            for line in fh:
                if line.startswith('flags'):
                    return ' avx2' in line and ' fma' in line
    except OSError:
        pass
    return False


def target_arch(cc='gcc'):
    """The compiler's target architecture ('x86_64', 'aarch64', ...), from `cc -dumpmachine`."""
    try:
        return subprocess.run([cc, '-dumpmachine'], capture_output=True, text=True).stdout.split('-')[0]
    except OSError:
        return platform.machine()


def target_os(cc='gcc'):
    """'darwin' when the compiler targets macOS, else 'linux'."""
    try:
        triple = subprocess.run([cc, '-dumpmachine'], capture_output=True, text=True).stdout
    except OSError:
        return 'darwin' if platform.system() == 'Darwin' else 'linux'
    return 'darwin' if 'darwin' in triple else 'linux'


def engine_variants(cc='gcc'):
    """Engine variants to build for the compiler's target: avx2 + portable on x86-64, neon +
    portable on ARM64, else portable."""
    arch = target_arch(cc)
    if arch in ('x86_64', 'amd64'):
        return ('avx2', 'portable')
    if arch in ('aarch64', 'arm64'):
        return ('neon', 'portable')
    return ('portable',)


def compiler_id(cc='gcc'):
    """(name, major version) of the C compiler: ('gcc', 9), ('clang', 17), or (None, None)."""
    try:
        out = subprocess.run([cc, '-dM', '-E', '-x', 'c', os.devnull], capture_output=True, text=True).stdout
    except OSError:
        return None, None
    macros = dict(l.split()[1:3] for l in out.splitlines() if l.startswith('#define ') and len(l.split()) >= 3)
    if '__clang_major__' in macros:
        return 'clang', int(macros['__clang_major__'])
    if '__GNUC__' in macros:
        return 'gcc', int(macros['__GNUC__'])
    return None, None


def c_std_flags(cc='gcc'):
    """gcc < 5 compiles C as gnu90 by default; CORE-MATH needs C99 or later."""
    name, major = compiler_id(cc)
    return ['-std=gnu11'] if name == 'gcc' and major < 5 else []


def undefined_symbols(path):
    """Symbols a library imports (needs nm; empty if nm is missing), without symbol versions."""
    nm = shutil.which('nm')
    if not nm:
        return []
    if platform.system() == 'Darwin':                    # Mach-O: nm -u, C names with a leading _
        out = subprocess.run([nm, '-u', path], capture_output=True, text=True).stdout
        return [l.strip()[1:] for l in out.splitlines() if l.strip().startswith('_')]
    out = subprocess.run([nm, '-D', '--undefined-only', path], capture_output=True, text=True).stdout
    return [l.split()[-1].split('@')[0] for l in out.splitlines() if l.split()]


def undefined_builtins(path):
    """Compiler builtins a library leaves undefined (needs nm; empty if nm is missing)."""
    return sorted({s for s in undefined_symbols(path) if s.startswith('__builtin_')})


def check_engine(path, cc='gcc'):
    """Stop if a compiled engine cannot be used: an undefined __builtin_* symbol (a compiler builtin
    the C compiler did not know), or, when built for this machine, any load error."""
    name, major = compiler_id(cc)
    who = f'{name} {major}' if name else cc
    missing = undefined_builtins(path)
    if missing:
        raise GPlaceExit(f'The engine {path} needs {", ".join(missing)}, which the C compiler ({who}) does not '
                         f'provide. Please report this; a newer compiler (e.g. gcc >= 10) avoids it.')
    if target_arch(cc) == platform.machine():
        import ctypes
        try:
            ctypes.CDLL(path)
        except OSError as e:
            raise GPlaceExit(f'The engine was built with {who} but cannot be loaded: {e}')


def variant_flags(variant, cc='gcc'):
    """Compiler flags of a variant; the ARM64 portable engine turns NEON off (XP_FORCE_SCALAR)."""
    if variant == 'portable' and target_arch(cc) in ('aarch64', 'arm64'):
        return ['-DXP_FORCE_SCALAR']
    return VARIANT_FLAGS[variant]


def link_engine(out, variant, ns, gsl_prefix, objs, link_args, cc='gcc', logger=None):
    """Compile the engine (xp.c) for one variant and state count, linked with GSL 1.16 and the
    CORE-MATH / glibc 2.19 objects, into the shared library out."""
    tmp = out + '.tmp'
    bind = [] if target_os(cc) == 'darwin' else ['-Wl,-Bsymbolic']     # macOS binds within the image anyway
    _run([cc] + FLAGS + c_std_flags(cc) + variant_flags(variant, cc)
         + ['-shared', '-DXP_CR', f'-DNS={ns}', '-I' + os.path.join(gsl_prefix, 'include')] + bind
         + list(link_args) + ['-o', tmp, os.path.join(CSRC, 'xp.c')] + list(objs)
         + [os.path.join(gsl_prefix, 'lib', 'libgsl.a'), os.path.join(gsl_prefix, 'lib', 'libgslcblas.a'), '-lm'],
         logger=logger)
    os.replace(tmp, out)


def download_gsl(tarball, logger=None, urls=None, sha256=None, timeout=60, rounds=2):
    """Download gsl-1.16.tar.gz to tarball from the first mirror (CONFIG.GSL_URLS) that answers with
    the expected SHA-256, trying the list `rounds` times (one server is sometimes unreachable)."""
    urls, sha256 = urls or CONFIG.GSL_URLS, sha256 or CONFIG.GSL_SHA256
    part, errors = tarball + '.part', []
    for _ in range(rounds):
        for url in urls:
            _say(logger, f'Downloading {url}')
            try:
                with urllib.request.urlopen(url, timeout=timeout) as response, open(part, 'wb') as fh:
                    shutil.copyfileobj(response, fh)
            except OSError as e:                         # URLError, timeouts, refused connections
                errors.append(f'{url}: {getattr(e, "reason", e)}')
                continue
            digest = _sha256(part)
            if digest == sha256:
                os.replace(part, tarball)
                return tarball
            errors.append(f'{url}: unexpected SHA-256 {digest}')
    if os.path.exists(part):
        os.remove(part)
    raise GPlaceExit('Could not download gsl-1.16.tar.gz:\n  ' + '\n  '.join(errors) + '\nDownload it from one of '
                     'these addresses yourself and set GPLACE_GSL_TARBALL=/path/to/gsl-1.16.tar.gz (gplace '
                     'build_engine: --gsl_tarball).')


def gsl_version(prefix):
    h = os.path.join(prefix, 'include', 'gsl', 'gsl_version.h')
    if not os.path.exists(h):
        return None
    with open(h) as fh:
        for line in fh:
            if line.startswith('#define GSL_VERSION'):
                return line.split('"')[1]
    return None


def build_gsl116(work_dir, tarball=None, jobs=8, logger=None, cc='gcc'):
    """Static, position-independent GSL 1.16 in <work_dir>/gsl116; returns the prefix. A cross
    compiler (cc for another target) builds GSL for that target."""
    host = target_arch(cc)
    prefix = os.path.join(work_dir, 'gsl116' if host == target_arch('gcc') else f'gsl116-{host}')
    if gsl_version(prefix) == '1.16' and os.path.exists(os.path.join(prefix, 'lib', 'libgsl.a')):
        return prefix
    os.makedirs(work_dir, exist_ok=True)
    if tarball is None:
        tarball = os.path.join(work_dir, 'gsl-1.16.tar.gz')
        if not os.path.exists(tarball):
            download_gsl(tarball, logger)
    digest = _sha256(tarball)
    if digest != CONFIG.GSL_SHA256:
        raise GPlaceExit(f'Unexpected SHA-256 for {tarball}: {digest} (expected {CONFIG.GSL_SHA256}).')
    src = os.path.join(work_dir, 'gsl-1.16')
    if os.path.exists(src):
        shutil.rmtree(src)
    with tarfile.open(tarball) as tf:
        tf.extractall(work_dir)
    _say(logger, 'Building GSL 1.16 (static, -fPIC, -ffp-contract=off); this takes a few minutes.')
    # GSL's configure checks call exit() without <stdlib.h>: gcc >= 14 and Apple clang reject that
    # by default, which would turn its IEEE checks off. On macOS (no -Bsymbolic): -fvisibility=hidden
    # keeps GSL's calls to its own functions and BLAS inside the engine (another BLAS loaded in the
    # process, e.g. numpy's, must never answer them), and -fno-builtin keeps every libm call a call
    # (bound to CORE-MATH, libm_bind.c): clang would otherwise fold constant arguments with the
    # system's libm or merge sin and cos into one call.
    cflags = ['-O2', '-ffp-contract=off', '-Wno-error=implicit-function-declaration']
    if target_os(cc) == 'darwin':
        cflags += ['-fvisibility=hidden', '-fno-builtin']
    env_cflags = 'CFLAGS=' + ' '.join(cflags)
    cross = [] if cc == 'gcc' else ['CC=' + cc, '--host=' + subprocess.run(
        [cc, '-dumpmachine'], capture_output=True, text=True).stdout.strip()]
    if target_os(cc) == 'darwin' and host in ('arm64', 'aarch64') and not cross:
        cross = ['--build=aarch64-apple-darwin']          # GSL 1.16's config.guess predates Apple Silicon
    _run(['./configure', '--prefix=' + prefix, '--disable-shared', '--with-pic', env_cflags] + cross,
         cwd=src, logger=logger)
    _run(['make', f'-j{jobs}'], cwd=src, logger=logger)
    _run(['make', 'install'], cwd=src, logger=logger)
    return prefix


def build_engine(out_dir=None, gsl_prefix=None, gsl_tarball=None, work_dir=None, cc='gcc', jobs=8):
    """Compile the engines into out_dir (default: $GPLACE_ENGINE_DIR or gplace/engine/lib): for every
    variant (engine_variants) one library per state count. Returns the list of paths."""
    logger = logging.getLogger('timestamp')
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
        gsl_prefix = build_gsl116(work_dir, gsl_tarball, jobs, logger, cc)
    libgsl = os.path.join(gsl_prefix, 'lib', 'libgsl.a')
    libcblas = os.path.join(gsl_prefix, 'lib', 'libgslcblas.a')
    for f in (libgsl, libcblas):
        if not os.path.exists(f):
            raise GPlaceExit(f'Missing {f}: GSL 1.16 must be built with --disable-shared --with-pic.')

    cm = os.path.join(CSRC, 'core-math')
    gl = os.path.join(CSRC, 'glibc')
    targets = []
    std = c_std_flags(cc)
    darwin = target_os(cc) == 'darwin'
    for variant in engine_variants(cc):
        flags = FLAGS + std + variant_flags(variant, cc)
        # macOS: the third-party functions are hidden inside the engine, so every call to them (GSL's
        # to hypot and log1p included) binds to these definitions, never to the system's
        dep = flags + (['-fvisibility=hidden'] if darwin else [])
        dep0 = FLAGS + std + (['-fvisibility=hidden'] if darwin else [])
        obj = os.path.join(work_dir, 'obj', f'{target_arch(cc)}-{variant}')
        os.makedirs(obj, exist_ok=True)
        o = lambda n: os.path.join(obj, n)
        _run([cc] + dep + ['-c', '-o', o('cr_exp.o'), os.path.join(cm, 'exp.c')], logger=logger)
        _run([cc] + dep + ['-I' + cm, '-c', '-o', o('cr_log.o'), os.path.join(cm, 'log.c')], logger=logger)
        _run([cc] + dep + ['-I' + os.path.join(cm, 'pow'), '-c', '-o', o('cr_pow.o'),
                             os.path.join(cm, 'pow', 'pow.c')], logger=logger)
        _run([cc] + dep + ['-c', '-o', o('cr_sin.o'), os.path.join(cm, 'sin.c')], logger=logger)
        _run([cc] + dep0 + ['-I' + gl, '-D__ieee754_hypot=hypot', '-c', '-o', o('hypot219.o'),
                             os.path.join(gl, 'e_hypot.c')], logger=logger)
        _run([cc] + dep0 + ['-I' + gl, '-D__log1p=log1p', '-c', '-o', o('log1p219.o'),
                             os.path.join(gl, 's_log1p.c')], logger=logger)
        libm_map = LIBM_MAP if variant == 'avx2' else LIBM_MAP_PORTABLE
        objs = [o('cr_exp.o'), o('cr_log.o'), o('cr_pow.o'), o('cr_sin.o'), o('hypot219.o'), o('log1p219.o')]
        if darwin:                                     # no --defsym: the bindings are defined in the engine
            _run([cc] + flags + ['-fno-builtin'] + (['-DXP_TRAP_FMA'] if variant == 'avx2' else [])
                 + ['-c', '-o', o('libm_bind.o'), os.path.join(CSRC, 'libm_bind.c')], logger=logger)
            objs.append(o('libm_bind.o'))
            libm_map = []
        for ns in sorted(ENGINE_NAMES, reverse=True):          # one engine per state count
            target = os.path.join(out_dir, engine_file(ns, variant))
            link_engine(target, variant, ns, gsl_prefix, objs, libm_map, cc, logger)
            check_engine(target, cc)
            _say(logger, f'Engine ({variant}, {ns} states) written to {target}')
            targets.append(target)
    return targets
