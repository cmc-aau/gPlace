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

"""Build the engine: libgplace_core.so (20 states, protein) and libgplace_core_nt.so (4 states, DNA).
The engine is the Rust crate gplace/engine/rust (cargo); this module compiles the third-party C it
links (GSL 1.16, CORE-MATH, glibc 2.19's hypot) and runs cargo once per variant and state count.
Variants, all with the same results: on x86-64 "avx2" (AVX2 kernels) and "portable" (scalar,
*_portable.so; CPUs without AVX2); on ARM64 "neon" (NEON kernels) and "portable" (scalar,
--cfg xp_force_scalar); elsewhere "portable" only. The loader picks one per CPU (core.engine_variant).

Bit-identity with GTDB-Tk's pplacer requires the same numerical libraries as its statically
linked binary: GSL 1.16 (whose eigen-solver differs from GSL 2.x in the last bits), glibc 2.19's
hypot() (called by GSL; glibc >= 2.35 changed it) and correctly-rounded exp/log/pow/sin (glibc
2.19's are; CORE-MATH provides the same correctly-rounded results, also for GSL's own calls). The C is
compiled with -ffp-contract=off, and Rust never fuses a multiply and an add, so no multiply-add is
fused implicitly.
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

CSRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc')        # third-party C
RUST_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rust')   # the engine (Rust)
RUST_TRIPLES = {'x86_64': 'x86_64-unknown-linux-gnu', 'aarch64': 'aarch64-unknown-linux-gnu'}
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


def undefined_builtins(path):
    """Compiler builtins a library leaves undefined (needs nm; empty if nm is missing)."""
    nm = shutil.which('nm')
    if not nm:
        return []
    out = subprocess.run([nm, '-D', '--undefined-only', path], capture_output=True, text=True).stdout
    return sorted({l.split()[-1] for l in out.splitlines() if l.split() and l.split()[-1].startswith('__builtin_')})


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


def cargo_path():
    """The Rust build tool: $GPLACE_CARGO, else cargo on PATH or in ~/.cargo/bin."""
    cargo = os.environ.get('GPLACE_CARGO') or shutil.which('cargo')
    if not cargo and os.path.exists(os.path.expanduser('~/.cargo/bin/cargo')):
        cargo = os.path.expanduser('~/.cargo/bin/cargo')
    if not cargo:
        raise GPlaceExit('The gPlace engine is written in Rust: install a Rust toolchain (cargo; '
                         'https://rustup.rs or conda-forge "rust") or set GPLACE_CARGO.')
    return cargo


def rust_flags(variant, ns, cc='gcc'):
    """RUSTFLAGS of an engine: the state count (--cfg xp_dna for 4), AVX2 kernels for avx2, the
    scalar kernels for the ARM64 portable engine. Rust never fuses multiplies and adds."""
    flags = []
    if variant == 'avx2':
        flags += ['-C', 'target-feature=+avx2,+fma']
    if variant == 'portable' and target_arch(cc) in ('aarch64', 'arm64'):
        flags += ['--cfg', 'xp_force_scalar']
    if ns == 4:
        flags += ['--cfg', 'xp_dna']
    return flags


def build_rust_engine(out, variant, ns, gsl_prefix, objs, link_args, work_dir, cc='gcc', logger=None, _relinked=False):
    """cargo build of the engine crate (gplace/engine/rust) for one variant and state count,
    linked with GSL 1.16 and the CORE-MATH / hypot objects; the library is copied to out."""
    arch = target_arch(cc)
    target_dir = os.path.join(work_dir, 'rust', f'{arch}-{variant}-{ns}')
    # cargo only relinks when the crate or these variables change: the hash of the linked C (objects
    # and GSL archives) makes it relink after the C is recompiled with other code or another compiler
    deps = hashlib.sha256()
    for path in list(objs) + [os.path.join(gsl_prefix, 'lib', n) for n in ('libgsl.a', 'libgslcblas.a')]:
        with open(path, 'rb') as fh:
            deps.update(fh.read())
    env = dict(os.environ, XP_GSL_PREFIX=gsl_prefix, XP_DEPS_OBJS=' '.join(objs), XP_DEPS_HASH=deps.hexdigest(),
               XP_LINK_ARGS=' '.join(link_args), RUSTFLAGS=' '.join(rust_flags(variant, ns, cc)),
               CARGO_TARGET_DIR=target_dir)
    cmd = [cargo_path(), 'build', '--release', '--offline', '--manifest-path',
           os.path.join(RUST_DIR, 'Cargo.toml')]
    lib = os.path.join(target_dir, 'release', 'libgplace_core.so')
    if arch != target_arch('gcc'):                       # cross-compiling: Rust target + C linker
        triple = RUST_TRIPLES[arch]
        cmd += ['--target', triple]
        env[f'CARGO_TARGET_{triple.upper().replace("-", "_")}_LINKER'] = cc
        lib = os.path.join(target_dir, triple, 'release', 'libgplace_core.so')
    if logger:
        logger.debug(' '.join(cmd))
    p = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise GPlaceExit(f'cargo build failed ({p.returncode}):\n{p.stdout[-4000:]}')
    if undefined_builtins(lib) and not _relinked:
        # a link cached by an earlier build (e.g. before an update of the C, with source timestamps
        # older than that link): link again from scratch, once
        shutil.rmtree(target_dir, ignore_errors=True)
        return build_rust_engine(out, variant, ns, gsl_prefix, objs, link_args, work_dir, cc, logger, _relinked=True)
    tmp = out + '.tmp'
    shutil.copy2(lib, tmp)
    os.replace(tmp, out)


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
    cross = [] if cc == 'gcc' else ['CC=' + cc, '--host=' + subprocess.run(
        [cc, '-dumpmachine'], capture_output=True, text=True).stdout.strip()]
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
    for variant in engine_variants(cc):
        flags = FLAGS + std + variant_flags(variant, cc)
        obj = os.path.join(work_dir, 'obj', f'{target_arch(cc)}-{variant}')
        os.makedirs(obj, exist_ok=True)
        o = lambda n: os.path.join(obj, n)
        _run([cc] + flags + ['-c', '-o', o('cr_exp.o'), os.path.join(cm, 'exp.c')], logger=logger)
        _run([cc] + flags + ['-I' + cm, '-c', '-o', o('cr_log.o'), os.path.join(cm, 'log.c')], logger=logger)
        _run([cc] + flags + ['-I' + os.path.join(cm, 'pow'), '-c', '-o', o('cr_pow.o'),
                             os.path.join(cm, 'pow', 'pow.c')], logger=logger)
        _run([cc] + flags + ['-c', '-o', o('cr_sin.o'), os.path.join(cm, 'sin.c')], logger=logger)
        _run([cc] + FLAGS + std + ['-I' + gl, '-D__ieee754_hypot=hypot', '-c', '-o', o('hypot219.o'),
                             os.path.join(gl, 'e_hypot.c')], logger=logger)
        _run([cc] + FLAGS + std + ['-I' + gl, '-D__log1p=log1p', '-c', '-o', o('log1p219.o'),
                             os.path.join(gl, 's_log1p.c')], logger=logger)
        libm_map = LIBM_MAP if variant == 'avx2' else LIBM_MAP_PORTABLE
        objs = [o('cr_exp.o'), o('cr_log.o'), o('cr_pow.o'), o('cr_sin.o'), o('hypot219.o'), o('log1p219.o')]
        for ns in sorted(ENGINE_NAMES, reverse=True):          # one engine per state count
            target = os.path.join(out_dir, engine_file(ns, variant))
            build_rust_engine(target, variant, ns, gsl_prefix, objs, libm_map, work_dir, cc, logger)
            check_engine(target, cc)
            _say(logger, f'Engine ({variant}, {ns} states) written to {target}')
            targets.append(target)
    return targets
