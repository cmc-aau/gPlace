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

import argparse
import tempfile
from contextlib import contextmanager

from gplace.biolib_lite.custom_help_formatter import ChangeTempAction, CustomHelpFormatter


@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler='resolve', help=desc,
                            formatter_class=CustomHelpFormatter)


@contextmanager
def mutex_group(parser, required):
    group = parser.add_argument_group(f'mutually exclusive '
                                      f'{"required" if required else "optional"} '
                                      f'arguments')
    yield group.add_mutually_exclusive_group(required=required)


@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)


def __msa_file(group, required):
    group.add_argument('--msa_file', required=required,
                       help='FASTA file of aligned query sequences, protein or DNA as the reference '
                            'package (e.g. GTDB-Tk user_msa_file.fasta for one class-level tree)')


def __refpkg(group, required):
    group.add_argument('--refpkg', required=required,
                       help='pplacer reference package directory (e.g. gtdbtk.package.N.refpkg)')


def __out_file(group, required):
    group.add_argument('--out_file', required=required,
                       help='output jplace file (same content as pplacer\'s -o file)')


def __batchfile(group, required):
    group.add_argument('--batchfile', required=required,
                       help='tab-separated file with 3 columns: query alignment, reference package, '
                            'output jplace; alignments sharing a reference package are placed together')


def __out_dir(group, required):
    group.add_argument('--out_dir', type=str, default=None, required=required,
                       help='directory to output files')


def __jplace_a(group):
    group.add_argument('--jplace_a', required=True, help='first jplace file')


def __jplace_b(group):
    group.add_argument('--jplace_b', required=True, help='second jplace file')


def __cpus(group):
    group.add_argument('--cpus', default='1', type=str,
                       help='number of CPUs to use, or "auto": the value saved by `gplace benchmark` '
                            'on this machine, otherwise measured during the run (placement is '
                            'memory-bandwidth bound, so the fastest value is often below the CPU count)')


def __tiling(parser):
    with mutex_group(parser, required=False) as grp:
        grp.add_argument('--tile_edges', type=int, default=None,
                         help='place tile by tile, holding one tile of at most this many edges in '
                              'memory at a time (exact; for reference trees too large for memory, such '
                              'as the full GTDB bacterial tree). Default: no tiling')
        grp.add_argument('--max_memory', type=str, default=None,
                         help='memory gPlace may use, e.g. 64G or 500M (a plain number is GiB). If the '
                              'estimated peak fits, gPlace places without tiles; otherwise it picks the '
                              'largest tile size that fits, from the reference tree and the number of '
                              'queries (chosen separately for each reference package)')


def __fix_gamma(group):
    group.add_argument('--fix_gamma_underflow', action='store_true', default=False,
                       help='gamma-model reference packages only: avoid pplacer\'s rescaling underflow, '
                            'which gives likelihood -inf on large trees. Results then differ from pplacer '
                            '(slightly where pplacer is finite). Default: reproduce pplacer exactly')


def __test_models(group):
    group.add_argument('--models', action='store_true', default=False,
                       help='also place the bundled test sets for the other pplacer models (LG, JTT, '
                            'empirical frequencies, GTR, gamma) and check them against pplacer')


def __queries(group):
    group.add_argument('--queries', type=int, default=None,
                       help='number of queries to plan for (default: all sequences in --msa_file); '
                            'with --refpkg alone, gives the memory estimate without placing anything')


def __max_queries(group):
    group.add_argument('--max_queries', type=int, default=400,
                       help='number of queries sampled from --msa_file for the timing runs')


def __cpus_list(group):
    group.add_argument('--cpus_list', type=str, default=None,
                       help='comma-separated worker counts to time (default: 1/8, 1/4, ... of the available CPUs)')


def __no_save(group):
    group.add_argument('--no_save', action='store_true', default=False,
                       help='do not save the fastest value for `--cpus auto` on this machine')


def __msa_file_opt(group):
    group.add_argument('--msa_file', default=None,
                       help='query alignment to benchmark with (a sample of it is placed)')


def __help(group):
    group.add_argument('-h', '--help', action="help", help="show help message")


def __temp_dir(group):
    group.add_argument('--tmpdir', action=ChangeTempAction, default=tempfile.gettempdir(),
                       help="specify alternative directory for temporary files")


def __gsl_prefix(group):
    group.add_argument('--gsl_prefix', default=None,
                       help='prefix of an existing static, -fPIC GSL 1.16 build (lib/libgsl.a, '
                            'include/gsl); otherwise GSL 1.16 is built from source')


def __gsl_tarball(group):
    group.add_argument('--gsl_tarball', default=None,
                       help='path to gsl-1.16.tar.gz (offline builds); downloaded from ftp.gnu.org if unset')


def __build_dir(group):
    group.add_argument('--build_dir', default=None,
                       help='scratch directory for the build (default: $GPLACE_BUILD_DIR or ~/.cache/gplace/build)')


def __engine_dir(group):
    group.add_argument('--engine_dir', default=None,
                       help='where to write the engine libraries (default: $GPLACE_ENGINE_DIR or the package directory)')


def __jobs(group):
    group.add_argument('--jobs', default=8, type=int, help='parallel make jobs for building GSL')


def get_main_parser():
    # Setup the main, and sub parsers.
    main_parser = argparse.ArgumentParser(
        prog='gplace', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(title="Subcommands", help="Available subcommands", dest='subparser_name')

    with subparser(sub_parsers, 'place', 'Place aligned queries on a reference package '
                                         '(bit-identical to pplacer -c <refpkg>).') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __msa_file(grp, required=True)
            __refpkg(grp, required=True)
            __out_file(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __cpus(grp)
            __fix_gamma(grp)
            __temp_dir(grp)
            __help(grp)
        __tiling(parser)

    with subparser(sub_parsers, 'place_batch', 'Place many alignments, pooling the queries of each '
                                               'reference package into one pass.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __batchfile(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __cpus(grp)
            __fix_gamma(grp)
            __temp_dir(grp)
            __help(grp)
        __tiling(parser)

    with subparser(sub_parsers, 'compare', 'Compare two jplace files placement by placement.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __jplace_a(grp)
            __jplace_b(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __help(grp)

    with subparser(sub_parsers, 'benchmark', 'Find the fastest --cpus on this machine and the memory and '
                                             'time a job needs (e.g. for a SLURM job).') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __refpkg(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __msa_file_opt(grp)
            __queries(grp)
            __max_queries(grp)
            __cpus_list(grp)
            __cpus(grp)
            __out_dir(grp, required=False)
            __no_save(grp)
            __help(grp)
        __tiling(parser)

    with subparser(sub_parsers, 'test', 'Place the bundled test genomes and check the output is '
                                        'identical to GTDB-Tk\'s pplacer.') as parser:
        with arg_group(parser, 'optional arguments') as grp:
            __out_dir(grp, required=False)
            __cpus(grp)
            __test_models(grp)
            __help(grp)
        __tiling(parser)

    with subparser(sub_parsers, 'check_install', 'Check the engine is installed and bit-exact.') as parser:
        with arg_group(parser, 'optional arguments') as grp:
            __help(grp)

    with subparser(sub_parsers, 'build_engine', 'Advanced: rebuild the placement engine (it is compiled '
                                                'automatically when gPlace is installed).') as parser:
        with arg_group(parser, 'optional arguments') as grp:
            __gsl_prefix(grp)
            __gsl_tarball(grp)
            __build_dir(grp)
            __engine_dir(grp)
            __jobs(grp)
            __help(grp)

    return main_parser
