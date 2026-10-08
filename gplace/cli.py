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
                       help='aligned query sequences, protein or DNA as the reference: FASTA, or '
                            'Stockholm if the name ends in .sto or .sth (.gz: gzipped) as in pplacer '
                            '(e.g. GTDB-Tk user_msa_file.fasta for one class-level tree)')


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


def __ranked_edges(group):
    group.add_argument('--ranked_edges', type=int, default=None, metavar='N',
                       help='untiled placement: keep only each query\'s first N edges of pplacer\'s ranking '
                            '(164 covers pplacer\'s default search) instead of its prefilter score for every '
                            'edge, scoring queries in batches of 20,000. Memory then barely grows with the '
                            'number of queries; results are identical (a query that needs more edges is '
                            'ranked again in full). Tiled placement always works this way. Default: keep '
                            'every score')


def __mrca_class(group):
    group.add_argument('--mrca_class', '--mrca-class', action='store_true', default=False,
                       help='reference packages with a taxonomy: classify each placement by the taxonomic '
                            'MRCA of the reference sequences below its edge, instead of by edge painting '
                            '(pplacer\'s --mrca-class)')


def __discard_nonoverlapped(group):
    group.add_argument('--discard_nonoverlapped', '--discard-nonoverlapped', action='store_true', default=False,
                       help='drop queries that have no residue in any column where a reference sequence '
                            'has one, instead of stopping (pplacer\'s --discard-nonoverlapped)')


def __search(parser):
    """pplacer's search settings, under gPlace's and pplacer's spellings (defaults: pplacer's)."""
    with arg_group(parser, 'pplacer settings (default: pplacer\'s; results stay identical to '
                           'pplacer run with the same settings)') as grp:
        grp.add_argument('--ml_tolerance', '--ml-tolerance', type=float, default=0.01, metavar='FLOAT',
                         help='tolerance of the first-stage branch length optimisation (the kept '
                              'placements are refined to 1e-5)')
        grp.add_argument('--start_pend', '--start-pend', type=float, default=0.1, metavar='FLOAT',
                         help='starting pendant branch length')
        grp.add_argument('--max_pend', '--max-pend', type=float, default=2.0, metavar='FLOAT',
                         help='maximum pendant branch length')
        grp.add_argument('--max_strikes', '--max-strikes', type=int, default=6, metavar='INT',
                         help='stop trying edges after this many strikes (edges more than strike_box below '
                              'the best so far); 0 tries every edge of the tree')
        grp.add_argument('--strike_box', '--strike-box', type=float, default=3.0, metavar='FLOAT',
                         help='size of the strike box in log-likelihood units')
        grp.add_argument('--max_pitches', '--max-pitches', type=int, default=40, metavar='INT',
                         help='try at most this many edges (plus one) per query')
        grp.add_argument('--keep_at_most', '--keep-at-most', type=int, default=7, metavar='INT',
                         help='maximum number of placements kept per query')
        grp.add_argument('--keep_factor', '--keep-factor', type=float, default=0.01, metavar='FLOAT',
                         help='keep placements whose likelihood is at least this fraction of the best '
                              'one\'s')
        grp.add_argument('--no_pre_mask', '--no-pre-mask', action='store_true', default=False,
                         help='do not pre-mask: keep every alignment column (query residues where no '
                              'reference has one then count), do not check that queries overlap the '
                              'reference sequences, and merge only queries with identical whole sequences '
                              'into one record')
        grp.add_argument('--always_refine', '--always-refine', action='store_true', default=False,
                         help='refine the model before placing: CAT models get new site rate categories '
                              'on the reference tree (pplacer always does this when the query alignment '
                              'holds the reference sequences); no effect on gamma models')


def __reference_files(parser):
    """pplacer's reference files (-t/-r/-s/-d), with or without --refpkg."""
    with arg_group(parser, 'reference files (instead of, or replacing those of, --refpkg; pplacer\'s '
                           '-t/-r/-s/-d)') as grp:
        grp.add_argument('--tree', '-t', default=None,
                         help='reference tree (Newick)')
        grp.add_argument('--ref_aln', '-r', default=None,
                         help='reference alignment (FASTA). Without it (and without --refpkg) the query '
                              'alignment must contain every reference sequence')
        grp.add_argument('--stats', '-s', default=None,
                         help='model: a RAxML 7 info file, a PhyML 3.0 stats file or a JSON model file '
                              '(as pplacer reads them; gamma rates). A --refpkg phylo_model takes precedence')
        grp.add_argument('--ref_dir', '-d', default=None,
                         help='directory containing the --tree, --ref_aln and --stats files')
    with arg_group(parser, 'model settings (pplacer\'s; used only for what a --stats file does not set)') as grp:
        grp.add_argument('--model', '-m', default='LG',
                         help='substitution model name. pplacer v1.1.alpha19 always takes it from the '
                              'statistics file or the reference package, so this has no effect (accepted '
                              'for pplacer command lines)')
        grp.add_argument('--model_freqs', '--model-freqs', action='store_true', default=False,
                         help='protein models from a --stats file: the model\'s amino acid frequencies '
                              'instead of the reference alignment\'s')
        grp.add_argument('--gamma_cats', '--gamma-cats', type=int, default=1, metavar='INT',
                         help='gamma categories, when the --stats file gives none (PhyML without gamma, '
                              'JSON with a null ras_model)')
        grp.add_argument('--gamma_alpha', '--gamma-alpha', type=float, default=1.0, metavar='FLOAT',
                         help='gamma shape parameter, when the --stats file gives none')


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
                       help='path to gsl-1.16.tar.gz (offline builds); downloaded from a GNU mirror if unset')


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

    with subparser(sub_parsers, 'place', 'Place aligned queries on a reference package or reference '
                                         'files (bit-identical to pplacer).') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __msa_file(grp, required=True)
            __out_file(grp, required=True)
        with arg_group(parser, 'reference (--refpkg, or --tree with the files below)') as grp:
            __refpkg(grp, required=False)
        with arg_group(parser, 'optional arguments') as grp:
            __cpus(grp)
            __fix_gamma(grp)
            __mrca_class(grp)
            __discard_nonoverlapped(grp)
            __ranked_edges(grp)
            __temp_dir(grp)
            __help(grp)
        __tiling(parser)
        __search(parser)
        __reference_files(parser)

    with subparser(sub_parsers, 'place_batch', 'Place many alignments, pooling the queries of each '
                                               'reference package into one pass.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __batchfile(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __cpus(grp)
            __fix_gamma(grp)
            __mrca_class(grp)
            __discard_nonoverlapped(grp)
            __ranked_edges(grp)
            __temp_dir(grp)
            __help(grp)
        __tiling(parser)
        __search(parser)

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

    with subparser(sub_parsers, 'tog', 'Make a tree with each query placed as a pendant edge (same output '
                                       'as guppy tog).') as parser:
        parser.add_argument('placefiles', nargs='+', metavar='placefile', help='jplace file(s)')
        with arg_group(parser, 'guppy tog options') as grp:
            grp.add_argument('-o', dest='tog_out_file', metavar='FILE', help='Specify the filename to write to.')
            grp.add_argument('--out-dir', '--out_dir', dest='tog_out_dir', metavar='DIR',
                             help='Specify the directory to write files to.')
            grp.add_argument('--prefix', dest='tog_prefix', metavar='PREFIX',
                             help='Specify a string to be prepended to filenames.')
            grp.add_argument('--pp', action='store_true', help='Use posterior probability for the weight.')
            grp.add_argument('--xml', action='store_true', help='Write phyloXML (with colors) for all visualizations.')
            grp.add_argument('--node-numbers', '--node_numbers', dest='node_numbers', action='store_true',
                             help='Put the node numbers in where the bootstraps usually go.')
            __help(grp)

    with subparser(sub_parsers, 'taxit', 'Reference packages for pplacer: taxit create, check and info (same '
                                         'as taxtastic\'s taxit).') as parser:
        taxit_sub = parser.add_subparsers(title='taxit subcommands', dest='taxit_command')
        with subparser(taxit_sub, 'create', 'Create a reference package.') as tp:
            tp.add_argument('-c', '--clobber', action='store_true', dest='clobber', default=False,
                            help='Delete an existing reference package.')
            req = tp.add_argument_group('Required arguments')
            req.add_argument('-P', '--package-name', required=True, dest='package_name', metavar='PATH',
                             help='Name of refpkg to create')
            req.add_argument('-l', '--locus', required=True, dest='locus', metavar='LOCUS',
                             help='The locus described by the reference package')
            meta = tp.add_argument_group('Package Metadata')
            meta.add_argument('-a', '--author', dest='author', metavar='NAME', help='Person who created the reference package')
            meta.add_argument('-d', '--description', dest='description', metavar='TEXT', help='An arbitrary description field')
            meta.add_argument('-r', '--package-version', dest='package_version', metavar='VERSION',
                              help='Release version for the reference package')
            inf = tp.add_argument_group('Input files')
            inf.add_argument('-f', '--aln-fasta', dest='aln_fasta', metavar='FILE', help='Multiple alignment in fasta format')
            inf.add_argument('-i', '--seq-info', dest='seq_info', metavar='file',
                             help='CSV format file describing the aligned reference sequences, minimally containing '
                                  'the fields "seqname" and "tax_id"')
            inf.add_argument('-m', '--mask', dest='mask', metavar='FILE', help='Text file containing a mask')
            inf.add_argument('-M', '--model-file', dest='model', metavar='FILE',
                             help='File containing model information usually the .bestModel file')
            inf.add_argument('-p', '--profile', dest='profile', metavar='FILE', help='Alignment profile')
            inf.add_argument('-R', '--readme', dest='readme', metavar='FILE',
                             help='README file describing the reference package')
            inf.add_argument('-s', '--tree-stats', dest='tree_stats', metavar='FILE',
                             help='File containing tree statistics (for example RAxML_info.whatever")')
            inf.add_argument('-S', '--aln-sto', dest='aln_sto', metavar='FILE', help='Multiple alignment in Stockholm format')
            inf.add_argument('-t', '--tree-file', dest='tree', metavar='FILE', help='Phylogenetic tree in newick format')
            inf.add_argument('-T', '--taxonomy', dest='taxonomy', metavar='FILE',
                             help='CSV format file defining the taxonomy. Fields include "tax_id","parent_id","rank",'
                                  '"tax_name" followed by a column defining tax_id at each rank starting with root')
            ti = tp.add_argument_group('Tree information')
            ti.add_argument('--stats-type', choices=('PhyML', 'FastTree', 'RAxML'),
                            help='stats file type [default: attempt to guess from file contents]')
            ti.add_argument('--frequency-type', choices=('empirical', 'model'),
                            help='Residue frequency type from the model. Required for PhyML Amino Acid alignments.')
            rg = tp.add_argument_group('Taxonomic Rerooting')
            rg.add_argument('--no-reroot', action='store_false', dest='reroot', default=True,
                            help='Do not reroot the reference package using `rppr reroot`. [default: reroot if '
                                 '`rppr` is available and a taxonomy file is specified]')
            rg.add_argument('--rppr', default='rppr', help='Name of the rppr executable. [default: %(default)s]')
        with subparser(taxit_sub, 'check', 'Validate a reference package.') as tp:
            tp.add_argument('refpkg', metavar='REFPKG', help='Path to Refpkg to check')
        with subparser(taxit_sub, 'info', 'Show information about reference packages.') as tp:
            tp.add_argument('refpkg', metavar='refpkg', help='the reference package to operate on')
            tp.add_argument('-n', '--seq-names', action='store_true', default=False, help='print a list of sequence names')
            tp.add_argument('-t', '--tally', action='store_true', default=False,
                            help='print a tally of sequences representing each taxon at rank RANK')
            tp.add_argument('-l', '--lengths', action='store_true', default=False, help='print sequence lengths')

    return main_parser
