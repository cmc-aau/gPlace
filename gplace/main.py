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

import gc
import logging
import os
import shutil
import sys
import tempfile
import time
from collections import OrderedDict

from gplace.biolib_lite.common import make_sure_path_exists
from gplace.config.output import PATH_TEST_JPLACE, PATH_TEST_REPORT, PATH_BATCH_SUMMARY
from gplace.exceptions import GPlaceExit, BatchfileMalformed, FileNotFound

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tests', 'data')


class OptionsParser(object):

    def __init__(self, version, output_dir=None):
        """Initialization.

        Parameters
        ----------
        version : str
            The current version number (e.g. 0.1.0).
        """
        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')
        self.version = version

    @staticmethod
    def _check_file(path, what):
        if not path or not os.path.isfile(path):
            raise FileNotFound(f'{what} does not exist: {path}')

    @staticmethod
    def _check_refpkg(path):
        if not path or not os.path.isfile(os.path.join(path, 'CONTENTS.json')):
            raise FileNotFound(f'Not a reference package (no CONTENTS.json): {path}')

    def _place_group(self, refpkg, queries, cpus):
        """Run the engine on one reference package; returns (Placer, results, seconds)."""
        from gplace.engine.placer import Placer
        t0 = time.time()
        placer = Placer(refpkg, queries, cpus)
        cpu_txt = f'{placer.workers} CPUs' if placer.workers else f'--cpus auto on {placer.cpus} CPUs'
        self.logger.info(f'Placing {len(placer.query_names):,} queries on {refpkg} '
                         f'({placer.nnodes:,} nodes, {placer.nsites:,} sites, {cpu_txt}; '
                         f'estimated peak memory {placer.memory_estimate() / 2 ** 30:.1f} GiB).')
        results = placer.run()
        n_warn = sum(1 for _, _, w in results if w)
        if n_warn:
            self.warnings.info(f'{n_warn} queries had locations skipped after a GSL error (as pplacer does).')
        return placer, results, time.time() - t0

    def place(self, options):
        """Place the queries of one alignment on one reference package."""
        from gplace.engine.core import read_fasta
        from gplace.engine.jplace import write_jplace
        self._check_file(options.msa_file, 'Query alignment')
        self._check_refpkg(options.refpkg)
        out_dir = os.path.dirname(os.path.abspath(options.out_file))
        make_sure_path_exists(out_dir)
        queries = read_fasta(options.msa_file)
        self.logger.info(f'Read {len(queries):,} sequences from {options.msa_file}')
        placer, results, secs = self._place_group(options.refpkg, queries, options.cpus)
        write_jplace(options.out_file, placer.jplace_tree(), results, ' '.join(sys.argv))
        self.logger.info(f'Wrote {len(results):,} placements to {options.out_file} ({secs:.1f} s).')

    def place_batch(self, options):
        """Place many alignments; queries sharing a reference package are placed in one pass."""
        from gplace.engine.core import read_fasta
        from gplace.engine.jplace import write_jplace
        self._check_file(options.batchfile, 'Batchfile')
        rows = []
        for i, line in enumerate(open(options.batchfile), 1):
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) != 3:
                raise BatchfileMalformed(f'Line {i} of {options.batchfile} does not have 3 tab-separated '
                                         f'columns (msa_file, refpkg, out_file).')
            rows.append(cols)
        if not rows:
            raise BatchfileMalformed(f'No entries in {options.batchfile}.')
        groups = OrderedDict()
        for msa, refpkg, out in rows:
            self._check_file(msa, 'Query alignment')
            self._check_refpkg(refpkg)
            groups.setdefault(os.path.realpath(refpkg), []).append((msa, refpkg, out))
        make_sure_path_exists(options.out_dir)
        summary = []
        for gi, (key, entries) in enumerate(groups.items(), 1):
            self.logger.info(f'Reference package {gi}/{len(groups)}: {entries[0][1]} '
                             f'({len(entries)} alignment file(s))')
            pooled, owner = [], {}
            for msa, _, out in entries:
                for nm, seq in read_fasta(msa):
                    if nm in owner:
                        raise BatchfileMalformed(f'Query {nm} appears in more than one alignment for '
                                                 f'{entries[0][1]}; names must be unique per reference package.')
                    owner[nm] = out
                    pooled.append((nm, seq))
            placer, results, secs = self._place_group(entries[0][1], pooled, options.cpus)
            tree = placer.jplace_tree()
            by_out = OrderedDict((out, []) for _, _, out in entries)
            for r in results:
                by_out[owner[r[0]]].append(r)
            for msa, refpkg, out in entries:
                make_sure_path_exists(os.path.dirname(os.path.abspath(out)))
                write_jplace(out, tree, by_out[out], ' '.join(sys.argv))
                summary.append((refpkg, msa, out, len(by_out[out]), f'{secs:.1f}'))
            placer.release()
            del placer, results, pooled, owner, by_out
            gc.collect()
        path = os.path.join(options.out_dir, PATH_BATCH_SUMMARY)
        with open(path, 'w') as fh:
            fh.write('refpkg\tmsa_file\tout_file\tqueries\trefpkg_seconds\n')
            for row in summary:
                fh.write('\t'.join(map(str, row)) + '\n')
        self.logger.info(f'Done: {sum(r[3] for r in summary):,} queries in {len(summary)} jplace files. '
                         f'Summary: {path}')

    def compare(self, options):
        """Compare two jplace files."""
        from gplace.engine.jplace import compare_jplace
        self._check_file(options.jplace_a, 'jplace file')
        self._check_file(options.jplace_b, 'jplace file')
        r = compare_jplace(options.jplace_a, options.jplace_b)
        n = max(1, r['genomes'])
        self.logger.info(f'Genomes: {r["genomes"]:,} (only in A: {r["only_a"]}, only in B: {r["only_b"]})')
        self.logger.info(f'Same reference tree string: {r["same_tree"]}')
        self.logger.info(f'Best edge identical:   {r["best_edge"]:,} ({100 * r["best_edge"] / n:.4f}%)')
        self.logger.info(f'Edge list identical:   {r["edge_list"]:,} ({100 * r["edge_list"] / n:.4f}%)')
        self.logger.info(f'All values identical:  {r["identical"]:,} ({100 * r["identical"] / n:.4f}%)')
        self.logger.info('Largest |difference| on the best placement: ' +
                         ', '.join(f'{k} {v:.3g}' for k, v in r['max_delta'].items()))
        for name, a, b in r['examples']:
            self.logger.info(f'  {name}\n      A {a}\n      B {b}')

    def benchmark(self, options):
        """Time placement at several worker counts; estimate memory and time for the full job."""
        import math
        import random
        import numpy as np
        from gplace.engine import resources
        from gplace.engine.core import read_fasta, informative
        self._check_refpkg(options.refpkg)
        n_avail = resources.available_cpus()
        GiB = 2 ** 30

        def slurm(mem_bytes, workers, seconds):
            lines = [f'#SBATCH --cpus-per-task={workers}',
                     f'#SBATCH --mem={math.ceil(mem_bytes * 1.10 / GiB)}G']
            if seconds is not None:
                lines.append(f'#SBATCH --time={resources.slurm_time(seconds * 1.5 + 600)}')
            return lines

        if options.msa_file is None:
            if not options.queries:
                raise GPlaceExit('Give --msa_file to benchmark, or --queries N for a memory estimate only.')
            nn, ni, nl, L, sites = resources.refpkg_dimensions(options.refpkg)
            c = resources.parse_cpus(options.cpus)
            w = (resources.load_tuned_cpus() or resources.heuristic_cpus(n_avail)) if c == 'auto' else c
            mem = resources.estimate_memory(nn, ni, sites, options.queries, nl, L, w)
            self.logger.info(f'Reference package: {nn:,} nodes, {sites:,} informative sites (upper bound).')
            self.logger.info(f'Estimated peak memory for {options.queries:,} queries with {w} workers: '
                             f'{mem / GiB:.1f} GiB')
            for line in slurm(mem, w, None):
                self.logger.info(f'  {line}')
            return

        self._check_file(options.msa_file, 'Query alignment')
        queries = read_fasta(options.msa_file)
        q_full = options.queries or len(queries)
        rnd = random.Random(1)
        sample = queries if len(queries) <= options.max_queries else rnd.sample(queries, options.max_queries)
        cands = sorted({int(x) for x in options.cpus_list.split(',')}) if options.cpus_list \
            else resources.default_candidates(n_avail)
        cands = [c for c in cands if 1 <= c <= n_avail] or [1]
        self.logger.info(f'Benchmarking on {len(sample):,} of {len(queries):,} queries; worker counts '
                         f'{", ".join(map(str, cands))} ({n_avail} CPUs available).')
        from gplace.engine.placer import Placer
        placer = Placer(options.refpkg, sample, resources.heuristic_cpus(n_avail))
        placer.cpus = max(cands)
        t0 = time.time(); placer.cache(); t_cache = time.time() - t0
        # prefilter time = fixed per-edge work + work per query: fit from two sweeps
        n1, n2 = len(sample), max(1, len(sample) // 4)
        w_pref = resources.heuristic_cpus(n_avail)
        t2 = placer.prefilter_seconds(n2, w_pref)
        t1 = placer.prefilter_seconds(n1, w_pref)          # leaves the full sample's scores in place
        per_q = max(0.0, (t1 - t2) / max(n1 - n2, 1))
        t_pref_full = max(t1 - per_q * n1, 0.0) + per_q * q_full
        rows = []
        for c in cands:
            t0 = time.time(); placer.place_all(workers=c); dt = time.time() - t0
            rows.append((c, dt, len(sample) / dt))
            self.logger.info(f'  {c:>4} workers: {dt:7.1f} s  ({len(sample) / dt:6.2f} queries/s)')
        top = max(r[2] for r in rows)                  # differences under 3% are noise:
        best, _, best_rate = min((r for r in rows if r[2] >= 0.97 * top), key=lambda r: r[0])

        # full-job estimates: informative sites of the whole alignment, Q-proportional prefilter
        L = placer.aln_len
        qa = np.frombuffer(''.join(s.upper() for _, s in queries).encode(), np.uint8).reshape(len(queries), L)
        full_mask = informative(qa).any(axis=0)
        nn, ni, nl, _, _ = resources.refpkg_dimensions(options.refpkg)
        sites_full = int((full_mask & self._ref_informative(options.refpkg, L)).sum())
        mem = resources.estimate_memory(nn, ni, sites_full, q_full, nl, L, best)
        t_full = t_cache + t_pref_full + q_full / best_rate
        placer.release()
        self.logger.info(f'Fastest: --cpus {best} ({best_rate:.2f} queries/s; '
                         f'{best_rate / max(r[2] for r in rows if r[0] == max(cands)):.2f}x the speed with {max(cands)}).')
        self.logger.info(f'Estimated for {q_full:,} queries: peak memory {mem / GiB:.1f} GiB, '
                         f'time about {t_full / 60:.0f} min.')
        self.logger.info('Suggested SLURM request (+10% memory, +50% time margin):')
        for line in slurm(mem, best, t_full):
            self.logger.info(f'  {line}')
        if options.out_dir:
            make_sure_path_exists(options.out_dir)
            with open(os.path.join(options.out_dir, 'gplace.benchmark.tsv'), 'w') as fh:
                fh.write('workers\tseconds\tqueries_per_second\n')
                for c, dt, r in rows:
                    fh.write(f'{c}\t{dt:.2f}\t{r:.4f}\n')
        if not options.no_save:
            path = resources.save_tuned_cpus(best, {'refpkg': os.path.abspath(options.refpkg),
                                                    'queries_per_second': round(best_rate, 4)})
            self.logger.info(f'Saved: --cpus auto will use {best} on this machine ({path}).')

    @staticmethod
    def _ref_informative(refpkg, L):
        import json
        import numpy as np
        from gplace.engine.core import read_fasta, informative
        files = json.load(open(os.path.join(refpkg, 'CONTENTS.json')))['files']
        ref = read_fasta(os.path.join(refpkg, files['aln_fasta']))
        RA = np.frombuffer(''.join(s for _, s in ref).encode(), np.uint8).reshape(len(ref), L)
        return informative(RA).any(axis=0)

    def run_test(self, options):
        """Place the bundled test genomes and compare with GTDB-Tk's pplacer output."""
        from gplace.engine.core import read_fasta
        from gplace.engine.jplace import write_jplace, compare_jplace
        out_dir = options.out_dir or tempfile.mkdtemp(prefix='gplace_test_')
        make_sure_path_exists(os.path.join(out_dir, os.path.dirname(PATH_TEST_JPLACE)))
        refpkg = os.path.join(TEST_DATA, 'mini.refpkg')
        queries = read_fasta(os.path.join(TEST_DATA, 'queries.fasta'))
        expected = os.path.join(TEST_DATA, 'expected_pplacer.jplace')
        placer, results, secs = self._place_group(refpkg, queries, options.cpus)
        out = os.path.join(out_dir, PATH_TEST_JPLACE)
        write_jplace(out, placer.jplace_tree(), results, ' '.join(sys.argv))
        r = compare_jplace(out, expected)
        with open(os.path.join(out_dir, PATH_TEST_REPORT), 'w') as fh:
            for k in ('genomes', 'best_edge', 'edge_list', 'identical', 'same_tree'):
                fh.write(f'{k}\t{r[k]}\n')
        self.logger.info(f'{r["identical"]}/{r["genomes"]} genomes: every printed value identical to '
                         f'GTDB-Tk\'s pplacer v1.1.alpha19; reference tree string identical: {r["same_tree"]}.')
        if r['identical'] != r['genomes'] or not r['same_tree'] or r['only_a'] or r['only_b']:
            raise GPlaceExit(f'The test placements differ from pplacer (see {out_dir}).')
        self.logger.info(f'Test passed ({secs:.1f} s). Output: {out_dir}')

    def check_install(self, options):
        """Check the engine is present and bit-exact."""
        from gplace.engine import core, selftest
        from gplace.engine.build import cpu_has_avx2
        ok = True
        self.logger.info(f'CPU AVX2/FMA support: {"yes" if cpu_has_avx2() else "NO"}')
        ok &= cpu_has_avx2()
        path = core.engine_path()
        if not os.path.exists(path):
            raise GPlaceExit(f'Engine not found at {path}. Build it with: gplace build_engine')
        L = core.lib()
        self.logger.info(f'Engine: {path}')
        self.logger.info(f'        {L.xp_engine_version().decode()}')
        n, bad = selftest.check_libm()
        self.logger.info(f'libm vs GTDB-Tk pplacer binary (hypot/exp/log): {n - bad}/{n} identical')
        ok &= bad == 0
        si, se, sv = selftest.check_wag_eigen()
        self.logger.info(f'WAG model vs GTDB-Tk pplacer binary: input {"identical" if si else "DIFFERENT"}, '
                         f'eigenvalues {"identical" if se else "DIFFERENT"}, '
                         f'eigenvectors {"identical" if sv else "DIFFERENT"}')
        ok &= si and se and sv
        nk, fails = selftest.check_kernels()
        self.logger.info(f'SIMD kernels vs scalar references: {nk - len(fails)}/{nk} identical')
        for f in fails[:10]:
            self.logger.error(f'  mismatch: {f}')
        ok &= not fails
        if not ok:
            raise GPlaceExit('The engine is not bit-exact on this system; placements would not be '
                             'guaranteed identical to pplacer.')
        self.logger.info('Installation OK: the engine reproduces GTDB-Tk\'s pplacer arithmetic bit for bit.')

    def build_engine(self, options):
        """Compile libgplace_core.so, then run the installation checks."""
        from gplace.engine.build import build_engine
        build_engine(out_dir=options.engine_dir, gsl_prefix=options.gsl_prefix,
                     gsl_tarball=options.gsl_tarball, work_dir=options.build_dir, jobs=options.jobs)
        if options.engine_dir:
            os.environ['GPLACE_ENGINE_DIR'] = options.engine_dir
        self.check_install(options)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        if options.subparser_name == 'place':
            self.place(options)
        elif options.subparser_name == 'place_batch':
            self.place_batch(options)
        elif options.subparser_name == 'compare':
            self.compare(options)
        elif options.subparser_name == 'benchmark':
            self.benchmark(options)
        elif options.subparser_name == 'test':
            self.run_test(options)
        elif options.subparser_name == 'check_install':
            self.check_install(options)
        elif options.subparser_name == 'build_engine':
            self.build_engine(options)
        else:
            self.logger.error('Unknown gPlace command: "' + options.subparser_name + '"\n')
            sys.exit(1)
        return 0
