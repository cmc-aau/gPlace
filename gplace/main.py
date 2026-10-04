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

    @staticmethod
    def _tiling_options(options):
        """(tile_edges, max_memory in bytes) from --tile_edges / --max_memory (None if unset)."""
        from gplace.engine import resources
        tile_edges = getattr(options, 'tile_edges', None)
        raw = getattr(options, 'max_memory', None)
        if tile_edges is not None and tile_edges < 2:
            raise GPlaceExit('--tile_edges must be at least 2.')
        try:
            max_memory = None if raw is None else resources.parse_memory(raw)
        except ValueError as e:
            raise GPlaceExit(f'--max_memory: {e}')
        return tile_edges, max_memory

    def _search_options(self, options):
        """core.Search from the command line, checked as pplacer checks its preferences."""
        from gplace.engine.core import Search
        search = Search.from_options(options)
        changed = search.changed()
        if changed:
            self.logger.info('pplacer search settings: ' + ', '.join(f'{k}={v}' for k, v in changed)
                             + ' (the others are pplacer\'s defaults).')
        if search.no_pre_mask and getattr(options, 'discard_nonoverlapped', False):
            self.logger.info('--discard_nonoverlapped has no effect with --no_pre_mask (as in pplacer).')
        return search

    def _log_plan(self, plan, max_memory):
        """Report how --max_memory was applied (see tiled.plan_tiling)."""
        GiB = 2 ** 30
        if max_memory is None:
            return
        if plan.tile_edges is None:
            self.logger.info(f'--max_memory {max_memory / GiB:.1f} GiB: the untiled placement fits '
                             f'(estimated peak {plan.untiled_memory / GiB:.1f} GiB); no tiles needed.')
        else:
            self.logger.info(f'--max_memory {max_memory / GiB:.1f} GiB: the untiled placement would need '
                             f'about {plan.untiled_memory / GiB:.1f} GiB; using --tile_edges '
                             f'{plan.tile_edges:,} (estimated peak {plan.memory / GiB:.1f} GiB).')

    def _place_group(self, refpkg, queries, cpus, tile_edges=None, max_memory=None, fix_gamma=False,
                     mrca_class=False, discard_nonoverlapped=False, search=None, always_refine=False,
                     ranked_edges=None, fit_tile_edges=None):
        """Run the engine on one reference package; returns (Placer, results, seconds, classifier):
        results has one entry per query (in input order; identical sequences are placed once),
        classifier is None unless the refpkg has a taxonomy (see taxonomy.Classifier).
        fit_tile_edges: without --tile_edges / --max_memory, place on tiles of at most this many
        edges when the untiled estimate exceeds the memory available (gplace test)."""
        import math
        from gplace.engine import resources, taxonomy
        from gplace.engine.factory import new_placer, tiled_placer
        from gplace.engine.tiled import plan_for_placer
        from gplace.engine.reference import Reference
        t0 = time.time()
        reference = Reference.of(refpkg)
        placer = new_placer(reference, queries, cpus, fix_gamma, discard_nonoverlapped, search, always_refine,
                            ranked_edges)
        if placer.discarded:
            self.logger.info(f'Discarded {len(placer.discarded):,} queries that do not overlap the reference '
                             f'sequences (--discard_nonoverlapped), e.g. {placer.discarded[0]}.')
        n_all, n_rep = len(placer.all_names), len(placer.query_names)
        if n_rep < n_all:
            self.logger.info(f'{n_all:,} queries, {n_rep:,} distinct sequences'
                             f'{"" if placer.search.no_pre_mask else " after pre-masking"}: identical '
                             f'sequences are placed once and written as one record (as pplacer does).')
        classifier = None
        if reference.refpkg and reference.tax_equipped():
            classifier = taxonomy.Classifier(reference.refpkg, placer.tree, mrca_class)
            self.logger.info(f'Reference package has a taxonomy: placements are classified by '
                             f'{classifier.method} ("classification" field).')
            for rank in classifier.missing_ranks:
                self.logger.warning(f'rank {rank} not represented in the lineage of any sequence in '
                                    f'reference package {reference.describe()}.')
        elif mrca_class:
            self.logger.info('--mrca_class has no effect: the reference has no taxonomy.')
        model = placer.model
        self.logger.info(f'Model: {model.describe()} ({model.alphabet.name})'
                         + (' [--fix_gamma_underflow]' if fix_gamma and model.ras == 'gamma' else ''))
        if fix_gamma and model.ras != 'gamma':
            self.logger.info('--fix_gamma_underflow has no effect: the model is not gamma.')
        if fit_tile_edges and tile_edges is None and max_memory is None:
            avail, need = resources.available_memory(), placer.memory_estimate()
            if avail and need > avail:
                self.logger.info(f'Untiled this would need about {need / 2 ** 30:.1f} GiB, more than the '
                                 f'{avail / 2 ** 30:.1f} GiB available: placing tile by tile (tiles of at most '
                                 f'{fit_tile_edges:,} edges; identical results).')
                tile_edges = fit_tile_edges
        plan = plan_for_placer(placer, tile_edges, max_memory)
        self._log_plan(plan, max_memory)
        if ranked_edges is not None and plan.tile_edges is None:
            self.logger.info(f'Keeping the first {min(ranked_edges, placer.nnodes - 1):,} ranked edges of each '
                             f'query (--ranked_edges); queries are scored in batches of '
                             f'{resources.QUERY_BATCH:,}.')
        if plan.tile_edges is not None:
            placer = tiled_placer(placer, plan.tile_edges)
        cpu_txt = f'{placer.workers} CPUs' if placer.workers else f'--cpus auto on {placer.cpus} CPUs'
        tile_txt = (f'{len(placer.tiles):,} tiles of up to {placer.max_edges:,} edges; '
                    if plan.tile_edges is not None else '')
        mem = placer.memory_estimate()
        self.logger.info(f'Placing {len(placer.query_names):,} queries on {reference.describe()} '
                         f'({placer.nnodes:,} nodes, {placer.nsites:,} sites, {cpu_txt}; {tile_txt}'
                         f'estimated peak memory {mem / 2 ** 30:.1f} GiB).')
        avail = resources.available_memory() if tile_edges is None and max_memory is None else None
        if avail and mem > avail:
            self.logger.warning(f'The estimated peak memory ({mem / 2 ** 30:.1f} GiB) is more than the '
                                f'{avail / 2 ** 30:.1f} GiB available; use --max_memory to place tile by '
                                f'tile within a memory limit.')
        results = placer.expand(placer.run())
        unplaced = [r[0] for r in results if not r[1]]
        if unplaced:                            # pplacer_run: "not placed and will be omitted"
            msg = (f'{len(unplaced)} queries have no placement left after --keep_at_most/--keep_factor and '
                   f'are omitted from the jplace (as pplacer does), e.g. {unplaced[0]}.')
            self.logger.warning(msg)
            self.warnings.info(msg)
            results = [r for r in results if r[1]]
        n_warn = sum(1 for _, _, w in results if w)
        if n_warn:
            self.warnings.info(f'{n_warn} queries had locations skipped after a GSL error (as pplacer does).')
        n_inf = sum(1 for _, pl, _ in results if any(not math.isfinite(p[3]) for p in pl))
        if n_inf:
            msg = (f'{n_inf} queries have a -inf likelihood (pplacer: "encountered -inf for final '
                   f'likelihood"); their placements are meaningless.')
            if model.ras == 'gamma' and not fix_gamma:
                msg += (' With gamma models this is pplacer\'s rescaling underflow on large trees; '
                        '--fix_gamma_underflow avoids it (results then differ from pplacer).')
            self.logger.warning(msg)
            self.warnings.info(msg)
        return placer, results, time.time() - t0, classifier

    def _reference(self, options):
        """reference.Reference from --refpkg and/or pplacer's -t/-r/-s/-d files."""
        from gplace.engine.reference import Reference, ModelSettings
        if not options.refpkg and not options.tree:
            raise FileNotFound('Give a reference package (--refpkg) or a reference tree (--tree, with '
                               '--stats and optionally --ref_aln), as pplacer needs -c or -t.')
        if options.refpkg:
            self._check_refpkg(options.refpkg)
        settings = ModelSettings.from_options(options)
        ref = Reference(options.refpkg, options.tree, options.ref_aln, options.stats, options.ref_dir, settings)
        for what, key, given in (('--tree', 'tree', options.tree), ('--ref_aln', 'aln_fasta', options.ref_aln),
                                 ('--stats', 'tree_stats', options.stats)):
            if given:
                self._check_file(ref.files[key], f'{what} file')
        given = [o for o, v, d in (('--model', settings.model, 'LG'), ('--model_freqs', settings.model_freqs, False),
                                   ('--gamma_cats', settings.gamma_cats, 1), ('--gamma_alpha', settings.gamma_alpha, 1.0))
                 if v != d]
        if ref.has('phylo_model'):
            if given or options.stats:
                self.logger.info(f'{", ".join(given + (["--stats"] if options.stats else []))}: no effect, the '
                                 f'reference package\'s phylo_model sets the model (as in pplacer).')
        elif '--model' in given:
            self.logger.info('--model: no effect, the statistics file names the model (as in pplacer v1.1.alpha19).')
        return ref

    def place(self, options):
        """Place the queries of one alignment on one reference (package)."""
        from gplace.engine.core import read_alignment
        from gplace.engine.jplace import write_jplace
        tiling = self._tiling_options(options)
        search = self._search_options(options)
        self._check_file(options.msa_file, 'Query alignment')
        reference = self._reference(options)
        out_dir = os.path.dirname(os.path.abspath(options.out_file))
        make_sure_path_exists(out_dir)
        queries = read_alignment(options.msa_file)
        self.logger.info(f'Read {len(queries):,} sequences from {options.msa_file}')
        placer, results, secs, classify = self._place_group(
            reference, queries, options.cpus, *tiling, fix_gamma=options.fix_gamma_underflow,
            mrca_class=options.mrca_class, discard_nonoverlapped=options.discard_nonoverlapped,
            search=search, always_refine=options.always_refine, ranked_edges=options.ranked_edges)
        records = placer.records([r[0] for r in results])
        write_jplace(options.out_file, placer.jplace_tree(), results, ' '.join(sys.argv), records, classify)
        self.logger.info(f'Wrote {len(results):,} placements ({len(records):,} records) to '
                         f'{options.out_file} ({secs:.1f} s).')

    def place_batch(self, options):
        """Place many alignments; queries sharing a reference package are placed in one pass."""
        from gplace.engine.core import read_alignment
        from gplace.engine.jplace import write_jplace
        tiling = self._tiling_options(options)
        search = self._search_options(options)
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
                for nm, seq in read_alignment(msa):
                    if nm in owner:
                        raise BatchfileMalformed(f'Query {nm} appears in more than one alignment for '
                                                 f'{entries[0][1]}; names must be unique per reference package.')
                    owner[nm] = out
                    pooled.append((nm, seq))
            placer, results, secs, classify = self._place_group(
                entries[0][1], pooled, options.cpus, *tiling, fix_gamma=options.fix_gamma_underflow,
                mrca_class=options.mrca_class, discard_nonoverlapped=options.discard_nonoverlapped,
                search=search, always_refine=options.always_refine, ranked_edges=options.ranked_edges)
            tree = placer.jplace_tree()
            by_out = OrderedDict((out, []) for _, _, out in entries)
            for r in results:
                by_out[owner[r[0]]].append(r)
            for msa, refpkg, out in entries:
                make_sure_path_exists(os.path.dirname(os.path.abspath(out)))
                write_jplace(out, tree, by_out[out], ' '.join(sys.argv),
                             placer.records([r[0] for r in by_out[out]]), classify)
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
        from gplace.engine.core import read_alignment, informative
        from gplace.engine.tiled import plan_tiling, load_tree
        self._check_refpkg(options.refpkg)
        tile_edges, max_memory = self._tiling_options(options)
        n_avail = resources.available_cpus()
        GiB = 2 ** 30

        def slurm(mem_bytes, workers, seconds):
            mem_gib = math.ceil(max_memory / GiB) if max_memory else math.ceil(mem_bytes * 1.10 / GiB)
            lines = [f'#SBATCH --cpus-per-task={workers}', f'#SBATCH --mem={mem_gib}G']
            if seconds is not None:
                lines.append(f'#SBATCH --time={resources.slurm_time(seconds * 1.5 + 600)}')
            return lines

        width = resources.refpkg_width(options.refpkg)

        def job_plan(sites, n_queries, n_leaves, L, workers, log=True):
            plan = plan_tiling(load_tree(options.refpkg), sites, n_queries, n_leaves, L, workers,
                               tile_edges, max_memory, width=width)
            if log:
                self._log_plan(plan, max_memory)
            txt = (f' (tiled: {plan.n_tiles:,} tiles of up to {plan.max_edges:,} edges)'
                   if plan.tile_edges is not None else '')
            return plan, txt

        if options.msa_file is None:
            if not options.queries:
                raise GPlaceExit('Give --msa_file to benchmark, or --queries N for a memory estimate only.')
            nn, ni, nl, L, sites = resources.refpkg_dimensions(options.refpkg)
            c = resources.parse_cpus(options.cpus)
            w = (resources.load_tuned_cpus() or resources.heuristic_cpus(n_avail)) if c == 'auto' else c
            plan, txt = job_plan(sites, options.queries, nl, L, w)
            self.logger.info(f'Reference package: {nn:,} nodes, {sites:,} informative sites (upper bound).')
            self.logger.info(f'Estimated peak memory for {options.queries:,} queries with {w} workers: '
                             f'{plan.memory / GiB:.1f} GiB{txt}')
            for line in slurm(plan.memory, w, None):
                self.logger.info(f'  {line}')
            return

        self._check_file(options.msa_file, 'Query alignment')
        queries = read_alignment(options.msa_file)
        q_full = options.queries or len(queries)
        rnd = random.Random(1)
        sample = queries if len(queries) <= options.max_queries else rnd.sample(queries, options.max_queries)
        cands = sorted({int(x) for x in options.cpus_list.split(',')}) if options.cpus_list \
            else resources.default_candidates(n_avail)
        cands = [c for c in cands if 1 <= c <= n_avail] or [1]
        w_pref = resources.heuristic_cpus(n_avail)

        # the full job: the informative sites of the whole alignment set its memory and tiling
        L = len(queries[0][1])
        qa = np.frombuffer(''.join(s.upper() for _, s in queries).encode(), np.uint8).reshape(len(queries), L)
        full_mask = informative(qa).any(axis=0)
        del qa
        nn, ni, nl, _, _ = resources.refpkg_dimensions(options.refpkg)
        sites_full = int((full_mask & self._ref_informative(options.refpkg, L)).sum())
        plan, _ = job_plan(sites_full, q_full, nl, L, w_pref)
        if plan.tile_edges is None and tile_edges is None and max_memory is None:
            need = resources.estimate_memory(nn, ni, sites_full, len(sample), nl, L, w_pref, width=width)
            avail = resources.available_memory()
            if avail and need > avail:
                raise GPlaceExit(f'Benchmarking this reference package untiled needs about {need / GiB:.1f} GiB, '
                                 f'more than the {avail / GiB:.1f} GiB available; give --max_memory to '
                                 f'benchmark the tiled placement.')
        self.logger.info(f'Benchmarking on {len(sample):,} of {len(queries):,} queries; worker counts '
                         f'{", ".join(map(str, cands))} ({n_avail} CPUs available).')
        if plan.tile_edges is None:
            rows, seconds = self._time_untiled(options.refpkg, sample, q_full, cands, w_pref)
        else:
            rows, seconds = self._time_tiled(options.refpkg, sample, q_full, cands, w_pref, plan.tile_edges)
        top = max(r[2] for r in rows)                  # differences under 3% are noise:
        best, _, best_rate = min((r for r in rows if r[2] >= 0.97 * top), key=lambda r: r[0])
        t_full, detail = seconds(best_rate)
        self.logger.info(f'Fastest: --cpus {best} ({best_rate:.2f} queries/s; '
                         f'{best_rate / max(r[2] for r in rows if r[0] == max(cands)):.2f}x the speed with {max(cands)}).')
        plan, txt = job_plan(sites_full, q_full, nl, L, best, log=False)
        mem = plan.memory
        self.logger.info(f'Estimated for {q_full:,} queries: peak memory {mem / GiB:.1f} GiB{txt}, '
                         f'time about {t_full / 60:.0f} min{detail}.')
        self.logger.info('Suggested SLURM request (' + ('--max_memory' if max_memory else '+10% memory')
                         + ', +50% time margin):')
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

    def _time_untiled(self, refpkg, sample, q_full, cands, w_pref):
        """Untiled benchmark: cache the sample's vectors, time the prefilter at two sample sizes
        and the placement at each worker count. Returns (rows, seconds(queries_per_second))."""
        from gplace.engine.factory import new_placer
        placer = new_placer(refpkg, sample, w_pref, discard_nonoverlapped=True)
        placer.cpus = max(cands)
        t0 = time.time(); placer.cache(); t_cache = time.time() - t0
        # prefilter time = fixed per-edge work + work per query: fit from two sweeps
        n1, n2 = len(sample), max(1, len(sample) // 4)
        t2 = placer.prefilter_seconds(n2, w_pref)
        t1 = placer.prefilter_seconds(n1, w_pref)          # leaves the full sample's scores in place
        per_q = max(0.0, (t1 - t2) / max(n1 - n2, 1))
        t_pref_full = max(t1 - per_q * n1, 0.0) + per_q * q_full
        rows = []
        for c in cands:
            t0 = time.time(); placer.place_all(workers=c); dt = time.time() - t0
            rows.append((c, dt, len(sample) / dt))
            self.logger.info(f'  {c:>4} workers: {dt:7.1f} s  ({len(sample) / dt:6.2f} queries/s)')
        placer.release()
        return rows, lambda rate: (t_cache + t_pref_full + q_full / rate, '')

    def _time_tiled(self, refpkg, sample, q_full, cands, w_pref, tile_edges):
        """Tiled benchmark: place the sample tile by tile once, timing every phase, then time the
        optimisations of one rebuilt tile at each worker count, so the tree is not rebuilt for
        every count. The rates in rows are placement-only queries/s. Returns (rows,
        seconds(queries_per_second)); see TiledPlacer.job_seconds for the time model."""
        from gplace.engine.factory import new_placer, tiled_placer
        tp = tiled_placer(new_placer(refpkg, sample, w_pref, discard_nonoverlapped=True), tile_edges)
        tp.cpus = max(cands)
        tp.record_work = True
        self.logger.info(f'Placing the sample tile by tile ({len(tp.tiles):,} tiles of up to '
                         f'{tp.max_edges:,} edges), timing each phase.')
        probe = len(sample) // 4 if len(sample) >= 8 else 0
        tp.cache()
        tp.prefilter(probe=probe)
        tp.place_all()
        ti, n_ev, rates = tp.eval_rates(cands)
        epq = tp.evals_per_query()
        batch = min(q_full, tp.query_batch)
        score_s = (tp.scoring_seconds(batch, max(1, probe), 16 * tp.workers)
                   if batch > max(1, probe) and batch > len(sample) else None)
        tp.release()
        self.logger.info(f'Timing the placement optimisations at each worker count ({n_ev:,} optimisations '
                         f'on tile {ti + 1}; {epq:.1f} per query):')
        rows = []
        for w, dt, eps in rates:
            rows.append((w, dt, eps / epq))
            self.logger.info(f'  {w:>4} workers: {dt:7.1f} s  ({eps / epq:6.2f} queries/s)')

        def seconds(rate):
            total, boundary, pref, place = tp.job_seconds(q_full, rate * epq, score_s)
            return total, (f' (boundary pass {boundary / 60:.0f} min, prefilter {pref / 60:.0f} min, '
                           f'placement {place / 60:.0f} min)')
        return rows, seconds

    @staticmethod
    def _ref_informative(refpkg, L):
        import json
        import numpy as np
        from gplace.engine.core import informative
        from gplace.engine.reference import Reference
        ref = Reference(refpkg).alignment()
        if ref is None:
            raise GPlaceExit(f'{refpkg} has no reference alignment (aln_fasta or aln_sto).')
        RA = np.frombuffer(''.join(s.upper() for _, s in ref).encode(), np.uint8).reshape(len(ref), L)
        return informative(RA).any(axis=0)

    def run_test(self, options):
        """Place the bundled test genomes and compare with GTDB-Tk's pplacer output; with --models,
        also the test sets of pplacer's other models (expected outputs from the stock pplacer)."""
        from gplace.engine import testsets
        from gplace.engine.jplace import write_jplace, compare_jplace
        from gplace.engine.core import read_fasta
        out_dir = options.out_dir or tempfile.mkdtemp(prefix='gplace_test_')
        make_sure_path_exists(os.path.join(out_dir, os.path.dirname(PATH_TEST_JPLACE)))
        tiling = self._tiling_options(options)
        t0 = time.time()
        sets = [(testsets.wag_set(), os.path.join(out_dir, PATH_TEST_JPLACE))]
        if getattr(options, 'models', False):
            work = os.path.join(out_dir, 'models')
            make_sure_path_exists(work)
            sets += [(ts, os.path.join(work, ts.name.split(' (')[0].replace(' ', '_').replace('+', '')
                                       .replace('--', '') + '.jplace')) for ts in testsets.model_sets(work)]
        failed, report = [], []
        for ts, out in sets:
            if len(sets) > 1:
                self.logger.info(f'Test set: {ts.name}')
            placer, results, _, classify = self._place_group(ts.refpkg, read_fasta(ts.queries), options.cpus,
                                                             *tiling, fix_gamma=ts.fix_gamma,
                                                             mrca_class=ts.mrca_class,
                                                             discard_nonoverlapped=ts.discard,
                                                             fit_tile_edges=ts.tile_edges)
            write_jplace(out, placer.jplace_tree(), results, ' '.join(sys.argv),
                         placer.records([r[0] for r in results]), classify)
            r = compare_jplace(out, ts.expected)
            ok = (r['identical'] == r['genomes'] and r['same_tree'] and r['same_records'] and r['same_fields']
                  and not r['only_a'] and not r['only_b'])
            report.append((ts.name, r))
            if not ok:
                failed.append(ts.name)
            if len(sets) > 1:
                self.logger.info(f'  {r["identical"]}/{r["genomes"]} identical'
                                 + ('' if ok else '  <-- DIFFERENT'))
        with open(os.path.join(out_dir, PATH_TEST_REPORT), 'w') as fh:
            fh.write('test_set\tgenomes\tbest_edge\tedge_list\tidentical\tsame_tree\n')
            for name, r in report:
                fh.write('\t'.join([name] + [str(r[k]) for k in ('genomes', 'best_edge', 'edge_list',
                                                                   'identical', 'same_tree')]) + '\n')
        r = report[0][1]
        self.logger.info(f'{r["identical"]}/{r["genomes"]} genomes: every printed value identical to '
                         f'GTDB-Tk\'s pplacer v1.1.alpha19; reference tree string identical: {r["same_tree"]}.')
        if len(sets) > 1:
            n_other = len(sets) - 1
            n_bad = sum(1 for f in failed if f != sets[0][0].name)
            self.logger.info(f'Other models: {n_other - n_bad}/{n_other} test sets identical to their expected '
                             f'output (the stock pplacer v1.1.alpha19; gPlace itself for --fix_gamma_underflow).')
        if failed:
            raise GPlaceExit(f'The test placements differ from pplacer for: {", ".join(failed)} (see {out_dir}).')
        self.logger.info(f'Test passed ({time.time() - t0:.1f} s). Output: {out_dir}')

    def check_install(self, options):
        """Check both engines (protein and DNA) are present and bit-exact."""
        from gplace.engine import core, selftest
        import platform
        from gplace.engine.build import cpu_has_avx2
        ok = True
        variant = core.engine_variant()
        self.logger.info(f'CPU: {platform.machine()}, AVX2/FMA: {"yes" if cpu_has_avx2() else "no"}; '
                         f'engine variant: {variant}'
                         + (' (set by GPLACE_ENGINE_VARIANT)' if os.getenv('GPLACE_ENGINE_VARIANT') else ''))
        for ns, what in ((20, 'protein'), (4, 'DNA')):
            path = core.engine_path(ns)
            if not os.path.exists(path):
                raise GPlaceExit(f'Engine not found at {path}. Build it with: gplace build_engine')
            L = core.lib(ns)
            self.logger.info(f'Engine ({what}): {path}')
            self.logger.info(f'        {L.xp_engine_version().decode()}')
            n, bad = selftest.check_libm(ns)
            self.logger.info(f'  libm vs GTDB-Tk pplacer binary (hypot/exp/log/pow/sin/log1p): {n - bad}/{n} identical')
            ok &= bad == 0
            if ns == 20:
                si, se, sv = selftest.check_wag_eigen()
                self.logger.info(f'  WAG model vs GTDB-Tk pplacer binary: input {"identical" if si else "DIFFERENT"}, '
                                 f'eigenvalues {"identical" if se else "DIFFERENT"}, '
                                 f'eigenvectors {"identical" if sv else "DIFFERENT"}')
                ok &= si and se and sv
            nk, fails = selftest.check_kernels(ns=ns)
            ng, gfails = selftest.check_gamma_kernels(ns=ns)
            self.logger.info(f'  kernels vs scalar references: {nk - len(fails)}/{nk} identical; gamma '
                             f'kernels vs references: {ng - len(gfails)}/{ng} identical')
            for f in (fails + gfails)[:10]:
                self.logger.error(f'  mismatch: {f}')
            ok &= not fails and not gfails
        bad_g = selftest.check_discrete_gamma()
        self.logger.info(f'Discrete gamma rates (GSL 1.16): {len(selftest.DISCRETE_GAMMA) - len(bad_g)}/'
                         f'{len(selftest.DISCRETE_GAMMA)} identical to the validated values')
        ok &= not bad_g
        if not ok:
            raise GPlaceExit('The engine is not bit-exact on this system; placements would not be '
                             'guaranteed identical to pplacer.')
        self.logger.info('Installation OK: the engines reproduce GTDB-Tk\'s pplacer arithmetic bit for bit.')

    def build_engine(self, options):
        """Compile the engines, then run the installation checks."""
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
        elif options.subparser_name == 'tog':
            from gplace.tools.tog import tog
            tog(options.placefiles, out_file=options.tog_out_file, out_dir=options.tog_out_dir,
                prefix=options.tog_prefix, pp=options.pp, xml=options.xml, node_numbers=options.node_numbers)
        elif options.subparser_name == 'taxit':
            from gplace.tools import taxit
            cmd = getattr(options, 'taxit_command', None)
            try:
                return self._taxit(taxit, cmd, options)
            except (ValueError, KeyError) as e:          # taxit leaves these uncaught: same last line
                print(f'{type(e).__name__}: {e}', file=sys.stderr)
                return 1
        else:
            self.logger.error('Unknown gPlace command: "' + options.subparser_name + '"\n')
            sys.exit(1)
        return 0

    @staticmethod
    def _taxit(taxit, cmd, options):
        if cmd == 'create':
            files = {k: getattr(options, k) for k in taxit.CREATE_FILES}
            st = taxit.create(options.package_name, options.locus, clobber=options.clobber,
                              author=options.author, description=options.description,
                              package_version=options.package_version, stats_type=options.stats_type,
                              frequency_type=options.frequency_type, reroot=options.reroot, rppr=options.rppr,
                              **files)
        elif cmd == 'check':
            st = taxit.check(options.refpkg)
        elif cmd == 'info':
            st = taxit.info(options.refpkg, seq_names=options.seq_names, tally=options.tally,
                            lengths=options.lengths)
        else:
            raise GPlaceExit('Use: gplace taxit {create,check,info} -h')
        return st
