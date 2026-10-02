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

"""CPU and memory planning: available CPUs, the memory models, and saved --cpus tuning."""

import json
import os
import platform
import re
import socket
import time

from gplace.config.common import CONFIG

NS = 20
GiB = 2 ** 30

# tiled placement (gplace.engine.tiled)
QUERY_BATCH = 20000         # queries scored per prefilter sweep of a tile
N_CAND = 160                # ranked edges kept per query (pplacer's search uses at most ~41)
MERGE_CELLS = 4_000_000     # candidate-merge chunk size (rows x columns)


def available_cpus():
    """CPUs this process may use (respects SLURM/cgroup CPU affinity)."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def default_candidates(n):
    """Worker counts tried by the benchmark / auto-calibration: n/8 ... n."""
    fr = (0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0)
    return sorted({max(1, int(round(n * f))) for f in fr})


def heuristic_cpus(n):
    """Fallback when nothing has been measured: placement is memory-bandwidth bound and on the
    reference 64-vCPU AMD EPYC server the optimum was 40 (5/8 of the vCPUs)."""
    return max(1, int(round(0.625 * n)))


def parse_cpus(value):
    """'auto' or a positive integer."""
    if isinstance(value, int):
        return value
    if str(value).lower() == 'auto':
        return 'auto'
    v = int(value)
    if v < 1:
        raise ValueError('--cpus must be >= 1 or "auto"')
    return v


def estimate_memory(n_nodes, n_internal, n_sites, n_queries, n_leaves, aln_len, workers, width=NS):
    """Peak memory (bytes) of one gPlace pass over a reference package.

    shared arrays : 2 x internal x sites x width doubles (distal/proximal; width = 20 for protein CAT,
                    4 for DNA CAT, rates x states for gamma) + 2 x internal x sites
                    int64 exponents + edges x queries doubles (prefilter scores)
    inputs        : alignments held as Python strings, bytes and masked arrays
    per worker    : private buffers and page tables for the shared arrays
    Calibrated on GTDB-Tk R226 class-level tree 2 (22,061 queries, 40 workers): measured cgroup
    peak 43.8 GiB (47.0 GB), model 44.7 GiB (48.0 GB); the model errs on the safe side."""
    shared = (2 * n_internal * n_sites * width * 8 + 2 * n_internal * n_sites * 8
              + max(n_nodes - 1, 1) * n_queries * 8)
    inputs = (n_leaves + n_queries) * aln_len * 3
    per_worker = 50e6 + 0.0005 * shared
    base = 400e6
    return int(shared + inputs + workers * per_worker + base)


def estimate_memory_tiled(n_slots, max_edges, n_tiles, n_sites, n_queries, n_leaves, aln_len, workers,
                          query_batch=QUERY_BATCH, n_cand=N_CAND, width=NS):
    """Peak memory (bytes) of one tiled pass (gplace.engine.tiled).

    shared arrays : one tile's vectors (2 x slots x sites x width doubles + int64 exponents; slots =
                    the tile's internal nodes and the roots of the tiles below it) + the prefilter
                    scores of one query batch against the largest tile
    boundaries    : distal and proximal vector of every tile root
    per query     : ranked candidate edges and the stored branch-length optimisations
    inputs, per worker, base : as in estimate_memory
    Calibrated on the GTDB R226 full bacterial tree (273,290 nodes, 5,036 sites; 24 tiles of
    <= 13,889 edges and 6,944 slots; 40 workers): 4,040 queries measured 14.9 GiB, model 16.0 GiB;
    127,343 queries measured 19.9 GiB (cgroup peak), model 20.5 GiB. Same tree with larger tiles:
    6 tiles of <= 55,595 edges, 127,343 queries measured 59.2 GiB, model 60.3 GiB; 5 tiles of
    <= 61,933 edges, 4,040 queries measured 54.3 GiB, model 56.1 GiB."""
    shared = (2 * n_slots * n_sites * width * 8 + 2 * n_slots * n_sites * 8
              + min(n_queries, query_batch) * max_edges * 8)
    boundaries = 2 * n_tiles * n_sites * (width + 1) * 8
    per_query = n_cand * 8 * 3 + 5000
    merge = MERGE_CELLS * 40
    inputs = (n_leaves + n_queries) * aln_len * 3
    per_worker = 50e6 + 0.0005 * shared
    base = 400e6
    return int(shared + boundaries + n_queries * per_query + merge + inputs + workers * per_worker + base)


def parse_memory(value):
    """'64G', '64GB', '64GiB', '500M', '1.5T' or a plain number (GiB) -> bytes. Units are binary
    (as in SLURM's --mem): K = 1024, M = 1024^2, G = 1024^3, T = 1024^4."""
    m = re.fullmatch(r'\s*([0-9]*\.?[0-9]+)\s*([KMGT]?)(I?B)?\s*', str(value), re.IGNORECASE)
    if not m:
        raise ValueError(f'invalid memory size "{value}" (examples: 64G, 500M, 1.5T)')
    unit = m.group(2).upper() or ('' if m.group(3) else 'G')
    size = float(m.group(1)) * 1024 ** ('KMGT'.index(unit) + 1 if unit else 0)
    if size < 1:
        raise ValueError(f'memory size must be positive: "{value}"')
    return int(size)


def available_memory():
    """Memory this process may use (bytes): the smallest cgroup limit on its path (SLURM jobs,
    systemd scopes, containers) or the machine's total memory. None if unknown."""
    limits = []
    try:
        with open('/proc/meminfo') as fh:
            for line in fh:
                if line.startswith('MemTotal:'):
                    limits.append(int(line.split()[1]) * 1024)
                    break
    except (OSError, ValueError):
        pass
    try:
        with open('/proc/self/cgroup') as fh:
            entries = [line.rstrip('\n').split(':', 2) for line in fh]
    except OSError:
        entries = []
    for _, controllers, path in (e for e in entries if len(e) == 3):
        if controllers == '':                                   # cgroup v2
            root, files = '/sys/fs/cgroup', ('memory.max',)
        elif 'memory' in controllers.split(','):                # cgroup v1
            root, files = '/sys/fs/cgroup/memory', ('memory.limit_in_bytes',)
        else:
            continue
        parts = [p for p in path.split('/') if p]
        for i in range(len(parts), -1, -1):
            for name in files:
                try:
                    raw = open(os.path.join(root, *parts[:i], name)).read().strip()
                except OSError:
                    continue
                if raw.isdigit() and int(raw) < 2 ** 60:
                    limits.append(int(raw))
    return min(limits) if limits else None


def refpkg_width(refpkg_dir):
    """Doubles per site of the likelihood vectors of a refpkg's model (models.width_of)."""
    from gplace.engine import models
    return models.width_of(models.read_phylo_model(refpkg_dir))


def refpkg_dimensions(refpkg_dir):
    """(n_nodes, n_internal, n_leaves, aln_len, ref_informative_sites) without building anything."""
    import numpy as np
    from gplace.engine.core import read_fasta, informative
    from gplace.engine.newick import parse_newick
    files = json.load(open(os.path.join(refpkg_dir, 'CONTENTS.json')))['files']
    t = parse_newick(open(os.path.join(refpkg_dir, files['tree'])).read())
    ref = read_fasta(os.path.join(refpkg_dir, files['aln_fasta']))
    L = len(ref[0][1])
    RA = np.frombuffer(''.join(s for _, s in ref).encode(), np.uint8).reshape(len(ref), L)
    sites = int(informative(RA).any(axis=0).sum())
    n_leaves = len(t.leaves())
    return t.n_nodes, t.n_nodes - n_leaves, n_leaves, L, sites


# ------------------------------------------------------------------ saved tuning --
def machine_key():
    model = platform.processor() or ''
    try:
        with open('/proc/cpuinfo') as fh:
            for line in fh:
                if line.startswith('model name'):
                    model = line.split(':', 1)[1].strip()
                    break
    except OSError:
        pass
    return f'{socket.gethostname()}|{model}|{available_cpus()} cpus'


def _tuning_path():
    return os.path.join(CONFIG.CACHE_DIR, 'cpus_tuning.json')


def load_tuned_cpus():
    """--cpus chosen by `gplace benchmark` on this machine (same host, CPU model and CPU count)."""
    try:
        data = json.load(open(_tuning_path()))
        return int(data[machine_key()]['cpus'])
    except (OSError, ValueError, KeyError, TypeError):
        return None


def save_tuned_cpus(cpus, details):
    path = _tuning_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    try:
        data = json.load(open(path))
    except (OSError, ValueError):
        data = {}
    data[machine_key()] = dict(details, cpus=int(cpus), date=time.strftime('%Y-%m-%d %H:%M:%S'))
    with open(path, 'w') as fh:
        json.dump(data, fh, indent=2)
    return path


def slurm_time(seconds):
    s = int(seconds)
    return f'{s // 86400}-{(s % 86400) // 3600:02d}:{(s % 3600) // 60:02d}:{s % 60:02d}' if s >= 86400 \
        else f'{s // 3600:02d}:{(s % 3600) // 60:02d}:{s % 60:02d}'
