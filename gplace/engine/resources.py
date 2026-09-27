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

"""CPU and memory planning: available CPUs, the memory model, and saved --cpus tuning."""

import json
import os
import platform
import socket
import time

from gplace.config.common import CONFIG

NS = 20


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


def estimate_memory(n_nodes, n_internal, n_sites, n_queries, n_leaves, aln_len, workers):
    """Peak memory (bytes) of one gPlace pass over a reference package.

    shared arrays : 2 x internal x sites x 20 doubles (distal/proximal) + 2 x internal x sites
                    int64 exponents + edges x queries doubles (prefilter scores)
    inputs        : alignments held as Python strings, bytes and masked arrays
    per worker    : private buffers and page tables for the shared arrays
    Calibrated on GTDB-Tk R226 class-level tree 2 (22,061 queries, 40 workers): measured cgroup
    peak 43.8 GiB (47.0 GB), model 44.7 GiB (48.0 GB); the model errs on the safe side."""
    shared = (2 * n_internal * n_sites * NS * 8 + 2 * n_internal * n_sites * 8
              + max(n_nodes - 1, 1) * n_queries * 8)
    inputs = (n_leaves + n_queries) * aln_len * 3
    per_worker = 50e6 + 0.0005 * shared
    base = 400e6
    return int(shared + inputs + workers * per_worker + base)


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
