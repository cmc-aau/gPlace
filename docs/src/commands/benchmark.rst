.. _commands/benchmark:

benchmark
=========

The ``benchmark`` command finds the fastest ``--cpus`` on the current machine, and estimates the
memory and time a job needs, for example to write a SLURM request.

Placement is limited by memory bandwidth and shared cache rather than by the number of cores, so
the fastest worker count is often below the CPU count: on a 64-vCPU AMD EPYC server it was 32-40,
and 64 workers were 1.4-1.5x slower. ``benchmark`` samples queries from ``--msa_file``, places
them with several worker counts and picks the smallest count within 3% of the fastest. The
choice is saved, so ``--cpus auto`` uses it on this machine (same host, CPU model and number of
available CPUs).

The memory estimate is for the whole job (all of ``--msa_file``, or ``--queries``); it was
calibrated on GTDB-Tk R226 class-level tree 2 (22,061 queries: estimated 44.2 GiB and about
6 min, measured 43.8 GiB and 5.7 min). With ``--refpkg`` and ``--queries`` only, gPlace prints the memory estimate
without placing anything.

With ``--max_memory`` or ``--tile_edges``, the estimate is for what :ref:`place <commands/place>`
would run with that option: untiled if it fits, otherwise the tile size gPlace picks, and the
suggested ``--mem`` is the ``--max_memory`` limit.

When the job is tiled, the benchmark is tiled too, with the same tile size, so it runs in the
memory of the job itself (on the full GTDB bacterial tree an untiled benchmark would need about
230 GiB). The sample is placed tile by tile once, timing every phase; each tile is also scored
for a quarter of the sample, which separates the prefilter's fixed cost from its cost per
query. The worker counts are then compared on one tile, rebuilt once, by replaying its recorded
optimisations, so the tree is not rebuilt for every count. Finally that tile is scored with the
sample repeated to the job's batch size (up to 20,000 queries): at that size the queries no
longer fit in the CPU caches, and scoring each query costs more than in the small sample. The
job's time is the sample's fixed costs (tile rebuilds, pool start-ups, one supernode sweep per
20,000 queries, the baseball rounds) plus its per-query costs scaled to the job, and is printed
per phase. Without ``--max_memory``/``--tile_edges``, ``benchmark`` stops with a message if the
untiled benchmark would not fit in the memory available.

Predicted and measured times of tiled jobs (40 workers, 64-vCPU AMD EPYC):

.. list-table::
   :widths: 30 16 18 16 12
   :header-rows: 1

   * - Job
     - Benchmark took
     - Predicted
     - Measured
     - Error
   * - class-level tree 2, 22,061 queries, ``--max_memory 12G``
     - 4 min
     - 11 min
     - 11 min 51 s
     - -7%
   * - full bacterial tree, 4,040 queries, ``--max_memory 60G``
     - 20 min
     - 18 min
     - 19 min 35 s
     - -8%
   * - full bacterial tree, 127,343 queries, ``--max_memory 64G``
     - 20 min
     - 2 h 34 min
     - 2 h 19 min
     - +11%

``--cpus auto`` without a saved benchmark measures the speed of a few worker counts on the first
part of the prefilter (the work is real, so nothing is repeated) and keeps the fastest; small
jobs use 5/8 of the available CPUs. The CPUs available are those the process may use, so under
SLURM ``--cpus auto`` stays within ``--cpus-per-task``.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: benchmark
   :nodefaultconst:

Files output
------------

* ``gplace.benchmark.tsv`` (with ``--out_dir``): seconds and queries per second for each worker count.
* ``~/.cache/gplace/cpus_tuning.json`` (or ``$GPLACE_CACHE_DIR``): the value used by ``--cpus auto``.

Example
-------

.. code-block:: bash

    gplace benchmark --refpkg gtdbtk.package.2.refpkg --msa_file tree_2/user_msa_file.fasta --queries 22061

.. code-block:: text

    INFO: Benchmarking on 400 of 1,270 queries; worker counts 8, 16, 24, 32, 40, 48, 64 (64 CPUs available).
    INFO:      8 workers:     9.5 s  ( 42.01 queries/s)
    INFO:     16 workers:     5.3 s  ( 75.02 queries/s)
    INFO:     24 workers:     4.0 s  ( 99.01 queries/s)
    INFO:     32 workers:     3.8 s  (105.12 queries/s)
    INFO:     40 workers:     3.9 s  (101.94 queries/s)
    INFO:     48 workers:     4.7 s  ( 85.23 queries/s)
    INFO:     64 workers:     5.3 s  ( 75.24 queries/s)
    INFO: Fastest: --cpus 32 (105.12 queries/s; 1.40x the speed with 64).
    INFO: Estimated for 22,061 queries: peak memory 44.2 GiB, time about 6 min.
    INFO: Suggested SLURM request (+10% memory, +50% time margin):
    INFO:   #SBATCH --cpus-per-task=32
    INFO:   #SBATCH --mem=49G
    INFO:   #SBATCH --time=00:19:31

Memory only:

.. code-block:: bash

    gplace benchmark --refpkg gtdbtk.package.2.refpkg --queries 50000 --cpus 32

Tile size for the full GTDB bacterial tree in a 64 GB job:

.. code-block:: bash

    gplace benchmark --refpkg gtdb_r226_bac120.refpkg --queries 127343 --cpus 40 --max_memory 64G

.. code-block:: text

    INFO: --max_memory 64.0 GiB: the untiled placement would need about 490.1 GiB; using --tile_edges 56,726 (estimated peak 60.3 GiB).
    INFO: Reference package: 273,290 nodes, 5,036 informative sites (upper bound).
    INFO: Estimated peak memory for 127,343 queries with 40 workers: 60.3 GiB (tiled: 6 tiles of up to 55,595 edges)
    INFO:   #SBATCH --cpus-per-task=40
    INFO:   #SBATCH --mem=64G
