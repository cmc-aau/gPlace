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
