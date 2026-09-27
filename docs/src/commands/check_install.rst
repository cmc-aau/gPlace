.. _commands/check_install:

check_install
=============

The ``check_install`` command checks that the compiled engine reproduces pplacer's arithmetic bit for bit on this machine: the CPU features, the libm functions against values taken from GTDB-Tk's pplacer binary, the WAG eigendecomposition against the one inside that binary, and every SIMD kernel against its scalar reference.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: check_install
   :nodefaultconst:

Example
-------

.. code-block:: bash

    gplace check_install

.. code-block:: text

    [2026-09-26 20:19:08] INFO: CPU AVX2/FMA support: yes
    [2026-09-26 20:19:08] INFO: Engine: gplace/engine/lib/libgplace_core.so
    [2026-09-26 20:19:08] INFO:         gplace-core 1: pplacer v1.1.alpha19 exact (GSL 1.16, glibc-2.19 hypot, CORE-MATH exp/log, AVX2)
    [2026-09-26 20:19:08] INFO: libm vs GTDB-Tk pplacer binary (hypot/exp/log): 1500/1500 identical
    [2026-09-26 20:19:08] INFO: WAG model vs GTDB-Tk pplacer binary: input identical, eigenvalues identical, eigenvectors identical
    [2026-09-26 20:19:09] INFO: SIMD kernels vs scalar references: 254/254 identical
    [2026-09-26 20:19:09] INFO: Installation OK: the engine reproduces GTDB-Tk's pplacer arithmetic bit for bit.
