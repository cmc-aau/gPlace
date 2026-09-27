.. _validation:

Validation
==========

gPlace was validated on the 127,343 GlobDB genomes that GTDB-Tk (release R226) assigned to the
8 bacterial class-level reference packages, by comparing gPlace's jplace files with the jplace
files written by GTDB-Tk's pplacer v1.1.alpha19.

.. list-table::
   :widths: 10 12 12 14 14 14 12 12
   :header-rows: 1

   * - Tree
     - Genomes
     - Nodes
     - Best edge identical
     - Placement list identical
     - All printed values identical
     - gPlace time
     - Earlier exact build
   * - 1
     - 14,567
     - 40,100
     - 100%
     - 100%
     - 100%
     - 3 min 30 s
     - 12 min 54 s
   * - 2
     - 22,061
     - 43,130
     - 100%
     - 100%
     - 100%
     - 5 min 41 s
     - 18 min 49 s
   * - 3
     - 12,020
     - 33,840
     - 100%
     - 100%
     - 100%
     - 2 min 37 s
     - 10 min 41 s
   * - 4
     - 16,809
     - 34,020
     - 100%
     - 100%
     - 100%
     - 3 min 24 s
     - 12 min 54 s
   * - 5
     - 11,424
     - 34,722
     - 100%
     - 100%
     - 100%
     - 2 min 35 s
     - 10 min 05 s
   * - 6
     - 25,016
     - 39,238
     - 100%
     - 100%
     - 100%
     - 5 min 19 s
     - 19 min 39 s
   * - 7
     - 20,018
     - 37,070
     - 100%
     - 100%
     - 100%
     - 4 min 17 s
     - 17 min 15 s
   * - 8
     - 5,428
     - 13,714
     - 100%
     - 100%
     - 100%
     - 57 s
     - 3 min 52 s
   * - **All**
     - **127,343**
     -
     - **127,343**
     - **127,343**
     - **127,343**
     - **28 min 25 s**
     - 1 h 46 min

* *All printed values identical*: every number of every placement (edge, likelihood,
  like-weight ratio, distal and pendant lengths) is textually identical to GTDB-Tk's file, and the
  112 reference tree strings are byte-identical.
* One run of ``gplace place_batch`` over the 112 GTDB-Tk alignment files, 40 worker processes, on
  a 64-vCPU AMD EPYC server, under a hard 64 GB memory limit (``systemd-run MemoryMax=64G``,
  swap disabled). Peak memory: 43.8 GiB = 47.0 GB (tree 2, whole process group). pplacer itself needs ~104 GB to cache tree 2.
* *Earlier exact build*: the first bit-exact version (scalar C kernels, 64 threads), before the
  SIMD kernels, the process pool and the memory-lean three-taxon kernel.

How exactness was established
-----------------------------

Where results differed during development, GTDB-Tk's pplacer binary (statically linked, with
debug symbols) was traced under gdb and compared with gPlace evaluation by evaluation. This
found the two library differences (GSL 1.16 and glibc 2.19's ``hypot``) and the NaN-retry
behaviour described in :ref:`how_it_works`. ``gplace check_install`` and ``gplace test`` repeat
the key comparisons on every installation.
