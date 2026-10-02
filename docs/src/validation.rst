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

Tiled placement
---------------

Tiled runs (``--tile_edges`` / ``--max_memory``) on the same machine, 40 worker processes, each
under a hard 64 GB limit. Peak memory is the whole process group (cgroup peak).

.. list-table::
   :widths: 16 10 18 22 14 12 12
   :header-rows: 1

   * - Tree
     - Genomes
     - Option (tiles)
     - Compared with
     - Result
     - Peak memory (estimate)
     - Time (untiled)
   * - class-level 8
     - 5,428
     - ``--tile_edges 3000`` (5 tiles)
     - GTDB-Tk jplace files
     - all printed values identical
     - 4.0 GiB (5.1)
     - 2 min 47 s (1 min 48 s)
   * - class-level 2
     - 22,061
     - ``--max_memory 12G`` (6 tiles)
     - GTDB-Tk jplace files
     - all printed values identical
     - 10.0 GiB (11.1)
     - 11 min 51 s (5 min 41 s)
   * - full bacterial
     - 4,040
     - ``--max_memory 60G`` (5 tiles)
     - GTDB-Tk full-tree classification tree
     - same edge, pendant and distal lengths
     - 54.3 GiB (56.1)
     - 19 min 35 s (needs ~232 GiB)
   * - full bacterial
     - 127,343
     - ``--max_memory 64G`` (6 tiles)
     - GTDB-Tk full-tree classification tree (3,459 genomes in common)
     - same edge, pendant and distal lengths
     - 59.2 GiB (60.3)
     - 2 h 19 min (needs ~490 GiB)

* The full R226 bacterial tree (273,290 nodes) is not split by GTDB-Tk; its reference is the
  classification tree GTDB-Tk writes for a full-tree run of 4,040 genomes, in which each genome is
  grafted at its best placement. It records the edge and the pendant and distal lengths with 6
  significant digits, but not likelihoods or secondary placements. Every genome in it is on the
  same edge, with the same pendant length as printed and distal lengths within the printing
  precision (largest differences 4.6e-6 and 8.4e-7).
* The 127,343-genome run placed the 112 GTDB-Tk class-level alignment files on the full tree in
  one ``place_batch``. Its jplace files are also identical, in every printed value, to those of
  the independent tiling prototype that preceded this implementation (24 tiles of up to 13,889
  edges, 19.9 GiB peak, 2 h 38 min).
* *Time (untiled)*: the untiled run of the same job, where it fits in memory (tree 8 rerun
  alone; tree 2 from the table above).

Other models
------------

Every model pplacer v1.1.alpha19 supports (:ref:`models`) was compared with the stock pplacer
binary on the same inputs, untiled and tiled (``gplace place``; every placement record and the tree
string textually identical, including pplacer's ``-inf``/``-nan`` values):

.. list-table::
   :header-rows: 1

   * - Model
     - Reference package
     - Genomes
     - Result (untiled and tiled)
   * - LG, JTT; WAG/LG/JTT + empirical frequencies; WAG, LG, JTT + gamma (alpha 0.05-1.3, 4 and 8 rates)
     - bundled test package (``gplace test --models``)
     - 20 each
     - identical
   * - GTR, WAG, JTT + gamma (PhyML)
     - pplacer's own test data
     - 60 each
     - identical
   * - GTR + CAT (FastTree + taxit), GTR + gamma (raxml-ng)
     - simulated 16S-like DNA, 1,200 and 5,000 references
     - 300 / 1,000
     - identical
   * - WAG (GTDB-Tk), LG, JTT, WAG + empirical, LG + empirical
     - GTDB-Tk R226 class-level package 8
     - 343 each
     - identical
   * - WAG + gamma (alpha 0.8)
     - GTDB-Tk R226 class-level package 8
     - 343
     - identical (342 genomes -inf in both, see :ref:`models/gamma_underflow`)

The gamma rates use GSL's gamma functions, which call ``exp``, ``log``, ``pow`` and ``sin``. Their
values in GTDB-Tk's pplacer binary (glibc 2.19) were compared under gdb with CORE-MATH's
correctly-rounded functions on 8,088 inputs (the arguments GSL passes for discrete-gamma rates,
plus random ones): identical in every case, while the system glibc 2.39 differed in 7. gPlace binds
GSL's calls to CORE-MATH, and ``gplace check_install`` checks these values.

How exactness was established
-----------------------------

Where results differed during development, GTDB-Tk's pplacer binary (statically linked, with
debug symbols) was traced under gdb and compared with gPlace evaluation by evaluation. This
found the two library differences (GSL 1.16 and glibc 2.19's ``hypot``) and the NaN-retry
behaviour described in :ref:`how_it_works`. ``gplace check_install`` and ``gplace test`` repeat
the key comparisons on every installation.
