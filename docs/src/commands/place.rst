.. _commands/place:

place
=====

The ``place`` command places the aligned queries of one alignment file on one pplacer reference package, or on reference files given as with pplacer's ``-t``/``-r``/``-s``, and writes a jplace file. It computes exactly what GTDB-Tk's ``pplacer -m WAG -j <cpus> -c <refpkg> -o <out.json> <msa>`` computes: same placements, same printed numbers, same tree string. The model comes from the reference package (WAG, LG, JTT or GTR; CAT or gamma; see :ref:`models`); ``--fix_gamma_underflow`` avoids pplacer's gamma underflow on large trees (:ref:`models/gamma_underflow`).

pplacer's other placement options are available under gPlace's names, and under pplacer's own
spelling, so a pplacer command line can be reused (``-c <refpkg>`` becomes ``--refpkg``, the query
file ``--msa_file``, ``-o`` ``--out_file``). With the same options gPlace gives the same output as
pplacer (see :ref:`validation/options`).

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: place
   :nodefaultconst:

Taxonomy, duplicates and non-overlapping queries
------------------------------------------------

gPlace follows pplacer here too:

* **Classification.** If the reference package lists a ``taxonomy`` (pplacer's tax table) and a
  ``seq_info`` file, e.g. one built with ``taxit create --taxonomy ... --seq-info ...``, every
  placement gets a ``classification`` field: the taxon of its edge by *edge painting* (the most
  specific rank at which exactly one taxon is cut by the edge), or with ``--mrca_class`` the
  taxonomic MRCA of the reference sequences below the edge. GTDB-Tk's reference packages have no
  taxonomy, so their output has no such field.
* **Duplicates.** Queries whose sequences are identical after pre-masking are placed once and
  written as one record listing all their names.
* **Non-overlapping queries.** A query with no residue in any column where a reference sequence
  has one stops the run (``Sequence ... doesn't overlap any reference sequences.``), unless
  ``--discard_nonoverlapped`` drops it.

Query files
-----------

As in pplacer, the query alignment is read as Stockholm if its name ends in ``.sto`` or ``.sth``,
otherwise as FASTA; a file ending in ``.gz`` is read gunzipped (``queries.sto.gz``). Stockholm
files are read as pplacer reads them: one alignment, ``#=`` markup ignored, interleaved blocks
joined, ``.`` read as a gap, lowercase letters kept as residues (insert columns are not removed).
pplacer orders the sequences of a Stockholm file by name, descending, and so does gPlace; the
order matters when the file holds the reference sequences and the model uses empirical
frequencies.

pplacer's settings
------------------

gPlace uses pplacer's defaults, which GTDB-Tk does not change. Each can be set as in pplacer:

======================== ======================== ========= ==================================================
gPlace                   pplacer                  default   meaning
======================== ======================== ========= ==================================================
``--ml_tolerance``       ``--ml-tolerance``       0.01      tolerance of the first-stage branch length
                                                            optimisation (kept placements: 1e-5)
``--start_pend``         ``--start-pend``         0.1       starting pendant branch length
``--max_pend``           ``--max-pend``           2.0       maximum pendant branch length
``--max_strikes``        ``--max-strikes``        6         edges tried stop after this many strikes; ``0``
                                                            tries every edge of the tree (slow on large trees)
``--strike_box``         ``--strike-box``         3.0       a strike: an edge this many log-likelihood units
                                                            below the best so far
``--max_pitches``        ``--max-pitches``        40        edges tried stop after this many (plus one)
``--keep_at_most``       ``--keep-at-most``       7         placements kept per query
``--keep_factor``        ``--keep-factor``        0.01      placements kept if their likelihood is at least
                                                            this fraction of the best one's
``--no_pre_mask``        ``--no-pre-mask``        off       see below
``--always_refine``      ``--always-refine``      off       see `Model refinement`_
======================== ======================== ========= ==================================================

``--start_pend`` must be positive and below ``--max_pend`` (pplacer's check). A query left with
no placement (``--keep_at_most 0``, or ``--keep_factor`` above 1) is omitted from the jplace, as
pplacer does, with a warning.

**Pre-masking.** By default pplacer, and gPlace, keep only the alignment columns where some
reference sequence and some query have a residue, stop at queries without a residue in any
column where a reference sequence has one, and place identical masked sequences once. With
``--no_pre_mask`` every column is kept, no query is checked or dropped (``--discard_nonoverlapped``
has no effect), only queries with identical whole sequences share a record, and a query without
any residue stops the run (``sequence '...' has no informative sites.``). Placements change for
queries with residues in columns where no reference sequence has one: those residues then count
towards the query's likelihood.

Reference files without a reference package
-------------------------------------------

Instead of ``--refpkg``, the reference can be given as files, as with pplacer's ``-t``, ``-r``
and ``-s``; with ``--refpkg`` they replace the package's files:

* ``--tree`` (``-t``): the reference tree, Newick.
* ``--ref_aln`` (``-r``): the reference alignment, FASTA. Without it, the query alignment must
  contain every reference sequence (see `Model refinement`_); otherwise every reference
  sequence must be absent from it.
* ``--stats`` (``-s``): the model, from a RAxML 7 info file (including the format RAxML writes with
  ``-f e``), a PhyML 3.0 stats file, or a JSON model file (``datatype``, ``subs_model``,
  ``ras_model`` ``gamma`` or ``null``, ``gamma``, ``subs_rates``; protein: WAG only). The model is
  always gamma; files are read exactly as pplacer v1.1.alpha19 reads them, including its limits:
  RAxML 8 info files and PhyML 3.3 stats files do not parse (``taxit create`` reads them; use a
  reference package). A reference package's ``phylo_model`` takes precedence over ``--stats``.
* ``--ref_dir`` (``-d``): a directory the three file names are relative to.

What a statistics file does not set comes from pplacer's model settings: ``--model_freqs``
(protein: the model's amino acid frequencies instead of the reference alignment's, which are
pplacer's default for statistics files), ``--gamma_cats`` and ``--gamma_alpha`` (when the file has
no gamma model: PhyML without gamma, JSON with ``ras_model`` ``null``; defaults 1 and 1.0).
``--model`` (``-m``) is accepted but, as in pplacer v1.1.alpha19, has no effect: the file names
the model. A reference package without ``phylo_model`` but with ``tree_stats`` is read the same
way (pplacer's deprecated path).

.. code-block:: bash

    gplace place --tree hiv.phy_phyml_tree.txt --ref_aln hiv.fasta --stats hiv.phy_phyml_stats.txt \
        --msa_file queries.fasta --out_file queries.jplace

Model refinement
----------------

When the query alignment contains the reference sequences (all of them; some but not all stops
the run), pplacer uses them as the reference alignment and refines the model before placing;
``--always_refine`` refines it in any case. For CAT models (FastTree reference packages, such as
GTDB-Tk's) the refinement chooses a new rate category for every site: the category with the
highest likelihood of the reference alignment on the reference tree, penalised as in pplacer
(log-likelihood + 2 log(rate) - 3 rate). Placements then differ from those with the package's
categories. For gamma models refinement changes nothing (as in pplacer).

Memory limits and large trees
-----------------------------

By default gPlace holds the likelihood vectors of the whole reference tree in memory, which is
the fastest way to place. For trees whose vectors do not fit, such as the full GTDB bacterial
tree (about 490 GiB untiled for 127,343 genomes), gPlace can place **tile by tile**: the tree is
cut into connected tiles and only one tile is held in memory at a time. The results are
identical to untiled placement (see :ref:`how_it_works`); tiled runs take about 1.7-2.5x longer.

Give at most one of:

* ``--max_memory 64G``: the memory gPlace may use (``K``, ``M``, ``G``, ``T`` are binary units,
  as in SLURM; a plain number is GiB). gPlace estimates the peak memory from the reference tree
  and the number of queries, keeping 5% of the limit as headroom. If untiled placement fits,
  it places without tiles; otherwise it uses the largest tiles that fit, because larger tiles
  are faster. If even tiles of 1,000 edges do not fit, it stops and prints the limit needed.
* ``--tile_edges 20000``: tiles of at most this many edges. Memory then depends on the tile
  size rather than the tree size: about 8 GiB of likelihood vectors per 10,000 edges (5,036
  sites), plus the alignments and the queries' working data.

``gplace benchmark --refpkg <refpkg> --queries <N> --max_memory 64G`` prints the tile size and
peak memory gPlace would use, without placing anything. Without either option, gPlace warns
when the estimated peak is more than the memory available to it (the machine's memory, or a
SLURM or cgroup limit).

Files output
------------

* :ref:`jplace file <files/jplace>` (``--out_file``)
* :ref:`gplace.log <files/gplace.log>` (next to the jplace file)

Example
-------

Input
^^^^^

.. code-block:: bash

    gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4

For a GTDB-Tk class-level tree, pass e.g. ``--refpkg release226/split/class_level/pplacer/gtdbtk.package.2.refpkg``.
On the full GTDB bacterial tree within a 64 GB job:

.. code-block:: bash

    gplace place --msa_file gtdbtk.bac120.user_msa.fasta --refpkg release226/pplacer/gtdb_r226_bac120.refpkg \
        --out_file gplace.jplace --cpus 40 --max_memory 64G

Output
^^^^^^

.. code-block:: text

    [2026-10-01 19:23:40] INFO: gPlace v0.1.0
    [2026-10-01 19:23:40] INFO: gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4
    [2026-10-01 19:23:42] INFO: Read 20 sequences from user_msa_file.fasta
    [2026-10-01 19:23:44] INFO: Model: WAG + CAT (20 rates), model frequencies (protein)
    [2026-10-01 19:23:44] INFO: Placing 20 queries on mini.refpkg (493 nodes, 4,891 sites, 4 CPUs; estimated peak memory 0.9 GiB).
    [2026-10-01 19:23:47] INFO: Cached likelihood vectors for 246 internal nodes (0.4 GB).
    [2026-10-01 19:23:48] INFO: Scored 20 queries against 492 edges.
    [2026-10-01 19:23:51] INFO: Wrote 20 placements to gplace.jplace (8.0 s).
