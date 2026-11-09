.. _advanced/pplacer_options:

pplacer's options
=================

gPlace takes pplacer's command line (``gplace place -c <refpkg> -o <out.jplace> <queries>``, see
:ref:`commands/place`) and uses pplacer's default settings, so a plain ``pplacer -c <refpkg>`` run
needs nothing on this page. pplacer's other options are available too, under gPlace's name and under
pplacer's own spelling, and give the same output as pplacer (see :ref:`validation/options`).
They work the same way with :ref:`commands/place_batch`, except reference files without a
reference package and query files holding the reference sequences, which need
:ref:`commands/place`.


Search settings
---------------

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
``--no_pre_mask``        ``--no-pre-mask``        off       see `Pre-masking`_
``--always_refine``      ``--always-refine``      off       see `Model refinement`_
======================== ======================== ========= ==================================================

``--start_pend`` must be positive and below ``--max_pend`` (pplacer's check). A query left with
no placement (``--keep_at_most 0``, or ``--keep_factor`` above 1) is omitted from the jplace, as
pplacer does, with a warning.


Posterior probabilities
-----------------------

======================== ======================== ========= ==================================================
gPlace                   pplacer                  default   meaning
======================== ======================== ========= ==================================================
``-p``, ``--pp``         ``-p``                   off       calculate posterior probabilities
``--pp_rel_err``         ``--pp-rel-err``         0.01      relative error of the integration
``--unif_prior``         ``--unif-prior``         off       uniform prior on the pendant branch length
``--inform_prior``       ``--inform-prior``       off       exponential prior with a mean for each edge
``--prior_lower``        ``--prior-lower``        0         added to each edge's mean (``--inform_prior``)
======================== ======================== ========= ==================================================

With ``-p``, gPlace computes what pplacer's ``-p`` computes for every kept placement: its
*marginal likelihood*, the likelihood integrated over the attachment point along the edge and over
the pendant branch length, weighted by a prior on the pendant branch length; and its *posterior
probability*, the marginal likelihoods normalised over the query's placements. Each jplace row
gets the fields ``marginal_like`` and ``post_prob``. The placements stay sorted by like-weight
ratio, as in pplacer, so the highest posterior probability is not always first.

As in pplacer, ``-p`` also sets ``--keep_at_most 20``, ``--keep_factor 0.001`` and
``--max_strikes 20`` where it appears on the command line: these options override it when they
come after ``-p``, and are overridden when they come before. The wider search can change the
best placement, so ``-p`` results differ from likelihood placement on the same data.

The prior is exponential, with the average branch length of the reference tree as its mean,
unless ``--unif_prior`` (uniform) or ``--inform_prior`` (for each edge, its length plus the average
distance from its lower end to the leaves below it, plus ``--prior_lower``). When an integration
does not reach ``--pp_rel_err`` even after narrowing its range ten times, pplacer warns and uses
the placement's likelihood as its marginal likelihood; gPlace does the same and says so in its
log.

Posterior probabilities cost much more time than likelihood placement (each kept placement
needs hundreds to thousands of likelihood evaluations) but no extra memory, and work with tiled
placement.


Pre-masking
-----------

By default pplacer, and gPlace, keep only the alignment columns where some reference sequence and
some query have a residue, stop at queries without a residue in any column where a reference
sequence has one, and place identical masked sequences once. With ``--no_pre_mask`` every column
is kept, no query is checked or dropped (``--discard_nonoverlapped`` has no effect), only queries
with identical whole sequences share a record, and a query without any residue stops the run
(``sequence '...' has no informative sites.``). Placements change for queries with residues in
columns where no reference sequence has one: those residues then count towards the query's
likelihood.


Taxonomy, duplicates and non-overlapping queries
------------------------------------------------

* **Classification.** If the reference package lists a ``taxonomy`` (pplacer's tax table) and a
  ``seq_info`` file, e.g. one built with ``taxit create --taxonomy ... --seq-info ...``, every
  placement gets a ``classification`` field: the taxon of its edge by *edge painting* (the most
  specific rank at which exactly one taxon is cut by the edge), or with ``--mrca_class`` the
  taxonomic MRCA of the reference sequences below the edge. Reference packages without a taxonomy,
  such as the GTDB ones, give no such field.
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
  reference package). An IQ-TREE ``.iqtree`` file or a RAxML-NG log, which pplacer cannot read, gives the
  model ``taxit create`` writes from it (see :ref:`advanced/own_trees`). A reference
  package's ``phylo_model`` takes precedence over ``--stats``.
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
the GTDB ones) the refinement chooses a new rate category for every site: the category with the
highest likelihood of the reference alignment on the reference tree, penalised as in pplacer
(log-likelihood + 2 log(rate) - 3 rate). Placements then differ from those with the package's
categories. For gamma models refinement changes nothing (as in pplacer).


Options gPlace does not implement
---------------------------------

These pplacer options write extra files or change what pplacer computes in ways gPlace does not
reproduce. Given to ``gplace place``, they stop the run with a message (with their default value,
which changes nothing, they are accepted):

* ``--map-identity``, ``--map-mrca``: MAP sequences and percent identity (``--map-mrca-min`` only
  applies to them);
* ``--groups`` (above 1): placing the queries in groups;
* ``--fig-cutoff`` (above 0), ``--fig-eval-all``, ``--fig-eval-discrepancy-tree``, ``--fig-tree``:
  fig-based search;
* ``--fantasy`` (above 0): fantasy baseball (``--fantasy-frac`` only applies to it);
* ``--write-masked``, ``--write-pre-masked``: masked alignments;
* ``--pretend``, ``--check-like``: input checks and reference likelihood reports.

``--mmap-file``, ``--verbosity`` and ``--timing`` are accepted and change nothing: gPlace needs no
scratch file, and logs its progress and timings in ``gplace.log``.
