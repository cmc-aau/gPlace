.. _advanced/own_trees:

Using your own tree
===================

gPlace places on any reference tree pplacer can use, not only GTDB-Tk's. The simplest way is to
build a reference package with :ref:`gplace taxit create <commands/taxit>` (the same package
taxtastic's ``taxit`` builds, usable by pplacer too) from your alignment, your tree and the log or
statistics file of the program that made the tree (FastTree, RAxML, RAxML-NG, IQ-TREE or PhyML):

.. code-block:: bash

    gplace taxit create -l bac -P bac.refpkg --aln-fasta ref.fasta --tree-file tree.nwk --tree-stats fasttree.log
    gplace place --refpkg bac.refpkg --msa_file queries.fasta --out_file queries.jplace

The reference files can also be given directly, as with pplacer's ``-t``/``-r``/``-s``: see
:ref:`advanced/pplacer_options`. The models pplacer, and so gPlace, can place with are listed in
:ref:`models`.


Trees fitted with IQ-TREE or RAxML-NG
-------------------------------------

A reference tree from IQ-TREE or RAxML-NG can be used in two ways, with the same placements:

.. code-block:: bash

    # a reference package, usable by pplacer too
    gplace taxit create -l bac -P bac.refpkg --aln-fasta ref.fasta --tree-file run.treefile --tree-stats run.iqtree
    gplace place --refpkg bac.refpkg --msa_file queries.fasta --out_file queries.jplace

    # the files directly (RAxML-NG: --tree run.raxml.bestTree --stats run.raxml.log)
    gplace place --tree run.treefile --ref_aln ref.fasta --stats run.iqtree --msa_file queries.fasta \
        --out_file queries.jplace

``taxit create`` (gPlace's, or taxtastic's: the same package) keeps from the fit the substitution
matrix, the gamma model (categories and alpha, as printed in the file), GTR's rates and whether the
frequencies are counted from the alignment. pplacer places with that, and so does gPlace, with
identical results (validated on 13 IQ-TREE 3 and RAxML-NG 2 fits, protein and DNA). The rest of the
fitted model is not used; gPlace says so in a warning that names the fit and what pplacer uses
instead:

* invariant sites (``+I``) and free rates (``+R``);
* ML-estimated frequencies (``+FO``): counted from the reference alignment for GTR and when taxit marks
  them empirical (RAxML-NG), else the matrix's own;
* RAxML-NG's counted frequencies (``+FC``): taxit records them as model frequencies, so a protein
  matrix uses its own;
* median gamma rates (pplacer uses mean rates) and partitions.

pplacer cannot use a fit without gamma rate variation (``LG``, ``LG+F``, ``+I+G4``, which taxit
records without rates, ``+R4``) or with a matrix it does not have (``Q.pfam``, ``HKY`` and the other
special cases of GTR): it stops with ``Not_found`` or ``unknown model``. gPlace stops in the same
cases, with the fitted model and what to refit with (e.g. ``LG+G4``).


.. _advanced/full_model:

The complete fitted model (--full_model)
----------------------------------------

``--full_model`` places with everything an IQ-TREE or RAxML-NG fit specifies instead of the part pplacer
reads. The placements are then the fitting program's model, not pplacer's, so they differ from pplacer's
(even for LG with gamma rates: other tables, other gamma algorithm). gPlace computes them with pplacer's
algorithm; the likelihoods equal the fitting program's own (below).

.. code-block:: bash

    # IQ-TREE: the checkpoint <prefix>.ckp.gz is read next to the .iqtree file
    gplace place --tree run.treefile --ref_aln ref.fasta --stats run.iqtree --full_model \
        --msa_file queries.fasta --out_file queries.jplace
    # RAxML-NG: <prefix>.raxml.bestModel is read next to the log
    gplace place --tree run.raxml.bestTree --ref_aln ref.fasta --stats run.raxml.log --full_model \
        --msa_file queries.fasta --out_file queries.jplace
    # with a reference package: give the model file
    gplace place --refpkg bac.refpkg --full_model run.ckp.gz --msa_file queries.fasta --out_file queries.jplace

IQ-TREE
^^^^^^^

The fitted parameters come from the fit's checkpoint (``<prefix>.ckp.gz``, IQ-TREE's full precision), the
model's name from ``<prefix>.iqtree``:

* every reversible DNA model (JC, F81, K2P, HKY, TN, K3P, TIM, TVM, SYM, GTR, 6-digit codes) and every
  reversible protein matrix of IQ-TREE 3 (LG, WAG, JTT, Q.pfam, Q.yeast, Q.insect, ..., mtREV, HIVb, FLU,
  FLAVI, Blosum62, Poisson; ``gplace/engine/data/iqtree_protein_models.json``, from IQ-TREE's source);
* the fitted frequencies (model, ``+F``, ``+FO``, ``+FQ``, ...);
* rate variation: ``+G`` (mean or median rates, ``--gamma-median``, computed with IQ-TREE's own gamma
  algorithm), ``+I``, ``+I+G``, ``+R``, ``+I+R``;
* protein mixtures: C10 to C60 (alone or with a matrix, ``LG+C20``), EX2, EX3, EHO, UL2, UL3, LG4M, LG4X and
  ``MIX{...}``, crossed with ``+G``/``+R`` or paired with them (``*G4``, ``*R4``).

For 2,940 placed queries over 87 fits (78 matrix / rate models and 9 mixtures; 4 fits with median gamma
rates), IQ-TREE's log-likelihood of the tree with the query attached where gPlace placed it, with the model
and branch lengths fixed, equals gPlace's to IQ-TREE's printed precision (at most 1e-8 for the DNA fits,
1e-6 for the protein fits; relative 9e-12); for the mixtures, IQ-TREE also reproduces each fit's
log-likelihood exactly from the model as gPlace reads it. A mixture costs about as many times more as it
has rate categories (C10 with 4 gamma rates: 40 categories).

RAxML-NG
^^^^^^^^

RAxML-NG writes the fitted model to ``<prefix>.raxml.bestModel`` (its own model syntax with every fitted
value), which RAxML-NG reads back as it is; gPlace reads it the same way, and the log (``<prefix>.raxml.log``)
for RAxML-NG's version, log-likelihood and counts. The file has 6 decimals by default, which can move the
likelihood a little from the fit's (by up to 2e-5 in our fits; 8e-3 for 190 user-given PROTGTR rates):
``raxml-ng --precision 12`` writes the model and the tree with 12 decimals.

* every DNA model (JC, K80, F81, HKY, TN93, TN93ef, K81, K81uf, TPM2/3, TPM2uf/3uf, TIM1/2/3, TIM1uf/2uf/3uf,
  TVM, TVMef, SYM, GTR, their aliases and ``DNA`` symmetry strings such as ``DNA010020``) and every protein
  matrix of RAxML-NG 2 (LG, WAG, JTT, the Q matrices, DAYHOFF, ..., FLAVI; ``gplace/engine/data/
  raxmlng_protein_models.json``, from RAxML-NG's source), fitted protein rates (``PROTGTR``);
* frequencies: the model's, ``+FU{...}`` (also what ``+FO`` fits are written as), ``+FE``, and ``+FC``,
  counted from the reference alignment as RAxML-NG counts them (checked against the counts in its log);
* rate variation: ``+G`` with mean or median rates (``+G4m``, ``+G4a``), computed with RAxML-NG's own gamma
  algorithm, ``+I`` (also ``+IC``, counted from the alignment), ``+R``, ``+I+G``, ``+I+R``;
* the mixtures LG4M and LG4X (also with ``+FC``, which RAxML-NG fits only with ``--force``).

For 3,380 placed queries over 89 RAxML-NG 2.0.3 fits (every DNA model and protein matrix with ``+G4``, and
the frequency, rate and mixture variants), RAxML-NG's log-likelihood of the tree with the query attached
where gPlace placed it (``raxml-ng --evaluate --opt-model off --opt-branches off``, the bestModel values)
equals gPlace's to the precision of the jplace file (12 significant digits): at most 5e-9 for the DNA fits,
6.8e-7 for the protein fits (relative 5e-12). The ``+FC`` frequencies and ``+IC`` proportions gPlace
counts equal RAxML-NG's to 5e-13.

Not supported yet
^^^^^^^^^^^^^^^^^

With ``--full_model``:

* invariable sites inside a mixture (e.g. ``LG+C20+I+G4``, ``LG4M+I``);
* partitioned models (IQ-TREE ``-p``, ``-q``, ``-Q``; RAxML-NG partition files): linked or proportional
  branch lengths are feasible (each site would carry every partition's rate categories, so cost grows with
  their total); unlinked branch lengths have no single placement position;
* IQ-TREE's fitted protein exchangeabilities (GTR20), user-defined matrices (``-mdef``, model files) and
  DNA mixtures; RAxML-NG's custom protein symmetries (``PROT...``) and frequencies of 1e-6 or less
  (RAxML-NG takes such states out of the rate matrix);
* models outside reversible DNA and protein substitution: codon, binary, multi-state, genotype and
  morphology models, non-reversible models (Lie Markov, UNREST, NONREV, NQ matrices), heterotachy (``+H``,
  GHOST), PoMo, ascertainment bias (``+ASC``) and site frequency models (PMSF, ``--site-freq``).
