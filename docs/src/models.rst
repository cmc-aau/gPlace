.. _models:

Models
======

gPlace reproduces every evolutionary model pplacer v1.1.alpha19 supports, bit for bit. Like
pplacer, it reads the model from the reference package's ``phylo_model`` JSON file (pplacer's
``-m`` option is ignored for reference packages, so gPlace has no model option either).

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Part of the model
     - ``phylo_model`` field
     - Values
   * - Substitution model
     - ``subs_model``
     - ``WAG``, ``LG``, ``JTT`` (protein; the exchangeabilities and frequencies inside pplacer)
       or ``GTR`` (DNA; ``subs_rates``: ``ac ag at cg ct gt``)
   * - Frequencies
     - ``empirical_frequencies``
     - ``false``: the model's own frequencies; ``true``: counted from the reference alignment.
       GTR always uses empirical frequencies, as in pplacer
   * - Rate variation across sites
     - ``ras_model``
     - ``Price-CAT`` (FastTree CAT: one rate per site, ``Price-CAT.Rates`` and
       ``SiteCategories``) or ``gamma`` (``gamma.alpha``, ``gamma.n_cats``: every site under all
       discrete-gamma rates, averaged)

GTDB-Tk's reference packages are WAG + Price-CAT with model frequencies. Reference packages built
with ``taxit create`` from FastTree logs (CAT), RAxML/PhyML statistics (gamma), or IQ-TREE and RAxML-NG
fits (gamma) work as they are; for IQ-TREE and RAxML-NG see :ref:`commands/place` (which parts of a
fitted model pplacer uses). ``place --full_model`` places with an IQ-TREE or RAxML-NG fit's complete model
instead (the program's matrices, frequencies, +G/+I/+R and protein mixtures; results then differ from
pplacer's).

What the parts mean
-------------------

* **Substitution model**: how easily each pair of letters replaces each other (190
  exchangeabilities for amino acids, 6 for nucleotides). WAG, LG and JTT are published tables
  estimated from large protein databases; GTR's six rates are fitted to the reference alignment
  when the tree is built.
* **Frequencies**: the long-run frequency of each letter: the model's published values, or values
  counted from the reference alignment (useful for unusual compositions).
* **Rate variation**: some alignment columns evolve much faster than others. CAT gives every column
  one of about 20 fixed rates, chosen when the tree was built; gamma assumes rates follow a gamma
  distribution with shape alpha, cut into 4 (or more) equally likely categories, and averages each
  column over all of them. Gamma costs 4x the memory and time of CAT (n_cats x in general).

Use the model the reference tree was built with: its branch lengths are only meaningful under that
model.

How the gamma categories are computed
-------------------------------------

The rate of each gamma category is the mean (or median) of the gamma distribution over that category.
Programs compute it with different numerical routines, which agree to about 1e-8 to 1e-7 (relative) and
much less for very small alpha. gPlace uses the routine of the program the model comes from, so that its
likelihoods are that program's:

* reference packages and ``--stats`` files (pplacer's models): pplacer's routine (GSL 1.16), for results
  identical to pplacer's;
* ``place --full_model`` with an IQ-TREE fit: IQ-TREE's (Yang's algorithm, as in IQ-TREE's
  ``RateGamma``; mean or median rates);
* ``place --full_model`` with a RAxML-NG fit: RAxML-NG's (the same algorithm, as in its coraxlib library;
  it rescales the rates differently, which changes the last digits).

The IQ-TREE and RAxML-NG rates are bit-identical to the programs' own code.

DNA
---

GTR reference packages use a second engine compiled for 4 states (``libgplace_core_nt.so``), with
pplacer's nucleotide alphabet: ``A C G T``, the IUPAC ambiguity codes ``R Y M K S W H B V D N``,
and ``X``, ``-``, ``?`` as fully ambiguous. ``U`` is not accepted (as in pplacer). Queries must be
aligned to the reference alignment, as for protein.

Memory
------

The likelihood vectors hold *states x rates* numbers per site: 20 for protein CAT, 4 for DNA CAT,
80 for protein with 4 gamma rates. ``--max_memory`` and ``gplace benchmark`` take this into
account; tiled placement works for every model. Measured on the GTDB-Tk R226 class-level reference
package 8 (13,714 nodes, 343 genomes): WAG + CAT 33.1 GiB with pplacer, 11.0 GiB with gPlace;
WAG + gamma 125.7 GiB with pplacer, 42.2 GiB with gPlace untiled and 6.3 GiB with tiles of at most
2,000 edges.

.. _models/gamma_underflow:

pplacer's gamma underflow on large trees
----------------------------------------

To avoid numerical underflow, pplacer stores each site's likelihoods with a power-of-two exponent
and rescales a site when its numbers get small. With gamma, one exponent is shared by all rates,
and pplacer decides whether to rescale from the largest binary exponent (``frexp``) of the rates'
maxima. But ``frexp(0.0)`` returns exponent 0: once the slow rates of a fast-evolving site have
underflowed to exactly 0.0, which happens on trees of thousands of divergent sequences, they make
the site look large, the site is never rescaled again, the remaining rates underflow too, and the
site likelihood becomes 0. Every placement of a query with a residue at such a site then has
likelihood ``-inf`` (and like-weight ratio ``-nan``); pplacer only prints
``Warning: encountered -inf for final likelihood``.

On reference package 8 with WAG + gamma (alpha 0.8), pplacer gives ``-inf`` for 342 of 343 genomes.
gPlace reproduces this exactly by default and warns about it. With ``--fix_gamma_underflow``,
all-zero rates are ignored when the exponent is chosen: all 343 genomes get finite likelihoods,
and 332 of them land on the same edge as the WAG + CAT placement. Where pplacer's results are
finite, the fix leaves them identical or changes only their last digits (where the blocked
rescaling passed through subnormal numbers), never the chosen edges. CAT models are not affected.
