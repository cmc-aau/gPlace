******
gPlace
******

.. warning::

   **Work in progress.** gPlace is under active development and its results are not yet fully
   tested. It reproduces pplacer's placements on the data sets described in :ref:`validation`,
   but it has not been tested widely beyond them, and its newest features (``--full_model`` for
   IQ-TREE and RAxML-NG fits, macOS builds) have had the least use. Options and outputs may still
   change. Before relying on gPlace for published results, check it on your own data, for example
   by comparing its jplace files with pplacer's on a subset of your genomes
   (:ref:`commands/compare`). Please report problems on the
   `issue tracker <https://github.com/cmc-aau/gPlace/issues>`_.

**gPlace** places genomes on the GTDB reference trees **exactly as GTDB-Tk's pplacer does**,
while using a fraction of the memory and time. It is a re-implementation of
``pplacer -c <refpkg>`` (pplacer v1.1.alpha19, the build GTDB-Tk ships): the WAG + FastTree CAT
reference packages used by GTDB-Tk, and every other model pplacer supports (LG, JTT and GTR for
DNA; CAT or gamma rates; model or empirical frequencies, see :ref:`models`).

* **Identical output:** same edges, likelihoods, like-weight ratios and branch lengths.
  Every number in the jplace file is textually identical to pplacer's, as is the reference
  tree string.
* **Bounded memory:** under 44 GiB for the largest GTDB-Tk R226 class-level tree
  (~104 GB with pplacer).
* **Self-verifying:** ``gplace check_install`` compares the compiled engine bit for bit with
  values taken from GTDB-Tk's pplacer binary.

gPlace is open source and released under the
`GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

Please post questions and issues related to gPlace on the Issues section of the GitHub repository.


Running gPlace
==============

1. Install gPlace and build its engine (:ref:`installing`)

2. Access the help documentation :ref:`commands`, or view the program help menu: ``gplace -h``

* Note: Individual help can be accessed via the specific command, e.g.: ``gplace place -h``


Citing gPlace
=============
We encourage you to cite gPlace and the third-party software as described in :ref:`references`.


.. toctree::
   :caption: Getting started
   :maxdepth: 1

   announcements
   installing/index
   faq


.. toctree::
   :caption: Running gPlace
   :maxdepth: 1

   commands/index
   files/index


.. toctree::
   :caption: About
   :maxdepth: 1

   how_it_works
   models
   validation
   changelog
   references
