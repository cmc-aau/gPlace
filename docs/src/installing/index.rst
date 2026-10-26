.. _installing:

Installing gPlace
=================

gPlace runs on 64-bit Linux (x86-64 or ARM64) and macOS, with Python 3.12 or newer.

.. note::

   gPlace is not on PyPI or Bioconda yet: install it from GitHub as below.


Install from GitHub
-------------------

You need a C compiler and ``make``, which most Linux systems have. If yours does not, conda can
provide them (first command below); on macOS, run ``xcode-select --install``.

.. code-block:: bash

    conda create -n gplace -c conda-forge python=3.12 gcc make    # optional: a fresh environment
    conda activate gplace

    git clone https://github.com/cmc-aau/gPlace.git
    cd gPlace
    pip install .            # also builds gPlace's engine: a few minutes, the first time only

Then check that it works:

.. code-block:: bash

    gplace check_install     # should end with "Installation OK"
    gplace test              # places 20 bundled genomes and compares them with pplacer's output

The build downloads one library (GSL 1.16) from a GNU mirror. Without internet access, or to
reuse an existing build, see :doc:`source`.


What you need to run it
-----------------------

* **Memory:** depends on the tree. Up to 44 GiB for the GTDB-Tk class-level trees; for the full
  GTDB bacterial tree, set a limit such as ``--max_memory 64G`` (see :ref:`advanced/memory`).
  ``gplace benchmark`` estimates the memory for your data.
* **CPUs:** any number. ``--cpus auto`` picks a good number of workers for your computer.


Other ways to install
---------------------

.. toctree::
   :maxdepth: 1

   source
   pip
   bioconda
