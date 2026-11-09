.. _installing:

Installing gPlace
=================

gPlace runs on 64-bit Linux (x86-64 or ARM64) and macOS, with Python 3.12 or newer. It is on
Bioconda and PyPI as ``gtdb-gplace``; the command it installs is ``gplace``.

Install with conda
------------------

.. code-block:: bash

    conda create -n gplace -c conda-forge -c bioconda gtdb-gplace
    conda activate gplace

Bioconda's packages are prebuilt for Linux (x86-64 and ARM64) and macOS (Intel and Apple silicon):
nothing is compiled (see :doc:`bioconda`).

Install with pip
----------------

.. code-block:: bash

    conda create -n gplace -c conda-forge python=3.12    # optional: a fresh environment
    conda activate gplace

    python -m pip install gtdb-gplace

On Linux (x86-64 and ARM64), pip installs a prebuilt package: nothing is compiled. On macOS, pip
compiles gPlace's engine during the install, which needs a C compiler and ``make``
(``xcode-select --install``) and takes a few minutes (see :doc:`pip`); conda avoids this.

Check the installation
----------------------

.. code-block:: bash

    gplace check_install     # should end with "Installation OK"
    gplace test              # places 20 bundled genomes and compares them with pplacer's output

To upgrade later, the way you installed it:

.. code-block:: bash

    conda update -n gplace -c conda-forge -c bioconda gtdb-gplace    # conda
    python -m pip install --upgrade gtdb-gplace                       # pip


What you need to run it
-----------------------

* **Memory:** depends on the tree. Up to 44 GiB for the GTDB class-level trees; for the full
  GTDB bacterial tree, set a limit such as ``--max_memory 64G`` (see :ref:`advanced/memory`).
  ``gplace benchmark`` estimates the memory for your data.
* **CPUs:** any number. ``--cpus auto`` picks a good number of workers for your computer.


More details
------------

.. toctree::
   :maxdepth: 1

   bioconda
   pip
   source
