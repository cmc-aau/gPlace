Bioconda packages
=================

gPlace is on Bioconda as `gtdb-gplace <https://anaconda.org/bioconda/gtdb-gplace>`_; the command it
installs is ``gplace``.

.. code-block:: bash

    conda create -n gplace -c conda-forge -c bioconda gtdb-gplace           # the latest version
    conda create -n gplace -c conda-forge -c bioconda gtdb-gplace=0.1.0     # one version
    conda activate gplace
    gplace check_install

To upgrade an environment: ``conda update -n gplace -c conda-forge -c bioconda gtdb-gplace``.
``mamba`` takes the same commands. If conda is already set up for Bioconda (the ``conda-forge`` and
``bioconda`` channels, strict channel priority), the ``-c`` options can be left out.

The packages are **prebuilt** for Linux (x86-64 and ARM64) and macOS (Intel and Apple silicon):
nothing is compiled during the install, on macOS either. Bioconda builds them from the PyPI source
package, with GSL 1.16 compiled into the engine, and publishes a package only after its tests pass;
they run ``gplace check_install`` and ``gplace test``.

The recipe is `recipes/gtdb-gplace
<https://github.com/bioconda/bioconda-recipes/tree/master/recipes/gtdb-gplace>`_ in bioconda-recipes.
