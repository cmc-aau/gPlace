From source
===========

.. code-block:: bash

    git clone https://github.com/cmc-aau/gPlace.git
    cd gPlace
    python -m pip install .      # compiles the engine (builds GSL 1.16 once, ~1-3 minutes)

    gplace check_install         # bit-exactness checks against GTDB-Tk's pplacer binary
    gplace test                  # places the bundled genomes; compares with pplacer's output

Offline installs:

.. code-block:: bash

    GPLACE_GSL_TARBALL=/path/to/gsl-1.16.tar.gz python -m pip install .

Reusing an existing static GSL 1.16 build (configured with ``--disable-shared --with-pic``):

.. code-block:: bash

    GPLACE_GSL_PREFIX=/opt/gsl-1.16 python -m pip install .

The engine can be rebuilt later without reinstalling with :ref:`commands/build_engine`.
