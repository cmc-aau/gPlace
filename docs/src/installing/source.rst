Installing from source
======================

gPlace's placement engine is written in C (``gplace/engine/csrc/xp.c``) and links three third-party C
libraries that exactness depends on: GSL 1.16, CORE-MATH and glibc 2.19's ``hypot`` and ``log1p``.
Installing from source compiles them; the results are identical on every platform and with every
compiler.


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
    gplace check_install

The build downloads one library (GSL 1.16) from a GNU mirror; without internet access, see below.


Build requirements
------------------

* A C compiler (gcc 4.8 or newer, e.g. CentOS 7's system gcc, or clang) and ``make``. On macOS,
  Apple's command line tools provide both: ``xcode-select --install``. The build stops with a
  message if the compiler lacks a builtin the C code needs, and checks that each engine loads.
* Internet access to download ``gsl-1.16.tar.gz`` from a GNU mirror (ftpmirror.gnu.org, ftp.gnu.org or
  mirrors.kernel.org, tried in turn), or a local copy, or an existing static GSL 1.16 build (below).
* Python >= 3.12; ``numpy`` (>= 1.26) and ``tqdm`` (~= 4.67) are installed automatically.

Offline installs:

.. code-block:: bash

    GPLACE_GSL_TARBALL=/path/to/gsl-1.16.tar.gz python -m pip install .

Reusing an existing static GSL 1.16 build (configured with ``--disable-shared --with-pic``):

.. code-block:: bash

    GPLACE_GSL_PREFIX=/opt/gsl-1.16 python -m pip install .

The engine can be rebuilt later without reinstalling with :ref:`commands/build_engine`.


Environment variables
---------------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Variable
     - Meaning
   * - ``GPLACE_ENGINE_DIR``
     - Directory holding ``libgplace_core.so`` (default: the installed package directory).
   * - ``GPLACE_GSL_TARBALL``
     - Installation: local ``gsl-1.16.tar.gz`` (offline builds, SHA-256 verified).
   * - ``GPLACE_GSL_PREFIX``
     - Installation: existing static, ``-fPIC`` GSL 1.16 build.
   * - ``GPLACE_BUILD_DIR``
     - Installation: scratch directory for the GSL build (default: ``~/.cache/gplace/build``).
   * - ``GPLACE_BUILD_JOBS``
     - Installation: parallel ``make`` jobs (default: 8).


Hardware details
----------------

* **CPU:** x86-64 (AVX2 kernels with AVX2 and FMA, otherwise the portable engine) or ARM64 (NEON),
  on Linux or macOS (Apple Silicon: NEON; Intel Macs: AVX2).
* **Memory:** depends on the reference tree and the number of queries: about
  *2 x internal nodes x sites x 168 bytes + edges x queries x 8 bytes* (with ``--ranked_edges``
  the second term stops growing past 20,000 queries); 43.8 GiB for GTDB R226 class-level tree 2
  with 22,061 queries.
* **CPUs:** placement is limited by memory bandwidth; on a 64-vCPU AMD EPYC server, 40 worker
  processes were fastest.


Unit tests
----------

From a source checkout, after ``pip install .`` (or ``pip install -e .``):

.. code-block:: bash

    python -m unittest discover -s tests

Tests that need the compiled engine are skipped when it is missing.
