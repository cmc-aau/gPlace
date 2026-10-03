.. _installing:

Installing gPlace
=================

gPlace's placement engine is written in Rust (``gplace/engine/rust``) and links three third-party C
libraries that exactness depends on: GSL 1.16, CORE-MATH and glibc 2.19's ``hypot`` and ``log1p``. On x86-64 and
ARM64 Linux, ``pip install gtdb-gplace`` installs it prebuilt; installing from source compiles it
automatically.


Sources
-------

.. toctree::
   :maxdepth: 1

   source
   pip
   bioconda


Hardware requirements
---------------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Resource
     - Requirement
   * - CPU
     - x86-64 (AVX2 kernels with AVX2 and FMA, otherwise the portable engine) or ARM64 (NEON).
   * - Memory
     - Depends on the reference tree and the number of queries: about
       *2 x internal nodes x sites x 168 bytes + edges x queries x 8 bytes*.
       43.8 GiB for GTDB-Tk R226 class-level tree 2 with 22,061 queries. ``gplace benchmark`` prints the estimate for your data.
   * - CPUs
     - Any; placement is memory-bandwidth bound. On a 64-vCPU AMD EPYC server, 40 worker
       processes were fastest.


Build requirements
------------------

Only needed when installing from source (not with the prebuilt wheel):

* A Rust toolchain, 1.75 or newer (``cargo``; e.g. from https://rustup.rs or conda-forge ``rust``) for the engine,
  and a C compiler (``gcc``) and ``make`` for GSL 1.16, CORE-MATH and ``hypot``. ``GPLACE_CARGO``
  selects a ``cargo`` that is not on ``PATH``. The engine has no Rust dependencies, so the build
  works offline.
* Internet access to download ``gsl-1.16.tar.gz`` from ftp.gnu.org, or a local copy
  (``GPLACE_GSL_TARBALL``, SHA-256 verified), or an existing static GSL 1.16 build
  (``GPLACE_GSL_PREFIX``).


Python libraries
----------------

gPlace is designed for Python >=3.12 and requires the following libraries, which will be
automatically installed:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Library
     - Version
   * - `numpy <https://numpy.org/>`_
     - >= 1.26
   * - `tqdm <https://github.com/tqdm/tqdm>`_
     - ~= 4.67


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
     - Installation: local ``gsl-1.16.tar.gz`` (offline builds).
   * - ``GPLACE_GSL_PREFIX``
     - Installation: existing static, ``-fPIC`` GSL 1.16 build.
   * - ``GPLACE_BUILD_DIR``
     - Installation: scratch directory for the GSL build (default: ``~/.cache/gplace/build``).
   * - ``GPLACE_BUILD_JOBS``
     - Installation: parallel ``make`` jobs (default: 8).
