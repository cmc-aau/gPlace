.. _commands/build_engine:

build_engine
============

The engine is compiled automatically when gPlace is installed; ``build_engine`` is an advanced command to rebuild it without reinstalling (e.g. with a different GSL 1.16 build). It compiles ``libgplace_core.so``: it builds GSL 1.16 as a static, position-independent library (downloaded from ftp.gnu.org and checked against its SHA-256, or taken from ``--gsl_tarball`` / ``--gsl_prefix``), compiles glibc 2.19's ``hypot`` and CORE-MATH's correctly-rounded ``exp``/``log`` (C, with ``-ffp-contract=off``), builds the engine itself (the Rust crate ``gplace/engine/rust``, with ``cargo``; Rust 1.75 or newer) and links them together, and then runs :ref:`check_install <commands/check_install>`. No multiply-add is fused anywhere: the C is compiled not to, and Rust never does.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: build_engine
   :nodefaultconst:

Example
-------

.. code-block:: bash

    gplace build_engine --jobs 16
