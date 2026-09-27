.. _commands/build_engine:

build_engine
============

The engine is compiled automatically when gPlace is installed; ``build_engine`` is an advanced command to rebuild it without reinstalling (e.g. with a different GSL 1.16 build). It compiles ``libgplace_core.so``: it builds GSL 1.16 as a static, position-independent library (downloaded from ftp.gnu.org and checked against its SHA-256, or taken from ``--gsl_tarball`` / ``--gsl_prefix``), compiles the engine with glibc 2.19's ``hypot`` and CORE-MATH's correctly-rounded ``exp``/``log``, and then runs :ref:`check_install <commands/check_install>`. Every object is compiled with ``-ffp-contract=off``, so no multiply-add is fused.

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
