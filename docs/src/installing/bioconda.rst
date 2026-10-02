Bioconda
========

A Bioconda package is planned. A recipe needs ``rust``, ``gcc`` and ``make`` at build time and the
``gsl-1.16.tar.gz`` source as a second source; its build script runs:

.. code-block:: bash

    GPLACE_GSL_TARBALL=gsl-1.16.tar.gz python -m pip install . --no-deps

and its test section runs ``gplace check_install`` and ``gplace test``.
