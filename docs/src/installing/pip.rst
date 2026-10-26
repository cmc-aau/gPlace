pip
===

.. note::

   Not available yet: gPlace is not on PyPI. Until it is, install it from GitHub
   (:ref:`installing`). This page describes the planned package.

.. code-block:: bash

    python -m pip install gtdb-gplace
    gplace check_install

On x86-64 Linux, pip installs a **prebuilt wheel** (``py3-none-manylinux2014_x86_64``): no
compiler is needed, and the same wheel works on every distribution with glibc >= 2.17 (CentOS 7
and newer) and every Python >= 3.12. It is built in the manylinux2014 container with GSL 1.16
linked in, and published only after ``gplace check_install`` and ``gplace test`` pass inside that
container and on other distributions (Ubuntu, Debian 11, Debian 12).

Why one wheel is exact everywhere: the engine brings its own ``exp``, ``log``, ``hypot``, ``log1p`` and GSL;
from the system it only uses ``sqrt``, ``frexp`` and ``ldexp``, which IEEE 754 defines exactly.

On macOS and other platforms pip falls back to the source package and compiles the engine (this
needs a C compiler and ``make``: on macOS, ``xcode-select --install``; it takes a few minutes, mostly
building GSL 1.16); the environment variables listed in :ref:`installing` apply.
