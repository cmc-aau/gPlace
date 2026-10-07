pip
===

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

On other platforms pip falls back to the source package and compiles the engine (this needs
``gcc`` and ``make``); the environment variables listed in :ref:`installing` apply.
