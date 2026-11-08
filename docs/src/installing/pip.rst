pip packages
============

gPlace is on PyPI as `gtdb-gplace <https://pypi.org/project/gtdb-gplace/>`_; the command it
installs is ``gplace``.

.. code-block:: bash

    python -m pip install gtdb-gplace              # the latest version
    python -m pip install gtdb-gplace==0.1.1       # one version
    python -m pip install --upgrade gtdb-gplace    # upgrade
    gplace check_install

On Linux, pip installs a **prebuilt wheel** (``manylinux2014``, for x86-64 and ARM64): no compiler
is needed, and the same wheel works on every distribution with glibc >= 2.17 (CentOS 7 and newer)
and every Python >= 3.12. It is built in the manylinux2014 container with GSL 1.16 linked in, and
published only after ``gplace check_install`` and ``gplace test`` pass inside that container and on
other distributions (Ubuntu, Debian 11, Debian 12).

Why one wheel is exact everywhere: the engine brings its own ``exp``, ``log``, ``hypot``, ``log1p`` and GSL;
from the system it only uses ``sqrt``, ``frexp`` and ``ldexp``, which IEEE 754 defines exactly.

On macOS and other platforms, pip installs the source package and compiles the engine. This needs a C
compiler and ``make`` (on macOS, ``xcode-select --install``) and internet access to download GSL 1.16,
and takes a few minutes, mostly building GSL. The options for offline builds and the environment
variables in :doc:`source` apply.
