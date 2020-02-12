.. _installation:

Installation
============

Compiler
--------

Before installing ARTES, the GNU Fortran compiler needs to be installed, for example with Homebrew (Mac):

.. code-block:: console

    $ brew install gcc

Or with APT (Linux):

.. code-block:: console

    $ sudo apt-get install gfortran

ARTES
-----

ARTES is now compiled with help of the `Makefile`:

.. code-block:: console

    make all

Or with the ``linux`` flag on a Linux operating system:

.. code-block:: console

    make all linux=true

Note that during installation several files with opacities and libraries (14 MB) will be downloaded which requires the use of ``wget``.

.. tip::
   To see the different options of the `Makefile`, simply run ``make`` in the folder of the `Makefile`.

Library path
------------

Finally, the library path needs to be exported as environment variable. For example, on a Mac:

.. code-block:: console

    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:/path/to/ARTES/lib"

Or on Linux:

.. code-block:: console

    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/path/to/ARTES/lib"
