.. _installation:

Installation
============

GNU Fortran compiler
--------------------

Before installing ARTES, the GNU Fortran compiler needs to be installed, for example with `Homebrew <https://brew.sh/>`_ (Mac):

.. code-block:: console

    $ brew install gcc

Or with `APT <https://en.wikipedia.org/wiki/APT_(software)>`_ (Linux):

.. code-block:: console

    $ sudo apt-get install gfortran

Downloading from Github
-----------------------

Next, we will download ARTES by cloning the Github repository:

.. code-block:: console

    $ git clone git@github.com:tomasstolker/ARTES.git

New commits can be pulled from Github once a local copy of the repository exists:

.. code-block:: console

    $ git pull origin main

Compiling the source code
-------------------------

ARTES can now be compiled with help of the `Makefile` which is located in the main folder:

.. code-block:: console

    $ make all

On a Linux operating system, make sure to set the ``linux`` flag when running ``make all``:

.. code-block:: console

    $ make all linux=true

Note that during installation several files with opacities and libraries (14 MB) will be downloaded which requires the use of ``wget``.

.. tip::
   To see the different options of the `Makefile`, simply run ``make`` in the folder of the `Makefile`.

Setting the library path
------------------------

Finally, the library path needs to be exported as environment variable. For example, on a Mac:

.. code-block:: console

    $ export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:/path/to/ARTES/lib"

Or on Linux:

.. code-block:: console

    $ export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/path/to/ARTES/lib"

Setting the Python path
-----------------------

It is also recommended to add the `tools` folder to the ``PYTHONPATH`` such that the `tools/opacity.py` functions are easily accessible with Python:

.. code-block:: console

    $ export PYTHONPATH="$PYTHONPATH:/path/to/ARTES/tools"
