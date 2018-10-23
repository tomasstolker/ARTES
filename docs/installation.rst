.. _compilation:

Compilation
===========

ARTES has to be compiled with the GNU Fortran compiler which can for example be installed as: ::

    brew install gcc (Mac)
    sudo apt-get install gfortran (Linux)

Make sure to export the library path: ::

    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/ARTES/lib" (Mac)
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/ARTES/lib" (Linux)

ARTES can be compiled from the main folder where the Makefile is located: ::

    make install (Mac)
    make install linux=true (Linux)
