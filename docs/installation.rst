.. _installation:

Installation
============

The GNU Fortran compiler has to be installed, for example: ::

    brew install gcc (Mac)
    sudo apt-get install gfortran (Linux)

The library path should be exported as environment variable, for example: ::

    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/artes/lib" (Mac)
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/artes/lib" (Linux)

ARTES can then be compiled with help of the Makefile: ::

    make all (Mac)
    make all linux=true (Linux)

Note that during installation several files (14 MB) will be downloaded which are too large for Github.