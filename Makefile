.PHONY: help install clean docs

debug = false
linux = false

FC = gfortran
DBG = -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fdump-core -fbacktrace -finit-real=nan -fbounds-check -Wline-truncation -Wcharacter-truncation
DBG += -Wsurprising -Waliasing -Wunused-parameter -fall-intrinsics -ffree-form -fdump-fortran-optimized -ffpe-trap=invalid,zero
PAR = -fopenmp

ifeq ($(debug),true)
	FLAGS = -O0 $(PAR) $(DBG)
else
	FLAGS = -O3 $(PAR)
endif

ifeq ($(linux),true)
	LIBS = -L lib -l cfitsio
else
	LIBS = -L lib -l cfitsio.5
endif

help:
	@echo "install - install ARTES"
	@echo "clean - remove artifacts"
	@echo "docs - generate documentation"

install:
	@echo \####################################################
	@echo \ \ \ \ \ \ \ \ \ \ \  _ \   \  ___ \  _____ \ \ ___ \ \ ___ 
	@echo \ \ \ \ \ \ \ \ \ \   \/_\\  \  \| _ \\ \|_ \  _\| \| __\| \/ __\|
	@echo \ \ \ \ \ \ \ \ \  \/ _ \\ \ \| \  \/ \  \| \| \  \| _\| \ \\__ \\
	@echo \ \ \ \ \ \ \ \ \ \/_\/ \\_\\ \|_\|_\\ \  \|_\| \  \|___\| \|___\/
	@echo
	@echo Atmospheric Radiative Transfer for Exoplanet Science
	@echo
	@echo ----------------------------------------------------
	@echo
	@echo Creating folders...
	@echo
	mkdir -p input/
	mkdir -p output/
	mkdir -p bin/
	mkdir -p lib/
	@echo
	@echo ----------------------------------------------------
	@echo
	@echo Downloading data...
	@echo
	wget -q --show-progress -O lib/libcfitsio.5.dylib https://people.phys.ethz.ch/~stolkert/ARTES/libcfitsio.5.dylib
	wget -q --show-progress -O lib/libcfitsio.so.3 https://people.phys.ethz.ch/~stolkert/ARTES/libcfitsio.so.3
	wget -q --show-progress -O dat/molecules/molecules.tar.gz https://people.phys.ethz.ch/~stolkert/ARTES/molecules.tar.gz
	wget -q --show-progress -O bin/ComputePartMac https://people.phys.ethz.ch/~stolkert/ARTES/ComputePartMac
	wget -q --show-progress -O bin/ComputePartLinux https://people.phys.ethz.ch/~stolkert/ARTES/ComputePartLinux
	@echo
	@echo ----------------------------------------------------
	@echo
	@echo Compiling...
	@echo
	$(FC) $(FLAGS) -c src/ARTES.f90
	$(FC) $(FLAGS) -o bin/ARTES ARTES.o $(LIBS)
	@echo
	@echo ----------------------------------------------------
	@echo
	@echo Unpacking data...
	@echo
	tar zxf dat/molecules/molecules.tar.gz -C dat/molecules/
	rm -f dat/molecules/molecules.tar.gz
	@echo
	@echo ----------------------------------------------------
	@echo
	@echo Finished!
	@echo
	@echo \####################################################

clean:
	rm -f ARTES.o
	rm -rf docs/_build

docs:
	sphinx-apidoc -o docs/ src
	$(MAKE) -C docs clean
	$(MAKE) -C docs html