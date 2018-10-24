.PHONY: help start data compile end download install all clean docs

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
	@echo "all - download and install"
	@echo "download - download required files"
	@echo "install - install ARTES"
	@echo "clean - remove artifacts"
	@echo "docs - generate documentation"

start:
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
	mkdir -p bin/
	mkdir -p lib/
	mkdir -p dat/molecules/
	@echo
	@echo ----------------------------------------------------

data:
	@echo
	@echo Downloading data...
	@echo
	wget -q --show-progress -O lib/libcfitsio.5.dylib https://people.phys.ethz.ch/~stolkert/artes/libcfitsio.5.dylib
	wget -q --show-progress -O lib/libcfitsio.so.3 https://people.phys.ethz.ch/~stolkert/artes/libcfitsio.so.3
	wget -q --show-progress -O dat/molecules/molecules.tar.gz https://people.phys.ethz.ch/~stolkert/artes/molecules.tar.gz
	wget -q --show-progress -O bin/computepart_mac https://people.phys.ethz.ch/~stolkert/artes/computepart_mac
	wget -q --show-progress -O bin/computepart_linux https://people.phys.ethz.ch/~stolkert/artes/computepart_linux
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
	@echo Setting permissions...
	@echo
	chmod 700 bin/computepart_mac
	chmod 700 bin/computepart_linux
	@echo
	@echo ----------------------------------------------------

compile:
	@echo
	@echo Compiling...
	@echo
	$(FC) $(FLAGS) -c src/artes.f90
	$(FC) $(FLAGS) -o bin/artes artes.o $(LIBS)
	@echo
	@echo ----------------------------------------------------

end:
	@echo
	@echo Finished!
	@echo
	@echo \####################################################

download: start data end

install: start compile end

all: start data compile end

clean:
	rm -f artes.o
	rm -rf docs/_build

docs:
	sphinx-apidoc -o docs/ src
	$(MAKE) -C docs clean
	$(MAKE) -C docs html