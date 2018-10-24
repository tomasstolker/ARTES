.. _input:

Input files
===========

The working folder should contain the following files: ::

  artes.in (different name possible)
  atmosphere.in (optional, different name possible)
  atmosphere.fits (can be created with atmosphere.in)
  pressure_temperature.dat (optional, different name possible)
  opacity/[FITS]

artes.in
--------

This file contains the input parameters for ARTES. A full description of all possible keywords is provided in the :ref:`artes.in` section. Command line keywords can be included with the '-k' flag which will overrule the input file keyword.

atmosphere.in
-------------

The atmosphere.in input is used by tools/atmosphere.py to create the required atmosphere.fits file which contains the atmospheric structure, opacities, scattering matrices, etc. A complete overview of all possible parameters is provided in the :ref:`atmosphere.in` section. Here we show a simple example: ::

    [grid]
    radius: 1 ; [Rjup]
    radial: 100, 200 ; [km]
    theta: 60, 120 ; [deg]
    phi: 120, 240 ; [deg]

    [composition]
    gas: off
    fits01: gas.fits
    fits02: clouds.fits
    ; FITS number, density [g cm-3], rIn, rOut, thetaIn, thetaOut, phiIn, phiOut
    opacity01: 1, 1.e-3, 0, nr, 0, ntheta, 0, nphi
    opacity02: 2, 1.e-1, 0, nr, 1, 2, 0, nphi

The [grid] part contains the grid structure of the atmosphere. Radial, polar, and azimuthal grid boundaries can be added.

The [composition] part contains a list of all the opacity FITS files that are used. Numbers counting from one, with a leading zero for single digit numbers. The opacity keywords specify which opacity sources belong to which grid cells:

    FITS, density [g cm-3], rIn, rOut, thetaIn, thetaOut, phiIn, phiOut

FITS gives the corresponding FITS file number, density the grid cell density, rIn/rOut the inner and outer radial grid cell face and the same for theta and phi. The outer most boundaries are given by nr, ntheta, nphi.

Gas opacities are read automatically by setting 'gas: on'. A logg [cm/s2] and mean molecular weight [g/mol] have to be specified in case a pressure temperature profile is given to set up the radial density structure.

A self-luminous circumplanetary disk can be added as: ::

    ; FITS number, density [g cm-3], temperature, rIn [km], rOut [km], thetaIn, thetaOut, dust2gas, gasAbs [cm2 g-1]
    ring: 1, 1e-5, 100., 2e4, 1e5, 3, 4, 1e-2, 1e-2

Make sure to use the following keyword in artes.in: ::

  planet:ring=on

atmosphere.fits
---------------

This FITS file contains the atmospheric structure and scattering properties. It should have the following nine HDU extensions:

  0. 1D Radial boundaries [m]
  1. 1D Polar boundaries [deg]
  2. 1D Azimuthal boundaries [deg]
  3. 1D Wavelength points [micron]
  4. 3D Density [kg m-3]
  5. 3D Temperature [K]
  6. 4D Scattering opacity [m-1]
  7. 4D Absorption opacity [m-1]
  8. 6D Scattering matrix
  9. 4D Asymmetry parameter
  
To run ARTES, the atmosphere.fits and artes.in files are required. The atmosphere.fits file can be created with the tools/atmosphere.py script and an atmosphere.in input file.

pressure_temperature.dat
------------------------

A pressure-temperature profile can be provided in the folder where also the opacity folder is located. The profile is used by ARTES to compute the gas densities, mixing ratios, and absorption cross sections. The profile should be given in units of [bar] and [K] with increasing pressure.

Scattering properties
---------------------

Several type of opacities can be generated. The opacity and scattering matrices need to be provided in a FITS format in which the first extension contains the wavelength dependent extinction, absorption, and scattering opacity, and the second extension contains the wavelength-dependent, 16-element scattering matrices.

The tools folder contains several tools to create the required FITS files for different particle types:

   1. opacity_henyey.py: Henyey-Greenstein scattering phase function.

   2. opacity_rayleigh.py: Rayleigh scattering phase function.

   3. opacity_gas.py: Gas opacities with Rayleigh scattering cross-section and wavelength dependent absorption coefficients.

   4. opacity_molecules.py: Pressure temperature dependent gas opacities with equilibrium chemistry mixing ratios.

   5. opacity_mie.py: Mie or DHS opacities and scattering matrices. This wrapper calls ComputePart, a tool developed by `Michiel Min <http://www.michielmin.nl/>`__.

      In case a segmentation fault appears when running this routine, then try: ::
      
        ulimit -s unlimited

   6. opacity_isotropic.py: Isotropic scattering phase function.

All opacity FITS files should be located in the opacity folder.
