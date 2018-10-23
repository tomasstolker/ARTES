.. _input:

Input files
===========

Folder structure
----------------

All input files have to be located in the input folder of which the file structure is the following: ::

  ARTES/input/[atmosphere]/artes.in
  ARTES/input/[atmosphere]/atmosphere.in
  ARTES/input/[atmosphere]/atmosphere.fits
  ARTES/input/[atmosphere]/pressure_temperature.dat (optional)
  ARTES/input/[atmosphere]/opacity/[opacityFITS]

Where [atmosphere] is a user-defined name.

Settings: artes.in
------------------

This file contains the input parameters for ARTES. A full description of all keywords is provided in the artes.in template file. Command line keywords can be provided with the '-k' flag which will overrule the input file keyword.

Atmosphere: atmosphere.fits
---------------------------

This FITS file contains the atmospheric structure and scattering properties. It should contain the following nine HDU extensions:

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
  
To run ARTES, the atmosphere.fits and artes.in files are required. Additional tools are provided to help create the atmosphere.fits file.

Pressure/temperature profile
----------------------------

A pressure-temperature profile can be provided in the [atmosphere] folder which is used by ARTES to determine the gas densities, mixing ratios, and absorption cross sections. The profile should be given in units of [bar] and [K] with increasing pressure.

Opacities and scattering matrices
---------------------------------

Several type of opacities can be generated. The opacity and scattering matrices need to be provided in a FITS format in which the first extension contains the wavelength dependent extinction, absorption, and scattering opacity, and the second extension contains the wavelength-dependent, 16-element scattering matrices.

The python folder contains a few opacity tools:

   1. opacity_henyey.py
      Henyey-Greenstein scattering phase function.

   2. opacity_rayleigh.py
      Rayleigh scattering phase function.

   3. opacity_gas.py
      Gas opacities with Rayleigh scattering cross-section and wavelength dependent absorption coefficients.

   4. opacity_molecules.py
      Pressure temperature dependent gas opacities with equilibrium chemistry mixing ratios.

   5. opacity_mie.py
      Mie or DHS opacities and scattering matrices. This wrapper calls ComputePart (M. Min, SRON). Make sure that the ComputePart binary file is executable: ::

        chmod 700 bin/ComputePart[Mac/Linux]
        
      In case a segmentation fault appears when running this routine, then try: ::
      
        ulimit -s unlimited

   6. opacity_isotropic.py
      Isotropic scattering phase function.

All opacity FITS files should be located in the opacity folder.

Create automatic input
----------------------

The atmosphere.in file has to be located in the [atmosphere] folder and its content should look something like: ::

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
