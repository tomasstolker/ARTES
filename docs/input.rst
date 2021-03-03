.. _input:

Input files
===========

The working folder should contain the following files:

* **artes.in** - Text file that contains the input parameters for ARTES (see :ref:`artes.in` for an example). This input file may have a different filename.

* **atmosphere.in** - Configuration file (INI format) that contains the input parameters for *tools/atmosphere.py* (see :ref:`atmosphere.in` for an example), which can be used to create the required *atmosphere.fits*. This input file is not mandatory but only required if *atmosphere.py* is used. This input file may have a different filename.

* **atmosphere.fits** - FITS file that contains the 3D atmosphere structure that is used by ARTES. The file can be manually created as long as the correct data structure is used. Alternatively, *tools/atmosphere.py* can be used  with an input file such as *atmosphere.in* which helps to build *atmosphere.fits*.

* **pressure_temperature.dat** - Optional file with a 1D P-T structure. This input file may have a different filename.

* **opacity/[FITS]** - The *opacity* folder should contain the opacity files that will be used by *atmosphere.in* to build *atmosphere.fits*.

Further details on these input files is provided below.

artes.in
--------

This file contains the input parameters for ARTES. A full description of all possible keywords is provided in the :ref:`artes.in` section. Command line keywords can be included with the ``-k`` flag, which will overrule the input file keyword.

atmosphere.in
-------------

The input from the *atmosphere.in* file is used by *tools/atmosphere.py* to create the required *atmosphere.fits* file, which contains the atmospheric structure. A complete overview of all possible parameters is provided in the :ref:`atmosphere.in` section. Here we show a simple example:

.. code-block:: ini

    [grid]
    radius: 1 ; [Rjup]
    radial: 100, 200 ; [km]
    theta: 60, 120 ; [deg]
    phi: 120, 240 ; [deg]

    [composition]
    gas: off
    fits01: gas.fits
    fits02: clouds.fits
    ; FITS number, density (g cm-3), r_in (km), r_out (km), theta_in (deg), theta_out (deg), phi_in (deg), phi_out (deg)
    opacity01: 1, 1.e-3, 0, nr, 0, ntheta, 0, nphi
    opacity02: 2, 1.e-1, 0, nr, 1, 2, 0, nphi

The *[grid]* section contains the grid structure of the atmosphere. Radial, polar, and azimuthal grid boundaries can be added.

The *[composition]* section contains a list of all the opacity FITS files that are used. Numbers counting from one, with a leading zero for single digit numbers. The opacity keywords specify which opacity sources belong to which grid cells.

Specifically, *FITS number* is the number of the opacity FITS file, *density* is the density of the grid cell, *r_in* and *r_out* are the inner and outer radial boundary, *theta_in* and *theta_out* are the inner and outer polar boundary, and *phi_in* and *phi_out* are the inner and outer azimuthal boundary. The outer most boundaries may also be specified as *nr*, *ntheta*, and *nphi*.

Gas opacities are read automatically by setting *gas: on*. The surface gravity, $\log(g)$ (with g in cm s-2), and mean molecular weight (in g/mol) should be specified in case a P-T profile is used (e.g. as *pressure_temperature.dat*) to set up the radial density structure.

A self-luminous circumplanetary disk can be added as:

.. code-block:: ini

    ; FITS number, density (g cm-3), temperature, r_in (km), r_out (km), theta_in (deg), theta_out (deg), dust2gas, gas_abs (cm2 g-1)
    ring: 1, 1e-5, 100., 2e4, 1e5, 3, 4, 1e-2, 1e-2

Make sure to use the following keyword in *artes.in*:

.. code-block:: ini

  planet:ring=on

atmosphere.fits
---------------

This FITS file contains the atmospheric structure and scattering properties and should have the following file structure:

  0. 1D Radial boundaries (m)
  1. 1D Polar boundaries (deg)
  2. 1D Azimuthal boundaries (deg)
  3. 1D Wavelength points (um)
  4. 3D Density (kg m-3)
  5. 3D Temperature (K)
  6. 4D Scattering opacity (m-1)
  7. 4D Absorption opacity (m-1)
  8. 6D Scattering matrix
  9. 4D Asymmetry parameter

The radial boundaries are included as the primary HDU of the FITS file and the 9 following extensions are image HDUs.

To run ARTES, the *atmosphere.fits* and *artes.in* files are required. The *atmosphere.fits* file can be created with the *tools/atmosphere.py* script and an *atmosphere.in* input file.

Alternatively, the user could also manually create *atmosphere.fits*, for example by adopting the atmospheric structure from a different model and using ARTES for calculating the polarization observables.

.. important::
	The extension with the 3D density structure is no longer required by ARTES. The density is already included in the extensions with the scattering and absorption opacities, which are the product of the particle opacity and mass density. Therefore, the density array may simply contain zeros.

pressure_temperature.dat
------------------------

A pressure-temperature profile can be provided in the folder where also the opacity folder is located. The profile is used by ARTES to compute the gas densities, mixing ratios, and absorption cross sections. The profile should be given in units of bar and K with increasing pressure.

.. important::
   When using a P/T profile, the radii (in km) corresponding to the pressure layers are calculated with `tools/atmosphere.py`. Therefore, no values should be provided to the ``radial`` keyword in the `atmosphere.in` configuration file.

Scattering properties
---------------------

Several type of opacities can be generated. The opacity and scattering matrices need to be provided in a FITS format in which the first extension contains the wavelength dependent extinction, absorption, and scattering opacity, and the second extension contains the wavelength-dependent, 16-element scattering matrices.

The tools/opacity.py module contains several functions to create the required FITS files for different particle types:

   1. **opacity_henyey** - Henyey-Greenstein scattering phase function.

   2. **opacity_rayleigh** - Rayleigh scattering phase function.

   3. **opacity_gas** - Gas opacities with Rayleigh scattering cross-section and wavelength dependent absorption coefficients.

   4. **opacity_molecules** - Pressure temperature dependent gas opacities with equilibrium chemistry mixing ratios.

   5. **opacity_dhs** - DHS or Mie opacities and scattering matrices. This wrapper calls ``ComputePart``, a tool developed by `Michiel Min <http://www.exoclouds.com/>`_.

      In case a segmentation fault appears when running this routine, then try:

      .. code-block:: console

        $ ulimit -s unlimited

   6. **opacity_isotropic** - Isotropic scattering phase function.

All opacity FITS files should be located in the *opacity* folder.
