atmosphere.in
=============

This is an example of the available parameters in atmosphere.in: ::

    [grid]

    ; Planet radius [Rjup]
    radius: 1.

    ; Radial boundaries [km]
    radial: 100, 200

    ; Polar boundaries [deg]
    theta: 30, 75, 89.9, 90.1, 95, 150

    ; Azimuthal boundaries [deg]
    phi: 90, 180, 270

    [composition]

    ; Gas opacities (on/off)
    gas: on

    ; Mean molecular weight [g/mol]
    molweight: 2.3

    ; Log g [cm s-2]
    logg: 4.0

    ; FITS files with opacities and scattering matrices
    fits01: clouds.fits
    fits02: hazes.fits

    ; Atmospheric structure
    ; FITS number, density [g cm-3], rIn, rOut, thetaIn, thetaOut, phiIn, phiOut
    opacity01: 1, 1e-5, 0, 1, 2, 3, 2, nphi
    opacity02: 1, 2e-5, 0, 1, 4, ntheta, 0, 1
    opacity03: 2, 1e-6, 1, 2, 0, 1, 0, nphi

    ; Thermal circumplanetary disk
    ; FITS number, density [g cm-3], temperature, rIn [km], rOut [km], thetaIn, thetaOut, dust2gas, gasAbs [cm2 g-1]
    ring: 1, 1e-5, 100., 2e4, 1e5, 3, 4, 1e-2, 1e-2