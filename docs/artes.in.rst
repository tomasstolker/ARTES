.. _artes.in:

artes.in
========

Here you find a list with all possible keywords for the input file that is required by ARTES. This file can have any name but is called `artes.in` throughout the documentation.

.. important::
   There are no spaces allowed between the keywords and values.

.. code-block:: ini

    ----------------------------------------------------------------------
    * General

    * Write log file (on) or screen output (off)
    general:log=off

    * Send e-mail when finished
    general:email=

    * Use random seed (on/off)
    general:random=on

    ----------------------------------------------------------------------
    * Photons

    * Photon source (star, planet)
    photon:source=star

    * Photon stopping parameter (0-1)
    photon:fstop=1d-5

    * Minimum photon energy (0 < minimum <= 1), otherwise removed after scattering
    photon:minimum=1d-20

    * Weight cell luminosity, equal photons emitted from all cells (on/off)
    photon:weight=off

    * Photon scattering (on/off)
    photon:scattering=on

    * Thermal photon emission direction (isotropic/biased)
    photon:emission=isotropic

    * Biased emission asymmetry parameter
    * 0 <= bias < 1
    photon:bias=0.8

    * Modified random walk parameter
    * walk > 0 -> on
    * walk < 0 -> off
    photon:walk=-1

    * Photon number for debug output, 0=off
    photon:number=0

    ----------------------------------------------------------------------
    * Star

    * Stellar temperature [K]
    star:temperature=5800

    * Stellar radius [R_sun]
    star:radius=1

    * Stellar direction [degrees]
    * 0 <= star:theta <= 180
    * 0 <= star:phi < 360
    star:theta=90
    star:phi=0

    ----------------------------------------------------------------------
    * Planet

    * Surface albedo (black=0, lambert=1)
    planet:albedo=0

    * Oblateness (0-1)
    planet:oblateness=0

    * Orbital radius [au]
    planet:orbit=5

    * Circumplanetary ring (on/off)
    planet:ring=off

    * Minimum distance to a cell boundary [m]
    planet:grid=1d-6

    * Maximum optical depth, tau > 0
    * tau < 0, recommended: tau=30 -> reflected, tau=5 -> emission
    planet:tau=-1

    ----------------------------------------------------------------------
    * Detector

    * Observation mode (spectrum, phase, imaging)
    detector:type=imaging

    * Detector location [degrees]
    * 0 <= detector:theta <= 180
    * 0 <= detector:phi < 360
    detector:theta=90
    detector:phi=90

    * Number of detector pixels in x and y direction
    detector:pixel=25

    * Distance [pc]
    detector:distance=10

    * Rotation angle [deg]
    * 0 < detector:angle < 360
    * angle < 0 -> off
    detector:angle=-1

    ----------------------------------------------------------------------
    * Output

    * Debug errors (on/off)
    output:debug=off

    * Global energy flow (on/off)
    output:global=off

    * Latitudinal energy flow (on/off)
    output:latitudinal=off
