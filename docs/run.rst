.. _run:

Running ARTES
=============

To build the atmospheric structure with the `atmosphere.in` file:

.. code-block:: console

    $ python /path/to/tools/atmosphere.py atmosphere.in

And to run the radiative transfer:

.. code-block:: console

    $ /path/to/bin/artes [input] [photons] -o [output] -k [keyword]=[value]

Several plot scripts are available in the tools folder which can be run directly on the output folder. For example:

.. code-block:: console

    $ python /path/to/tools/plot_opacity.py [output]
    $ python /path/to/tools/plot_phase_curve.py [output]
