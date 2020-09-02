.. _opacities:

Creating opacities
==================

There are several functions available in `tools/opacity.py` that can be used for creating opacities and scattering matrices. To use these functions, Python needs to know where the code is located so it is recommended to add the `tools` folder to the ``PYTHONPATH`` (see :ref:`installation` section). For details on the functions and the parameters, please have a look at the docstrings in the `tools/opacity.py` file.

.. important::
   ARTES runs the radiative transfer at the wavelengths that are found in the opacity files. Therefore, it is important to make sure that the opacities have been calculated at the required wavelengths. Also, when using multiple opacity files, the same wavelengths should be used in all opacity files. For imaging (``detector:type=imaging``) and phase curves (``detector:type=phase``), the first wavelength of the opacity files is used in case multiple wavelengths are present.

Let's now create some opacities for Rayleigh scattering particles:

.. code-block:: python

   import opacity

   opacity.opacity_rayleigh(wavelength=(0.5, 5., 100),
                            output='rayleigh.fits',
                            albedo=0.5,
                            depolarization=0.,
                            mmw=2.)

A FITS files with two arrays is created. The primary array contains the wavelengths and the absorption, scattering, and extinction opacities. The second array contains the 4 by 4 scattering matrices for each wavelength and 180 scattering angles. This file can be placed in `opacity` folder of the working folder such that it is found when running `tools/atmosphere.py` (see :ref:`input` section).

Similarly we can create the scattering and absorption properties for icy ammonia particles. Here we need to provide the refractive indices, of which a few are available in the `dat/refractive` folder. We set the ``wavelength`` argument to ``None`` such that the wavelengths of the file with refractive indices are used.

.. code-block:: python

   import opacity

   opacity.opacity_dhs(ri_file='/path/to/ARTES/dat/refractive/ammonia_ice.dat')
                       wavelength=None,
                       output='dhs.fits',
                       density=1.,
                       fmax=0.2,
                       r_eff=1.0,
                       v_eff=0.1)
