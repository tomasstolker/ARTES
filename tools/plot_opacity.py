import os
import sys

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

import numpy as np
from astropy.io import fits

# ------------------------------------------------------------
# Input

opacity = os.path.join(sys.argv[1], "")

# ------------------------------------------------------------

count = 0
    
for file in os.listdir(opacity):
    if file.endswith('.fits') and not "gas_opacity_" in file[:-5]:
        count = count + 1

colors = iter(plt.cm.rainbow(np.linspace(0,1,count)))

j = 0

for file in os.listdir(opacity):
    if file.endswith('.fits'):

        fitsfile = opacity+file
        
        hdulist = fits.open(fitsfile)
        hdu0 = hdulist[0]
        hdu1 = hdulist[1]
        hdulist.close()

        naxis = hdu0.header['NAXIS1']

        wavelength = np.zeros(naxis)
        absorption = np.zeros(naxis)
        scattering = np.zeros(naxis)

        data = fits.getdata(fitsfile,0)

        for i in range(naxis):
            wavelength[i] = data[0,i]
            absorption[i] = data[2,i]
            scattering[i] = data[3,i]

        if not "gas_opacity_" in file[:-5]:

            c = next(colors)

            if np.size(wavelength) == 1:
                plt.plot(wavelength, scattering, color=c, marker='o', ms=5, mew=0, label=file[:-5])
                plt.plot(wavelength, absorption, color=c, marker='^', ms=5, mew=0)
            else:
                plt.plot(wavelength, scattering, color=c, ls='--', label=file[:-5])
                plt.plot(wavelength, absorption, color=c, ls='-')
                
                j += 1

        else:

            if np.size(wavelength) == 1:
                plt.plot(wavelength, scattering, color='gray', marker='o', ms=3, mew=0)
                plt.plot(wavelength, absorption, color='gray', marker='^', ms=3, mew=0)
            else:
                plt.plot(wavelength, scattering, color='gray', ls='--', lw=0.6)
                plt.plot(wavelength, absorption, color='gray', ls='-', lw=0.6)

plt.xlabel('Wavelength [micron]')
plt.ylabel('Opacity [cm$^2$/g]')

plt.xscale('log')
plt.yscale('log')

if count != j:
    plt.legend(loc='upper left')

plt.savefig(os.path.join(opacity, 'opacity.pdf'), bbox_inches='tight')
