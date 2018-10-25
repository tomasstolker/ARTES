import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from palettable.colorbrewer.qualitative import Paired_10

# ------------------------------------------------------------
# Input

opacity = os.path.join(sys.argv[1], "")
wavelength = int(sys.argv[2])

# ------------------------------------------------------------

angle = np.zeros(180)
for i in range(180):
    angle[i] = float(i)+0.5

# Phase function

j = 0

for file in os.listdir(opacity):
    if file.endswith(".fits") and not file.startswith("gas_opacity_"):

	fitsfile = opacity+file
	hdulist = fits.open(fitsfile)
	data = hdulist[1].data
	hdulist.close()
	plt.plot(angle, data[:,0,wavelength], '-', color=Paired_10.hex_colors[j], label=file[:-5])
	j += 1

if j > 0:

    plt.xlabel('Scattering angle [degrees]')
    plt.ylabel('Phase function')
    plt.xlim(0,180)
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.savefig(os.path.join(opacity, 'phase_function.pdf'), bbox_inches='tight')
    plt.clf()

# Single scattering polarization

j = 0

for file in os.listdir(opacity):
    if file.endswith(".fits") and not file.startswith("gas_opacity_"):

    	fitsfile = opacity+file
    	hdulist = fits.open(fitsfile)
    	data = hdulist[1].data
    	hdulist.close()
    	plt.plot(angle, -data[:,1,wavelength]/data[:,0,wavelength], '-', color=Paired_10.hex_colors[j], label=file[:-5])
        
	j += 1
    
if j > 0:
	
    plt.xlabel('Scattering angle [degrees]')
    plt.ylabel('Single scattering polarization')
    plt.xlim(0,180)
    plt.legend(loc='upper left')
    plt.savefig(os.path.join(opacity, 'polarization.pdf'), bbox_inches='tight')
