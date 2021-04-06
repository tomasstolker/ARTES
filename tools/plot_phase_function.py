import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from palettable.colorbrewer.qualitative import Paired_10


# ------------------------------------------------------------
# Input

op_folder = os.path.join(sys.argv[1], "")
n_wavel = int(sys.argv[2])

# ------------------------------------------------------------

angle = np.arange(0.5, 180., 1.)

# Phase function

j = 0

for op_file in os.listdir(op_folder):
    if op_file.endswith(".fits") and not op_file.startswith("gas_opacity_"):
    	fitsfile = op_folder + op_file
    	hdulist = fits.open(fitsfile)
    	data = hdulist[1].data
    	hdulist.close()
    	plt.plot(angle, data[:, 0, n_wavel], '-', color=Paired_10.hex_colors[j], label=op_file[:-5])
    	j += 1

if j > 0:
    plt.xlabel('Scattering angle (deg)')
    plt.ylabel('Phase function')
    plt.xlim(0., 180.)
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.savefig(os.path.join(op_folder, 'phase_function.pdf'), bbox_inches='tight')
    plt.clf()

# Single scattering polarization

j = 0

for op_file in os.listdir(op_folder):
    if op_file.endswith(".fits") and not op_file.startswith("gas_opacity_"):
    	fitsfile = op_folder + op_file
    	hdulist = fits.open(fitsfile)
    	data = hdulist[1].data
    	hdulist.close()
    	plt.plot(angle, -data[:,1,n_wavel]/data[:, 0, n_wavel], '-', color=Paired_10.hex_colors[j], label=op_file[:-5])
    	j += 1
    
if j > 0:
    plt.xlabel('Scattering angle (deg)')
    plt.ylabel('Single scattering polarization')
    plt.xlim(0., 180.)
    plt.legend(loc='upper left')
    plt.savefig(os.path.join(op_folder, 'polarization.pdf'), bbox_inches='tight')
