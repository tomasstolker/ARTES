import os
import sys

import matplotlib.pyplot as plt
import numpy as np


# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], '')

# ------------------------------------------------------------

plot_dir = output+'plot/'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

data = np.loadtxt(output+'output/spectrum.dat')
wavelength = data[:, 0]
flux = data[:, 1]

plt.plot(wavelength, flux, ls='-')

plt.xlim(np.amin(wavelength), np.amax(wavelength))
plt.ylim(1.e-27, 1.e-11)

plt.xscale('log')
plt.yscale('log')

plt.xlabel('Wavelength (um)', fontsize=12)
plt.ylabel('Flux (W/m$^2$/um)', fontsize=12)

plt.savefig(plot_dir+'spectrum_emission.pdf', bbox_inches='tight')
