import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")

# ------------------------------------------------------------

data = np.loadtxt(output+'output/spectrum.dat')
wavelength = data[:,0]
flux = data[:,1]

plt.plot(wavelength, flux, ls='-')

plt.xlim(np.amin(wavelength),np.amax(wavelength))
plt.ylim(1.e-27,1.e-11)

plt.xscale('log')
plt.yscale('log')

plt.xlabel('Wavelength [micron]', fontsize=12)
plt.ylabel('Flux [W/m$^2$/micron]', fontsize=12)

plt.savefig(output+'plot/spectrum_emission.pdf', bbox_inches='tight')