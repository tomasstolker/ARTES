import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")

# ------------------------------------------------------------

plotDir = output+'plot/'
dataDir = output+'output/'

if not os.path.exists(plotDir):
    os.makedirs(plotDir)

wavelength, stokesi, stokesq, stokesu, stokesv = np.loadtxt(dataDir+'spectrum.dat', unpack=True)
wavelength, norm1, norm2 = np.loadtxt(dataDir+'normalization.dat', unpack=True)

plt.plot(wavelength, stokesi, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]', fontsize=12)
plt.ylabel('Flux [W m$^{-2}$ micron$^{-1}$]', fontsize=12)
plt.savefig(plotDir+'spectrum_reflected_flux.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, stokesi/norm1, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]', fontsize=12)
plt.ylabel('Normalized Stokes I', fontsize=12)
plt.savefig(plotDir+'spectrum_reflected_flux_1.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, stokesi/norm2, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]', fontsize=12)
plt.ylabel('Normalized Stokes I', fontsize=12)
plt.savefig(plotDir+'spectrum_reflected_flux_2.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, np.sqrt(stokesq**2+stokesu**2+stokesv**2)/stokesi, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]', fontsize=12)
plt.ylabel('Degree of polarization', fontsize=12)
plt.savefig(plotDir+'spectrum_reflected_polarization.pdf', bbox_inches='tight')
