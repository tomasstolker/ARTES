import matplotlib.pyplot as plt
import numpy as np
import sys, os

# ------------------------------------------------------------
# Input

output = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))
plotDir = directory[:-6]+'output/'+output+'/plot/'
dataDir = directory[:-6]+'output/'+output+'/output/'

wavelength, stokesi, stokesq, stokesu, stokesv = np.loadtxt(dataDir+'spectrum.dat', unpack=True)
wavelength, norm1, norm2 = np.loadtxt(dataDir+'normalization.dat', unpack=True)

plt.plot(wavelength, stokesi, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]')
plt.ylabel('Flux [W m$^{-2}$ micron$^{-1}$]')
plt.savefig(plotDir+'spectrum_reflected_flux.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, stokesi/norm1, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]')
plt.ylabel('Normalized Stokes I')
plt.savefig(plotDir+'spectrum_reflected_flux_1.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, stokesi/norm2, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]')
plt.ylabel('Normalized Stokes I')
plt.savefig(plotDir+'spectrum_reflected_flux_2.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, np.sqrt(stokesq**2+stokesu**2+stokesv**2)/stokesi, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength [micron]')
plt.ylabel('Degree of polarization')
plt.savefig(plotDir+'spectrum_reflected_polarization.pdf', bbox_inches='tight')