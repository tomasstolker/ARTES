import os
import sys

import numpy as np
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")

# ------------------------------------------------------------

plot_dir = output+'plot/'
data_dir = output+'output/'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

wavelength, stokes_i, stokes_q, stokes_u, stokes_v = np.loadtxt(data_dir+'spectrum.dat', unpack=True)
wavelength, norm_1, norm_2 = np.loadtxt(data_dir+'normalization.dat', unpack=True)

plt.plot(wavelength, stokes_i, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength (um)', fontsize=12)
plt.ylabel('Flux (W m$^{-2}$ um$^{-1}$)', fontsize=12)
plt.savefig(plot_dir+'spectrum_reflected_flux.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, stokes_i/norm_1, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength (um)', fontsize=12)
plt.ylabel('Normalized Stokes I', fontsize=12)
plt.savefig(plot_dir+'spectrum_reflected_flux_1.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, stokes_i/norm_2, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength (um)', fontsize=12)
plt.ylabel('Normalized Stokes I', fontsize=12)
plt.savefig(plot_dir+'spectrum_reflected_flux_2.pdf', bbox_inches='tight')
plt.clf()

plt.plot(wavelength, np.sqrt(stokes_q**2+stokes_u**2+stokes_v**2)/stokes_i, ls='-')
plt.xlim(min(wavelength),max(wavelength))
plt.xlabel('Wavelength (um)', fontsize=12)
plt.ylabel('Degree of polarization', fontsize=12)
plt.savefig(plot_dir+'spectrum_reflected_polarization.pdf', bbox_inches='tight')
