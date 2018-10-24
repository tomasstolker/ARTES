import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")

# ------------------------------------------------------------

stokesFile = output+'output/phase.dat'
normFile = output+'output/normalization.dat'
plotDir = output+'plot/'

phase, stokesI, errorI, stokesQ, errorQ, stokesU, errorU, stokesV, errorV = np.loadtxt(stokesFile, unpack=True)
wavelength, norm1, norm2 = np.loadtxt(normFile, unpack=True)

stokesPI = np.sqrt(stokesQ**2+stokesU**2+stokesV**2)

# Lambertian surface
ss_albedo = 1.0
x = np.linspace(0,180,1e3)
y = (2./3.)*ss_albedo*(np.sin(x*np.pi/180.)+(np.pi-(x*np.pi/180.))*np.cos(x*np.pi/180.))/np.pi
plt.plot(x,y, ls='--')
ss_albedo = 0.5
x = np.linspace(0,180,1e3)
y = (2./3.)*ss_albedo*(np.sin(x*np.pi/180.)+(np.pi-(x*np.pi/180.))*np.cos(x*np.pi/180.))/np.pi
plt.plot(x,y, ls='--')

# Stokes I
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('Normalized Stokes I', fontsize=12)
plt.plot(phase, stokesI/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(0,1)
plt.savefig(plotDir+'phase_I.pdf', bbox_inches='tight')
plt.clf()

# Stokes Q
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('Normalized Stokes Q', fontsize=12)
plt.plot(phase, stokesQ/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Q.pdf', bbox_inches='tight')
plt.clf()

# Stokes U
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('Normalized Stokes U', fontsize=12)
plt.plot(phase, stokesU/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_U.pdf', bbox_inches='tight')
plt.clf()

# Stokes V
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('Normalized Stokes V', fontsize=12)
plt.plot(phase, stokesV/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_V.pdf', bbox_inches='tight')
plt.clf()

# -Q/I
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('-Q/I', fontsize=12)
plt.plot(phase[stokesI>0.], -stokesQ[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Qpol.pdf', bbox_inches='tight')
plt.clf()

# U/I
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('U/I', fontsize=12)
plt.plot(phase[stokesI>0.], stokesU[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Upol.pdf', bbox_inches='tight')
plt.clf()

# V/I
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('V/I', fontsize=12)
plt.plot(phase[stokesI>0.], stokesV[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Vpol.pdf', bbox_inches='tight')
plt.clf()

# Degree of polarization
plt.xlabel('Phase angle [degrees]', fontsize=12)
plt.ylabel('Degree of polarization', fontsize=12)
plt.plot(phase[stokesI>0.], stokesPI[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(0,1)
plt.savefig(plotDir+'phase_polarization.pdf', bbox_inches='tight')