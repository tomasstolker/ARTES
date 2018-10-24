import matplotlib.pyplot as plt
import numpy as np
import sys, os
from astropy.io import fits

# ------------------------------------------------------------
# Input

output = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))

stokesFile = directory[:-6]+'output/'+output+'/output/phase.dat'
normFile = directory[:-6]+'output/'+output+'/output/normalization.dat'
plotDir = directory[:-6]+'output/'+output+'/plot/'

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
plt.xlabel('Phase angle [degrees]')
plt.ylabel('Normalized Stokes I')
plt.plot(phase, stokesI/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(0,1)
plt.savefig(plotDir+'phase_I.pdf', bbox_inches='tight')
plt.clf()

# Stokes Q
plt.xlabel('Phase angle [degrees]')
plt.ylabel('Normalized Stokes Q')
plt.plot(phase, stokesQ/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Q.pdf', bbox_inches='tight')
plt.clf()

# Stokes U
plt.xlabel('Phase angle [degrees]')
plt.ylabel('Normalized Stokes U')
plt.plot(phase, stokesU/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_U.pdf', bbox_inches='tight')
plt.clf()

# Stokes V
plt.xlabel('Phase angle [degrees]')
plt.ylabel('Normalized Stokes V')
plt.plot(phase, stokesV/norm2, ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_V.pdf', bbox_inches='tight')
plt.clf()

# -Q/I
plt.xlabel('Phase angle [degrees]')
plt.ylabel('-Q/I')
plt.plot(phase[stokesI>0.], -stokesQ[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Qpol.pdf', bbox_inches='tight')
plt.clf()

# U/I
plt.xlabel('Phase angle [degrees]')
plt.ylabel('U/I')
plt.plot(phase[stokesI>0.], stokesU[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Upol.pdf', bbox_inches='tight')
plt.clf()

# V/I
plt.xlabel('Phase angle [degrees]')
plt.ylabel('V/I')
plt.plot(phase[stokesI>0.], stokesV[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plotDir+'phase_Vpol.pdf', bbox_inches='tight')
plt.clf()

# Degree of polarization
plt.xlabel('Phase angle [degrees]')
plt.ylabel('Degree of polarization')
plt.plot(phase[stokesI>0.], stokesPI[stokesI>0.]/stokesI[stokesI>0.], ls='-')
plt.xlim(0,180)
plt.ylim(0,1)
plt.savefig(plotDir+'phase_polarization.pdf', bbox_inches='tight')