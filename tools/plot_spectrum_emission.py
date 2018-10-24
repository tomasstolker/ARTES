import matplotlib.pyplot as plt
import numpy as np
import sys, os

# ------------------------------------------------------------
# Input

output = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))

data = np.loadtxt(directory[:-6]+'output/'+output+'/output/spectrum.dat')
wavelength = data[:,0]
flux = data[:,1]

plt.plot(wavelength, flux, ls='-')

plt.xlim(np.amin(wavelength),np.amax(wavelength))
plt.ylim(1.e-30,1.e-11)

plt.xscale('log')
plt.yscale('log')

plt.xlabel('Wavelength [micron]')
plt.ylabel('Flux [W/m$^2$/micron]')

plt.savefig(directory[:-6]+'output/'+output+'/plot/spectrum_emission.pdf', bbox_inches='tight')