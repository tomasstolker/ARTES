import matplotlib.pyplot as plt
import numpy as np
import sys, os

# ------------------------------------------------------------
# Input

output = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))

wavelength, extinction, absorption, scattering = np.loadtxt(directory[:-6]+'output/'+output+'/output/optical_depth.dat', unpack=True)

plt.plot(wavelength, absorption, '-', label='Absorption')
plt.plot(wavelength, scattering, '-', label='Scattering')

plt.xlabel('Wavelength [micron]')
plt.ylabel('Optical depth')

plt.xscale('log')
plt.yscale('log')

plt.xlim(min(wavelength),max(wavelength))

plt.legend(loc='upper left')

plt.savefig(directory[:-6]+'output/'+output+'/plot/optical_depth.pdf', bbox_inches='tight')
