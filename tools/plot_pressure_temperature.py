import matplotlib.pyplot as plt
import numpy as np
import sys, os

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))

pressure, temperature = np.loadtxt(directory[:-6]+'input/'+atmosphere+'/pressure_temperature.dat', unpack=True)

plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [bar]')

plt.gca().invert_yaxis()
plt.yscale('log')

plt.ylim(max(pressure), min(pressure))

plt.plot(temperature, pressure, ls='-')

plt.savefig(directory[:-6]+'input/'+atmosphere+'/plot/pressure_temperature.pdf', bbox_inches='tight')