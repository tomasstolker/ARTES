import sys
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Input

pt_profile = sys.argv[1]

# ------------------------------------------------------------

pressure, temperature = np.loadtxt('pressure_temperature.dat', unpack=True)

plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [bar]')

plt.gca().invert_yaxis()
plt.yscale('log')

plt.plot(temperature, pressure, ls='-')

plt.savefig('pressure_temperature.pdf', bbox_inches='tight')
