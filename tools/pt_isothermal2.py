import os
import numpy as np

# ------------------------------------------------------------
# Input

Tiso = 800. # [K]
Pmin = 1e-3 # [bar]
Pmax = 1e2 # [bar]
Plevels = 40

# ------------------------------------------------------------

pressure = np.logspace(np.log10(Pmin), np.log10(Pmax), Plevels) # [bar]
pressure *= 1e6 # [bar] -> [Ba]

temperature = np.zeros(pressure.shape)
temperature[:] = Tiso

pressure *= 1e-6 # [Ba] -> [bar]

z = np.array(zip(pressure, temperature))

f = open('pressure_temperature.dat', 'w')
f.write("# Pressure [bar] - Temperature [K]\n\n")
np.savetxt(f, z)
f.close()
