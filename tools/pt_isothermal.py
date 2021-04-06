import os
import numpy as np

# ------------------------------------------------------------
# Input

Tiso = 800.  # (K)
Pmin = 1e-3  # (bar)
Pmax = 1e2  # (bar)
Plevels = 40

# ------------------------------------------------------------

pressure = np.logspace(np.log10(Pmin), np.log10(Pmax), Plevels)  # (bar)
temperature = np.full(pressure.shape[0], Tiso)

np.savetxt('pressure_temperature.dat',
           np.column_stack((pressure, temperature)),
           header='Pressure (bar) - Temperature (K)')
