import os
import math
import numpy as np

# ------------------------------------------------------------
# Input

Teff = 800. # [K]
kappa = 1e-2 # [cm2 g-1]
log_g = 4. # [cm2 s-2]
Pmin = 1e-3 # [bar]
Pmax = 1e2 # [bar]
Plevels = 40

# ------------------------------------------------------------

g = 10.**log_g

pressure = np.logspace(np.log10(Pmin), np.log10(Pmax), Plevels) # [bar]
pressure *= 1e6 # [bar] -> [Ba]

temperature = np.zeros(pressure.shape)
tau = np.zeros(pressure.shape)

for i, item in enumerate(pressure):
    tau[i] = kappa*item/g
    temperature[i] = ( (3.*(Teff**4)/4.) * ( (2./3.) + tau[i] ) )**(1./4.)

pressure *= 1e-6 # [Ba] -> [bar]

z = np.array(zip(pressure, temperature))

f = open('pressure_temperature.dat', 'w')
f.write("# Pressure [bar] - Temperature [K]\n\n")
np.savetxt(f, z)
f.close()
