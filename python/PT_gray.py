import numpy as np
import math, os

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

Teff = 800. # [K]
kappa = 1e-2 # [cm2 g-1]
log_g = 4. # [cm2 s-2]
Pmin = 1e-3 # [bar]
Pmax = 1e2 # [bar]
Plevels = 40

# ------------------------------------------------------------

g = 10.**log_g

P = np.logspace(np.log10(Pmin),np.log10(Pmax),Plevels) # [bar]

P *= 1e6 # [bar] -> [Ba]

T = np.zeros(np.size(P))
tau = np.zeros(np.size(P))

for i in range(np.size(P)):

    tau[i] = kappa*P[i]/g
    T[i] = ( (3.*(Teff**4)/4.) * ( (2./3.) + tau[i] ) )**(1./4.)

P *= 1e-6 # [Ba] -> [bar]

z = np.array(zip(P,T))

directory = os.path.dirname(os.path.abspath(__file__))

f = open(directory[:-6]+'input/'+atmosphere+'/pressure_temperature.dat', 'w')
f.write("# Pressure [bar] - Temperature [K]\n\n")
np.savetxt(f,z)
f.close()