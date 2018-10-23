import numpy as np
import os

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

Tiso = 800. # [K]
Pmin = 1e-3 # [bar]
Pmax = 1e2 # [bar]
Plevels = 40

# ------------------------------------------------------------

P = np.logspace(np.log10(Pmin),np.log10(Pmax),Plevels) # [bar]
P *= 1e6 # [bar] -> [Ba]
T = np.zeros(np.size(P))

for i in range(np.size(P)):
    T[i] = Tiso

P *= 1e-6 # [Ba] -> [bar]

z = np.array(zip(P,T))

directory = os.path.dirname(os.path.abspath(__file__))

f = open(directory[:-6]+'input/'+atmosphere+'/pressure_temperature.dat', 'w')
f.write("# Pressure [bar] - Temperature [K]\n\n")
np.savetxt(f,z)
f.close()