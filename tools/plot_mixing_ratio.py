import matplotlib.pyplot as plt
import numpy as np
import sys, os, random
from scipy import interpolate

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))

molecules = [ 'H2O', 'CO2', 'N2O', 'CO', 'CH4', 'NO', 'NO2', 'NH3', 'HCN', 'C2H2', 'PH3', 'C2H4', 'H2', 'CS', 'Na', 'K' ]
molList = [ r'H$_2$O', r'CO$_2$', r'N$_2$O', 'CO', r'CH$_4$', 'NO', r'NO$_2$', r'NH$_3$', 'HCN', r'C$_2$H$_2$', r'PH$_3$', r'C$_2$H$_4$', r'H$_2$', 'CS', 'Na', 'K' ]

numMol = len(molList)

colors = [ 'black', 'gray', 'maroon', 'red', 'olive', 'yellow', 'green', 'lime', 'teal', 'aqua', 'navy', 'blue', 'purple', 'fuchsia', 'coral', 'sienna' ]

pressure, temperature = np.loadtxt(directory[:-6]+'input/'+atmosphere+'/pressure_temperature.dat', unpack=True)

fileNumber, Pgrid, Tgrid = np.loadtxt(directory[:-6]+'dat/molecules/PTgrid.dat', unpack=True)

data = np.loadtxt(directory[:-6]+'dat/molecules/mixingratios.dat')

VMRmol = np.zeros((numMol,len(Pgrid)))

for i in range(numMol):
    VMRmol[i,:] = data[:,i+2]

# Change zero mixing ratios to 1e-220 because of log10
for i in range(numMol):
    for j in range(len(Pgrid)):
        if VMRmol[i,j] == 0.:
            VMRmol[i,j] = 1e-200

VMR = np.zeros((numMol,len(pressure)))
for i in range(numMol):
    
    mT, mP = np.meshgrid(temperature, np.log10(pressure))
    mT = np.transpose(mT)
    mP = np.transpose(mP)
    
    VMR[i,:] = 10.**np.diag(interpolate.griddata((Tgrid, np.log10(Pgrid)), np.log10(VMRmol[i,:]), (mT, mP), method='linear'))

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

plt.xlabel('Volume mixing ratio')
plt.ylabel('Pressure [bar]')

plt.plot(VMR[0,:], pressure, '-', color=colors[0], label=molList[0])
plt.plot(VMR[1,:], pressure, '-', color=colors[1], label=molList[1])
plt.plot(VMR[2,:], pressure, '-', color=colors[2], label=molList[2])
plt.plot(VMR[3,:], pressure, '-', color=colors[3], label=molList[3])
plt.plot(VMR[4,:], pressure, '-', color=colors[4], label=molList[4])
plt.plot(VMR[5,:], pressure, '-', color=colors[5], label=molList[5])
plt.plot(VMR[6,:], pressure, '-', color=colors[6], label=molList[6])
plt.plot(VMR[7,:], pressure, '-', color=colors[7], label=molList[7])
plt.plot(VMR[8,:], pressure, '-', color=colors[8], label=molList[8])
plt.plot(VMR[9,:], pressure, '-', color=colors[9], label=molList[9])
plt.plot(VMR[10,:], pressure, '-', color=colors[10], label=molList[10])
plt.plot(VMR[11,:], pressure, '-', color=colors[11], label=molList[11])
plt.plot(VMR[12,:], pressure, '-', color=colors[12], label=molList[12])
plt.plot(VMR[13,:], pressure, '-', color=colors[13], label=molList[13])
plt.plot(VMR[14,:], pressure, '-', color=colors[14], label=molList[14])
plt.plot(VMR[15,:], pressure, '-', color=colors[15], label=molList[15])

plt.gca().invert_yaxis()

plt.yscale('log')
plt.xscale('log')

plt.xlim(5.e-11,1.e-2)
plt.ylim(max(pressure),min(pressure))

plt.legend(loc='center right')

plt.savefig(directory[:-6]+'/input/'+atmosphere+'/plot/mixing_ratios.pdf', bbox_inches='tight')
