import matplotlib.pyplot as plt
import numpy as np
import os, sys
from matplotlib.colors import LogNorm
from astropy.io import fits

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]
wavelength = int(sys.argv[2])
azimuth = int(sys.argv[3])

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))

hdulist = fits.open(directory[:-6]+'input/'+atmosphere+'/'+'atmosphere.fits')
r_face = hdulist[0].data # [m]
t_face = hdulist[1].data # [deg]
scattering = hdulist[6].data
absorption = hdulist[7].data
hdulist.close()

scattering = scattering[wavelength,azimuth,:,:]
absorption = absorption[wavelength,azimuth,:,:]

tau_scat = np.rot90(scattering, k=3)
tau_abs = np.rot90(absorption, k=3)

tau_scat_old = 0.
tau_abs_old = 0.

for i in range(np.size(tau_scat,0)):

    tau_scat[np.size(tau_scat,0)-i-1,:] = tau_scat_old + tau_scat[np.size(tau_scat,0)-i-1,:] * (r_face[np.size(r_face)-i-1]-r_face[np.size(r_face)-i-2])

    tau_abs[np.size(tau_abs,0)-i-1,:] = tau_abs_old + tau_abs[np.size(tau_abs,0)-i-1,:] * (r_face[np.size(r_face)-i-1]-r_face[np.size(r_face)-i-2])

    tau_scat_old = tau_scat[np.size(tau_scat,0)-i-1,:]
    tau_abs_old = tau_abs[np.size(tau_abs,0)-i-1,:]
    
r_face *= 1e-3 # [km]
r_face -= r_face[0]

r = np.zeros((np.size(r_face)-1))
t = np.zeros((np.size(t_face)-1))

for i in range(np.size(r_face)-1):
    r[i] = (r_face[i]+r_face[i+1])/2. - r_face[0]

for i in range(np.size(t_face)-1):
    t[i] = (t_face[i]+t_face[i+1])/2.

plt.figure(figsize=(8,5))

fig = plt.imshow(tau_scat, origin='lower', cmap='inferno', interpolation='nearest', aspect='auto', extent=[t_face[0],t_face[len(t_face)-1],r_face[0],r_face[len(r_face)-1]], norm=LogNorm(vmin=np.amin(tau_scat), vmax=np.amax(tau_scat)))

cb = plt.colorbar(fig)
cb.set_label('Scattering optical depth')

plt.xlabel('Latitude [deg]')
plt.ylabel('Height [km]')

plt.savefig(directory[:-6]+'input/'+atmosphere+'/'+'plot/tau_scattering.pdf', bbox_inches='tight')
plt.clf()

plt.figure(figsize=(8,5))

fig = plt.imshow(tau_abs, origin='lower', cmap='inferno', interpolation='nearest', aspect='auto', extent=[t_face[0],t_face[len(t_face)-1],r_face[0],r_face[len(r_face)-1]])

cb = plt.colorbar(fig)
cb.set_label('Absorption optical depth')

plt.xlabel('Latitude [deg]')
plt.ylabel('Height [km]')

plt.savefig(directory[:-6]+'input/'+atmosphere+'/'+'plot/tau_absorption.pdf', bbox_inches='tight')